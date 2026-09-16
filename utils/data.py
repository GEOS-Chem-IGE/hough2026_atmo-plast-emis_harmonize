"""Helper functions for loading and transforming data"""

import os
import re
from glob import glob

# isort: off
import xarray as xr
import pint_xarray
import cf_xarray.units
# isort: on

import pandas as pd

from utils.constants import (
    PLASTIC_DENSITY,
    PLASTIC_MW,
    PLASTIC_SHAPE_FACTORS,
    PLASTIC_SIZES,
    PLASTIC_SOURCES,
)
from utils.powerlaw import powerlaw_compute_bin_edges, powerlaw_compute_mass
from utils.utils import convert_molec_to_mass, convert_vvdry_to_ugm3

# Directories
PROCESSED_DIR = "data/processed"  # processed data files
RESULTS_DIR = "results"  # results
SIMS_DIR = "simulations"  # simulation outputs


def apply_scales(
    dataset: xr.Dataset, scales: xr.DataArray, label: str | None = None
) -> xr.Dataset:
    """Multiply simulation outputs by scaling factors

    Args:
        data: Dataset of simulation outputs
        scales: DataArray of scaling factors

    Returns:
        A new Dataset containing scaled variables
    """

    # Check inputs
    label = label or scales.attrs["label"]
    if not label and not scales.attrs.get("label"):
        raise ValueError("label must be provided or set as attribute of scales")
    label = label or scales.attrs["label"]

    # Scale emissions, concentration, and deposition
    vars_to_scale = [
        x
        for x in [
            "emission",
            "concentration",
            "total_deposition",
            "dry_deposition",
            "wet_loss",
        ]
        if x in dataset.data_vars
    ]
    with xr.set_options(keep_attrs=True):
        scaled = dataset[vars_to_scale] * scales.drop_attrs()
    for data in scaled.data_vars.values():
        long_name = data.attrs["long_name"].removeprefix("Mean ")
        data.attrs.update(long_name=f"Scaled {long_name.lower()}")

    # Copy over any other variables e.g. area
    for k, v in dataset.data_vars.items():
        if k not in scaled:
            scaled[k] = v

    # Copy any attributes
    scaled.attrs = dataset.attrs

    # Describe and label
    _, obs_label = scales.label.split("_")
    scaled.attrs["description"] = " ".join(
        [dataset.description, f"constrained by {obs_label} observations"]
    )
    scaled.attrs["label"] = label

    return scaled


def combine_land_sources(dataset: xr.Dataset) -> xr.Dataset:
    """Sum data from all terrestrial sources into a single 'land' source

    Args:
        dataset: Dataset to process
    """

    source_vars = [x for x, data in dataset.data_vars.items() if "source" in data.dims]
    land = (
        dataset[source_vars]
        .drop_sel(source="ocean")
        .sum(dim="source", keep_attrs=True)
        .expand_dims({"source": ["land"]})
    )
    merged = xr.concat(
        [dataset.sel(source="ocean"), land], dim="source", data_vars="minimal"
    ).transpose(*dataset.dims)

    return merged


def load_observations(
    filename: str, data_dir: str = RESULTS_DIR, label: str | None = None
) -> xr.Dataset:
    """Load processed observations as an xr.Dataset"""

    # Parse label from filename
    file_label = filename.removeprefix("obs_").removesuffix(".csv").replace("_", "-")
    if not label:
        label = file_label

    path = os.path.join(data_dir, filename)
    obs = pd.read_csv(path)

    # Exclude observations with no microplastics
    obs = obs.loc[obs["mass"].gt(0), :].reset_index(drop=True)

    # Check units
    units = {}
    for col in ["number", "mass"]:
        units[col] = {}
        for measure in ["concentration", "deposition"]:
            obs_units = (
                obs.loc[obs["measure"].eq(measure), f"{col}_units"].dropna().unique()
            )
            if len(obs_units) > 1:
                raise ValueError(
                    f"All {measure} observations must have same {col} units; got {obs_units}"
                )
            units[col][measure] = obs_units[0]
    size_units = obs["size_units"].dropna().unique()
    if len(size_units) > 1:
        raise ValueError(f"All observations must have same size units; got {size_units}")
    size_units = size_units[0]

    # Convert to xr.Dataset
    # fmt: off
    cols = [
        "author", "year", "doi", "measure", "lat", "lon", "size_min", "size_max",
        "number", "mass", "shape", "area"
    ]
    # fmt: on
    cols = [x for x in cols if x in obs.columns]
    obs = obs[cols].to_xarray().set_coords(["lat", "lon"])
    for measure in ["concentration", "deposition"]:
        mask = obs["measure"] == measure
        for variable in ["number", "mass"]:
            varname = f"{variable}_{measure}"
            obs[varname] = xr.where(mask, obs[variable], None)
            obs[varname].attrs.update(
                long_name=f"{variable.title()} {measure} of microplastics",
                units=units[variable][measure],
            )
    obs = obs.drop_vars(["measure", "number", "mass"])

    # Set attributes
    obs["lat"].attrs.update(standard_name="latitude", units="degrees_north", axis="Y")
    obs["lon"].attrs.update(standard_name="longitude", units="degrees_east", axis="X")
    obs["size_min"].attrs.update(long_name="Aerodynamic diameter", units=size_units)
    obs["size_max"].attrs.update(long_name="Aerodynamic diameter", units=size_units)
    obs.attrs["description"] = "Atmospheric microplastic observations"
    obs.attrs["label"] = label

    return obs


def load_sim_output(
    sim_dir: str,
    collection: str,
    source: str = "all",
    chunks: str | dict = "auto",
    reduce_ocean: bool = False,
) -> xr.Dataset:
    """Load raw outputs of a GEOS-Chem microplastic simulation

    Args:
        sim_dir: Directory containing the raw simulation outputs
        collection: One of ["emission", "concentration", "total_deposition",
            "dry_deposition", "wet_loss"]
        source: One of ["all", "ocean", "mmpw", "agricultural", "residential", "road"].
            "all" (default) loads diagnostics for all microplastic sources.
        chunks: Optional dictionary specifying chunks passed to xr.open_mfdataset().
            "auto" (default) sets chunks to a multiple of on-disk chunks.
        reduce_ocean: Whether to reduce ocean microplastics by a factor of 1e-9 (default:
            False). Used for simulations using Fu2023's emission factors, which result in
            ocean microplastic emissions that are orders of magnitude larger than
            emissions from other sources.

    Returns:
        xarray.Dataset including microplastic variable(s) of dimension:
            * [size, time, lat, lon] for "emission" and "dry_deposition" collections
            * [size, time, lev, lat, lon] for all other collections
            * [source, size, time, ...] if source was "all"

    Assumes the simulation outputs are organized in subdirs of sim_dir by source
    [1-ocen, 2-mmpw,  3-agri, 4-resi, 5-road] and timeperiod e.g. 2018-01-01_2019-01-01:

        [sim_dir]
        ├── [source]
        │   ├── [timeperiod]
        │   │   └── OutputDir
        │   │       └── *.nc4
        │   └── ...
        └── ...
    """

    # Validate collection
    globs = {
        "emission": "HEMCO_diagnostics.*.nc",
        "concentration": "GEOSChem.SpeciesConc.*.nc4",
        "total_deposition": "GEOSChem.[DW][re][yt][DL]*.nc4",
        "dry_deposition": "GEOSChem.DryDep.*.nc4",
        "wet_loss": "GEOSChem.WetLoss*.nc4",
    }
    collection = collection.lower()
    if collection not in globs:
        raise ValueError(
            f"Unrecognized collection '{collection}'. Must be one of {list(globs)}"
        )

    # If source is None, load and return data for all sources
    source = source.lower()
    if source == "all":
        data = []
        for src in PLASTIC_SOURCES:
            data.append(
                load_sim_output(
                    sim_dir=sim_dir,
                    collection=collection,
                    source=src,
                    chunks=chunks,
                    reduce_ocean=reduce_ocean,
                )
            )
        data = xr.concat(data, dim="source")
        return data

    # Validate source
    if source not in PLASTIC_SOURCES:
        raise ValueError(
            f"Unrecognized source '{source}'. Must be one of {PLASTIC_SOURCES}"
        )

    # List simulation output files
    # e.g. "data/1-ocen/{START DATE}_{END DATE}/OutputDir/HEMCO_diagnostics.*.nc"
    dirnames = {
        x: f"{i + 1}-" + x[0:4].replace("ocea", "ocen")
        for i, x in enumerate(PLASTIC_SOURCES)
    }
    date_glob = "20[12][8901]-[01][0-9]-[0-3][0-9]"
    file_glob = os.path.join(
        sim_dir,
        dirnames[source],
        date_glob + "_" + date_glob,  # simulation period
        "OutputDir",
        globs[collection],  # diagnostic collection
    )
    paths = sorted(glob(file_glob))
    if not any(paths):
        raise ValueError(f"No {source} {collection} data in {sim_dir}")

    # Load simulation output
    data = xr.open_mfdataset(paths, chunks=chunks)

    # Rename emission variables to conform to GEOS-Chem's standard naming scheme
    # e.g. EmisPST1_OCEN -> Emis_PST1
    # This will facilitate matching on and extracting parts of variable names
    if collection == "emission":
        emis_vars = [v for v in data.variables if v.startswith("EmisPST")]
        for var in emis_vars:
            data[var].attrs.update(
                long_name="Emission flux of species " + data[var].attrs["long_name"][0:4]
            )
        data = data.rename_vars(
            {v: re.sub(r"(PST\d)_[A-Z]+", r"_\1", v) for v in emis_vars}
        )

    # Define variables to reduce by a factor of 1e-9 if reduce_ocean is True
    reduce_vars = [
        "Emis_PST",
        "SpeciesConcVV_PST",
        "DryDep_PST",
        "WetLossConv_PST",
        "WetLossLS_PST",
    ]
    reduce_scale = 1e-9

    # Stack microplastic variables along a new "size" dimension
    # * Each microplastic field (e.g. DryDep) becomes one variable
    # * Some collections include multiple microplastic fields e.g. DryDep and DryDepVel
    size_pattern = r"(\w+_PST)(\d+)"
    size_mapping = {f"PST{i + 1}": x for i, x in enumerate(PLASTIC_SIZES)}
    mp_vars = [v for v in data.variables if re.match(size_pattern, v)]
    if any(mp_vars):
        mp_fields = set(re.sub(size_pattern, r"\1", v) for v in mp_vars)
        for field in sorted(mp_fields):
            # Stack the field variables along a new dimension "size"
            field_vars = sorted(v for v in mp_vars if v.startswith(field))
            size_dim = xr.Variable(
                dims="size",
                data=[
                    size_mapping[re.sub(size_pattern, r"PST\2", v)] for v in field_vars
                ],
                attrs={"long_name": "Aerodynamic diameter", "units": "µm"},
            )
            field_data = (
                data[field_vars].to_dataarray(dim="size").assign_coords(size=size_dim)
            )

            # Possibly reduce ocean microplastics
            if reduce_ocean and source == "ocean" and field in reduce_vars:
                print(f"Scaling ocean {field} by {reduce_scale}")
                with xr.set_options(keep_attrs=True):
                    field_data = reduce_scale * field_data

            # Set attributes
            # * Modify long name
            #   e.g. "Concentration of species PST1" -> "Concentration of microplastics"
            field_data.attrs = data[field_vars[0]].attrs
            field_data.attrs["long_name"] = (
                field_data.attrs["long_name"]
                .replace(" species ", " microplastics ")
                .removesuffix(" PST1")
            )

            # Replace the original field variables with the new field variable
            data = data.drop_vars(field_vars).assign({field: field_data})

    # Add a "source" dimension
    data = data.assign(
        source=xr.Variable(
            dims="source", data=[source], attrs={"long_name": "Emission source"}
        )
    ).set_coords("source")

    return data


def process_sim_output(
    sim_dir: str, label: str, reduce_ocean: bool = False
) -> xr.Dataset:
    """Postprocess outputs of GEOS-Chem microplastic simulations

    Args:
        sim_dir: Directory containing the simulation outputs
        label: Label to identify simulation
        reduce_ocean: Whether to reduce ocean microplastics by a factor of 1e-9 (default:
            False). Used for simulations using Fu2023's emission factors, which result in
            ocean microplastic emissions that are about nine orders of magnitude higher
            than emissions from any other source.

    Returns:
        xarray.Dataset with time-averaged atmospheric microplastic and other variables:
            * emission (kg/km2/year) [source, size, lat, lon]
            * concentration (µg/m3) [source, size, lev, lat, lon]
            * dry_deposition (kg/km2/year) [source, size, lat, lon]
            * wet_loss (kg/km2/year) [source, size, lat, lon]
            * total_deposition (kg/km2/year) [source, size, lat, lon]
            * air_volume (m3) [lev, lat, lon]
            * area: (m2) [lat, lon]
    """

    # Define units for final output
    emis_units = "kg/km2/year"
    conc_units = "µg/m3"
    depo_units = "kg/km2/year"

    # Compute MP emissions
    with xr.set_options(keep_attrs=True):
        emis_outputs = load_sim_output(
            sim_dir=sim_dir, collection="emission", reduce_ocean=reduce_ocean
        )
        emis = emis_outputs["Emis_PST"].pint.quantify().pint.to(emis_units)
        emis.attrs["averaging_method"] = "time-averaged"

    # Compute MP concentration
    with xr.set_options(keep_attrs=True):
        concentration = load_sim_output(
            sim_dir=sim_dir, collection="concentration", reduce_ocean=reduce_ocean
        )

        # Convert concentration mixing ratio to mass per volume
        conc = convert_vvdry_to_ugm3(
            spec=concentration["SpeciesConcVV_PST"],
            spec_mw=PLASTIC_MW,
            air_density=concentration["Met_AIRDEN"],
        ).pint.to(conc_units)

    # Compute MP total deposition
    with xr.set_options(keep_attrs=True):
        wet_outputs = load_sim_output(
            sim_dir=sim_dir, collection="wet_loss", reduce_ocean=reduce_ocean
        )

        # Compute total wet loss (convective + large-scale precipitation)
        wet_loss = (
            (wet_outputs["WetLossConv_PST"] + wet_outputs["WetLossLS_PST"])
            .assign_attrs(long_name="Wet loss flux of soluble microplastics")
            .sum(dim="lev")
            .pint.quantify()
        )

        # Compute wet loss per unit area
        area = wet_outputs["AREA"].pint.quantify()
        wet_loss = (wet_loss / area).pint.to(depo_units)

        # Load dry deposition
        dry_outputs = load_sim_output(
            sim_dir=sim_dir, collection="dry_deposition", reduce_ocean=reduce_ocean
        )

        # Convert dry deposition from (molecules/cm2/s) to (kg/km2/year)
        dry_dep = convert_molec_to_mass(
            spec=dry_outputs["DryDep_PST"], spec_mw=PLASTIC_MW
        ).pint.to(depo_units)

        # Compute total deposition (wet + dry)
        total_dep = (wet_loss + dry_dep).assign_attrs(
            long_name="Total deposition flux of microplastics"
        )

    # Compute dry air volume (m3)
    # This will allow computing atmospheric burden (kg) from concentration (µg/m3)
    # Note: we derive dry air volume from dry air mass + dry air density for consistency
    # with GEOS-Chem's internal conversions, which use dry air mass and density rather
    # than grid cell volume (which is not precisely defined). Deriving dry air volume from
    # dry air mass and density gives values close, but not identical, to the StateMet
    # "Met_AIRVOL" variable, which is computed by multiplying grid cell area by height.
    air_mass = concentration["Met_AD"].isel(source=0, drop=True)
    air_density = concentration["Met_AIRDEN"].isel(source=0, drop=True)
    air_volume = (
        air_mass.pint.quantify({"lev": None}) / air_density.pint.quantify({"lev": None})
    ).assign_attrs(long_name="Dry air volume", averaging_method="time-averaged")

    # Create final dataset
    result = xr.Dataset(
        {
            "emission": emis.pint.dequantify(),
            "concentration": conc.pint.dequantify(),
            "dry_deposition": dry_dep.pint.dequantify(),
            "wet_loss": wet_loss.pint.dequantify(),
            "total_deposition": total_dep.pint.dequantify(),
            "air_volume": air_volume.pint.dequantify(),
            "area": area.isel(source=0).pint.dequantify(),
            "ilev": concentration["ilev"],
            "lat_bnds": concentration["lat_bnds"].isel(source=0),
            "lon_bnds": concentration["lon_bnds"].isel(source=0),
        },
    )

    # Compute mean over study period
    result = result.mean(dim="time", keep_attrs=True)

    # Label
    name = "alternate" if label == "alt" else label
    result = result.assign_attrs(
        description=f"Simulated atmospheric microplastics from {name} configuration",
        label=label,
    )

    # Set order of coordinate variables (cosmetic only; dimension order unchanged)
    result = result.assign_coords(
        {
            "source": result["source"],
            "size": result["size"],
            "lev": result["lev"],
            "lat": result["lat"],
            "lon": result["lon"],
            "ilev": result["ilev"],
            "nb": result["nb"],
        }
    )

    return result


def spatial_integrate(
    dataset: xr.Dataset,
    varname: str,
    units: str | None = None,
    sum_dims: list[str] | str | None = None,
) -> xr.DataArray:
    """Multiply a variable by cell area and sum over all grid cells

    Args:
        data: Dataset with variable to integrate
        varname: Name of variable to integrate
        units: Optional units for output
        sum_dims: Optional dimensions to sum over in addition to [lat, lon]
    """

    # Check inputs
    varnames = [
        "emission",
        "concentration",
        "total_deposition",
        "dry_deposition",
        "wet_loss",
    ]
    if varname not in varnames:
        raise ValueError(f"Unrecognized varname: '{varname}'. Must be one of {varnames}")
    if varname == "concentration":
        if "air_volume" not in dataset.data_vars:
            raise ValueError(
                "data must have variable 'air_volume' to compute atmospheric burden"
            )
    elif "area" not in dataset.data_vars:
        raise ValueError(f"data must have variable 'area' to compute global {varname}")

    # Set target units
    if units is None:
        units = "Gg" if varname == "concentration" else "Gg/year"

    # Set dimensions over which to sum
    if sum_dims is None:
        sum_dims = []
    elif isinstance(sum_dims, str):
        sum_dims = [sum_dims]
    sum_dims = sum_dims + ["lat", "lon"]
    if varname == "concentration":
        sum_dims = sum_dims + ["lev"]

    # Multiply by cell area and sum
    if varname == "concentration":
        integrated = (
            (dataset[varname].pint.quantify() * dataset["air_volume"].pint.quantify())
            .sum(dim=sum_dims)
            .rename("burden")
        )
    else:
        integrated = (
            (dataset[varname].pint.quantify() * dataset["area"].pint.quantify())
            .sum(dim=sum_dims)
            .rename(varname)
        )

    # Convert units
    integrated = integrated.pint.to(units).pint.dequantify("cf")

    # Assign dataset label, if any
    label = dataset.attrs.get("label")
    if label:
        integrated.attrs["label"] = label

    return integrated
