"""Helper functions for loading and transforming data"""

import os

# isort: off
import xarray as xr
import pint_xarray
import cf_xarray.units
# isort: on

import pandas as pd

# Directories
PROCESSED_DIR = "data/processed"  # processed data files
RESULTS_DIR = "results"  # results
SIMS_DIR = "simulations"  # simulation outputs


def apply_scales(
    dataset: xr.Dataset,
    scales: xr.DataArray,
    label: str | None = None,
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

    # Infer observation type, units, and label from filename
    file_label, obs_type = filename.removeprefix("obs_").removesuffix(".csv").split("_")
    if obs_type not in ["concentration", "deposition"]:
        raise ValueError(f"Could not infer observation type from '{filename}'")
    units = {"concentration": "ug/m3", "deposition": "t/km2/yr"}
    units = units[obs_type]
    if not label:
        label = file_label

    # Load and rename columns
    path = os.path.join(data_dir, filename)
    obs = pd.read_csv(path).rename(columns={units: obs_type})

    # Exclude observations with no microplastics
    obs = obs.loc[obs[obs_type].gt(0)].reset_index(drop=True)

    # Get size range units
    size_units = obs["size_units"].dropna().unique()
    if len(size_units) > 1:
        raise ValueError("Column 'size_units' contains multiple units: {size_units}")
    size_units = size_units[0]

    # Convert to xr.Dataset
    # fmt: off
    cols = [
        "study", "doi", "lat", "lon", obs_type, "size_low", "size_high", "shape",
        "setting"
    ]
    # fmt: on
    obs = (
        obs[cols]
        .to_xarray()
        .set_coords(["lat", "lon"])
        .assign_attrs(description=f"Observed atmospheric microplastic {obs_type}")
    )

    # Set attributes
    obs[obs_type].attrs.update(long_name=f"Mass {obs_type}", units=units)
    obs["lat"].attrs.update(standard_name="latitude", units="degrees_north", axis="Y")
    obs["lon"].attrs.update(standard_name="longitude", units="degrees_east", axis="X")
    obs["size_low"].attrs.update(long_name="Aerodynamic diameter", units=size_units)
    obs["size_high"].attrs.update(long_name="Aerodynamic diameter", units=size_units)
    obs.attrs["label"] = label

    return obs


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
