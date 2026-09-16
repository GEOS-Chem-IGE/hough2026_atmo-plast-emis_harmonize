"""General helper functions"""

import os
import re

# isort: off
import xarray as xr
import pint_xarray
import cf_xarray.units
# isort: on

import numpy as np
from IPython.display import HTML
from matplotlib import colors
from pandas import DataFrame
from pandas.io.formats.style import Styler
from pint import Quantity
from pint_xarray import unit_registry as ureg

from utils.constants import AIR_MW, PLASTIC_DENSITY, PLASTIC_SHAPE_FACTORS

# Colormap for plots
WhBlYlRd = colors.ListedColormap(
    np.genfromtxt(os.path.join(os.path.dirname(__file__), "WhBlYlRd.txt")) / 255,
    name="WhBlYlRd",
)


def clean_styler(x: DataFrame | Styler) -> HTML:
    """Format a DataFrame as an HTML table without element IDs

    The HTML table returned by DataFrame.style includes IDs for table elements which
    change each time the table is generated, even if the table content is unchanged.
    This function removes the IDs from the table to ensure the output only changes if
    the table content changes.
    """

    # Get styled HTML string
    if isinstance(x, DataFrame):
        x = x.style
    html = x.to_html()

    # Remove cell IDs
    html = re.sub(r" id=\"\w+\"", "", html)

    # Remove cell classes
    html = re.sub(r" class=\".+?\"", "", html)

    # Remove whitespace
    html = re.sub(r"&nbsp;<", "<", html)
    html = re.sub(r" +>", ">", html)

    return HTML(html)


def convert_molec_to_mass(spec: xr.DataArray, spec_mw: Quantity) -> xr.DataArray:
    """Convert molecules to mass e.g. molec/cm3 -> µg/m3

    Args:
        spec: Species with units including molecules e.g. molec/cm3, molec/cm2/s
        spec_mw: Species molecular weight (g/mol)

    Returns:
        xr.DataArray of species with quantity expressed as mass e.g. µg/m3, µg/m2/s

    Conversion formula is:

        molecules species    1 mol               µg species
        ------------------ * ----------------- * -----------
        [other units]        6.02e23 molecules   mol species

        = quantity / Avogadro's number * species molecular weight

        = µg species / [other units]
    """

    # Quantify units
    spec = spec.pint.quantify()
    spec_mw = spec_mw.to("µg/mol")

    # Check units
    orig_units = spec.pint.units
    if orig_units._units["particle"] != 1:
        raise ValueError(f"spec_molec must have units 'particle / X'; got '{orig_units}'")

    # Convert
    avogadro = ureg("particle") / ureg("particle").to("mol")
    mass = spec * (1 / avogadro) * spec_mw

    return mass


def convert_vvdry_to_ugm3(
    spec: xr.DataArray, spec_mw: Quantity, air_density: xr.DataArray
) -> xr.DataArray:
    """Convert molar fraction (mol/mol dry air) to (µg/m3)

    Args:
        spec: Species with units of (mol/mol dry air)
        spec_mw: Species molecular weight (g/mol)
        air_density: Dry air density (kg/m3). Must be broadcastable with spec.

    Returns:
        xr.DataArray of species concentration (µg/m3)

    Based on gcpy's convert_to_ugm3 routine:
    https://github.com/geoschem/gcpy/blob/1.6.2/gcpy/plot/compare_single_level.py#L443-L446

    The conversion formula is:

        mol species   g/mol species   µg dry air
        ----------- * ------------- * ----------
        mol dry air   g/mol dry air   m3

        = volume mixing ratio * (MW species / MW dry air) * dry air density

        = µg species / m3
    """

    # Quantify units
    # * pint doesn't undersand "mol-1 dry air" so use "mol mol-1"
    # * pint doesn't know that dimension "lev" with unit "level" is unitless
    spec = spec.assign_attrs(
        units=spec.attrs["units"].replace("mol mol-1 dry", "mol mol-1")
    ).pint.quantify({"lev": None})
    air_density = air_density.pint.quantify({"lev": None}).pint.to("µg/m3")
    spec_mw = spec_mw.to("g/mol")

    # Convert vv to µg/m3
    spec_conc = spec * (spec_mw / AIR_MW) * air_density

    return spec_conc


def format_units_mpl(x: str) -> str:
    """Format units for display on matplotlib plots"""

    return re.sub(r"(-?\d+)", r"$^{\1}$", x)


def mass_to_particles(
    mass: xr.DataArray,
    shape_factor: float = PLASTIC_SHAPE_FACTORS["fragments"],
    density: Quantity = PLASTIC_DENSITY,
) -> xr.DataArray:
    """Convert mass to number of paticles

    Args:
        mass: Mass of particles. Must have coordinate "size" that stores volume-weighted
            mean particle size for each size bin.
        shape_factor: Multiplier to compute volume from cube of particle size (default:
            0.0968 i.e. fragments)
        density: Particle density (default: 1 g/cm3)

    Returns:
        Number of particles in each size bin
    """

    # Determine final units
    mass = mass.pint.quantify()
    for unit_str, unit_exp in mass.pint.units._units.items():
        if unit_exp > 0:
            base_unit = ureg(unit_str)
            if base_unit.dimensionality == "[mass]":
                mass_unit = base_unit.units
                break
    else:
        raise ValueError(f"No mass dimension found in unit {mass.pint.units}")
    units = ureg("particle").units * (mass.pint.units / mass_unit)

    volume = shape_factor * xr.DataArray(mass["size"]).pint.quantify() ** 3
    mass_per_particle = volume * density * ureg("1/particle")
    with xr.set_options(keep_attrs=True):
        particles = (mass * (1 / mass_per_particle)).pint.to(units)
        particles = particles.pint.dequantify().rename(mass.name)

    return particles


def particles_to_mass(
    particles: xr.DataArray,
    shape_factor: float = PLASTIC_SHAPE_FACTORS["fragments"],
    density: Quantity = PLASTIC_DENSITY,
) -> xr.DataArray:
    """Convert number of particles to mass

    Args:
        particles: Number of particles. Must have coordinate "size" that stores
            volume-weighted mean particle size for each size bin.
        shape_factor: Multiplier to compute volume from cube of particle size (default:
            0.0968 i.e. fragments)
        density: Particle density (default: 1 g/cm3)

    Returns:
        Mass of particles in each size bin
    """

    # Determine final units
    particles = particles.pint.quantify()
    for unit_str, unit_exp in particles.pint.units._units.items():
        if unit_exp > 0:
            base_unit = ureg(unit_str)
            if base_unit.dimensionality == "[substance]":
                particle_unit = base_unit.units
                break
    else:
        raise ValueError(f"No particles dimension found in unit {particles.pint.units}")
    units = ureg("Gg").units * (particles.pint.units / particle_unit)

    volume = shape_factor * xr.DataArray(particles["size"]).pint.quantify() ** 3
    mass_per_particle = volume * density * ureg("1/particle")
    with xr.set_options(keep_attrs=True):
        mass = (particles * mass_per_particle).pint.to(units)
        mass = mass.pint.dequantify().rename(particles.name)

    return mass
