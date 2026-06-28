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
from numpy.typing import ArrayLike
from pandas import DataFrame
from pandas.io.formats.style import Styler
from pint import Quantity
from pint.errors import DimensionalityError
from pint_xarray import unit_registry as ureg
from scipy import optimize

# Molecular weight of dry air
# https://github.com/geoschem/geos-chem/blob/14.1.1/Headers/physconstants.F90#L26
AIR_MW = 28.9644 * ureg("g/mol")

# Density of microplastic aerosols
PLASTIC_DENSITY = 1 * ureg("g/cm3")

# Molecular weight of microplastics
PLASTIC_MW = 29 * ureg("g/mol")

# Simulation tracer sizes (aerodynamic diameter [µm])
PLASTIC_SIZES = [0.3, 2.5, 7, 15, 35, 70]

# Shape factors for atmospheric microplastics
# Multiplier to compute volume from cube of particle size
PLASTIC_SHAPE_FACTORS = {
    # Ellipse volume = pi / 6 * length * width * height
    "fragments": np.pi / 6 * 0.68 * (0.4 * 0.68),  # l = 1; h = 0.68 l; w = 0.4 h
}

# Simulated atmospheric microplastic sources
PLASTIC_SOURCES = ["ocean", "mmpw", "agricultural", "residential", "road"]

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


def powerlaw_compute_bin_edges(
    xmin: float, vol_mean_sizes: list[float], alpha: float, figs: int | None = None
) -> list[float]:
    """
    Compute size bin edges assuming a power law size distribution

    Args:
        xmin: Minimum size (lower edge of first bin)
        vol_mean_sizes: Volume mean size of each bin
        alpha: Power law parameter alpha, the exponent (must be > 1)
        figs: Round bin edges to figs signficant figures (default: no rounding)

    Both xmin and alpha must have size 1.

    Returns:
        Edges of the size bins

    Derivation:
                 n(x) = C * x^-alpha
        N(xmin, xmax) = Integral[C * x^-alpha] dx
                      = C / (1 - alpha) * (xmax^(1 - alpha) - xmin^(1 - alpha))
                 v(x) = k * x^3
        V(xmin, xmax) = Integral[v(x) * n(x)] dx
                      = Integral[k * C * x^(3 - alpha)] dx
                      = k * C / (4 - alpha) * (xmax^(4 - alpha) - xmin^(4 - alpha))
               v(x_v) = V(xmin, xmax) / N(xmin, xmax)
            k * x_v^3 = k * (1 - alpha) / (4 - alpha)
                          * (xmax^(4 - alpha) - xmin^(4 - alpha))
                          / (xmax^(1 - alpha) - xmin^(1 - alpha))
                  x_v = [(1 - alpha) / (4 - alpha)
                          * (xmax^(4 - alpha) - xmin^(4 - alpha))
                          / (xmax^(1 - alpha) - xmin^(1 - alpha))]^(1/3)
    """

    alpha = float(alpha)  # numpy disallows negative integer powers
    if alpha <= 1:
        raise ValueError(f"alpha must be > 1; got {alpha}")
    if alpha >= 4:
        raise ValueError(f"alpha must be < 4; got {alpha}")

    def compute_upper_edge(x1: float) -> float:
        vol_mean = exp1 / exp4 * (x1**exp4 - x0**exp4) / (x1**exp1 - x0**exp1)
        return vol_mean - vol_mean_sizes[i] ** 3

    # Solve for bin upper edges using numerical root-finding
    edges = [xmin]
    exp1 = 1 - alpha
    exp4 = 4 - alpha
    for i in range(len(vol_mean_sizes)):
        x0 = edges[-1]  # lower edge of bin
        sol = optimize.root(compute_upper_edge, x0 + 1)  # initial guess = x0 + 1
        x1 = sol.x
        if figs:
            # round to figs significant figures
            mags = 10 ** (figs - 1 - np.floor(np.log10(x1)))
            x1 = np.round(x1 * mags) / mags
        edges = edges + x1.tolist()

    return edges


def powerlaw_compute_c(
    n: Quantity, xmin: Quantity, xmax: Quantity, alpha: float | ArrayLike
) -> Quantity:
    """Compute power law parameter C

    Power law number size distribution:
        n(x) = C * x^-alpha

    Args:
        n: Number of particles in size range
        xmin, xmax: Size range bounds
        alpha: Power law parameter alpha, the exponent (must be > 1)

    All arguments must have size 1 or n.

    Returns:
        Value(s) of C

    Derivation:
                 n(x) = C * x^-alpha
        N(xmin, xmax) = Integral[C * x^-alpha] dx
                      = C / (1 - alpha) * [xmax^(1 - alpha) - xmin^(1 - alpha)]
                    C = N * (1 - alpha) / [xmax^(1 - alpha) - xmin^(1 - alpha)]
    """

    alpha = np.asarray(alpha).astype(float)  # numpy disallows negative integer powers
    if np.any(alpha <= 1):
        raise ValueError(f"alpha must be > 1; got {alpha}")

    exp = 1 - alpha
    try:
        x_diff = xmax**exp - xmin**exp
    except DimensionalityError as exc:
        raise ValueError(
            "Cannot raise a Quantity to an array exponent. Use a single value for alpha"
            "or remove the units from xmin and xmax.",
        ) from exc
    scale = n * exp / x_diff

    return scale


def powerlaw_compute_mass(
    xmin: Quantity,
    xmax: Quantity,
    alpha: float | ArrayLike,
    c: Quantity,
    shape_factor: float | ArrayLike,
    density: Quantity,
) -> Quantity:
    """
    Compute mass in size range assuming power law particle size distribution

    Power law mass size distribution:
        M(x) = rho * k * C * x^(3 - alpha)

    Args:
        xmin, xmax: Size range bounds
        alpha: Power law parameter alpha, the exponent (must be > 1)
        c: Power law parameter C
        shape_factor: Multiplier to compute volume from cube of particle size
        density: Particle density

    All arguments must have size 1 or n.

    Returns:
        Mass of particles in size range(s) [xmin, xmax]

    Derivation (rho = density; k = shape factor):
                 m(x) = rho * k * x^3
                 n(x) = C * x^-alpha
        M(xmin, xmax) = Integral[m(x) * n(x)] dx
                      = Integral[rho * k * C * x^(3 - alpha)] dx
                      = rho * k * C / (4 - alpha) * [xmax^(4 - alpha) - xmin^(4 - alpha)]
    """

    alpha = np.asarray(alpha).astype(float)  # numpy disallows negative integer powers
    if np.any(alpha <= 1):
        raise ValueError(f"alpha must be > 1; got {alpha}")

    shape_factor = np.asarray(shape_factor)

    # Ensure output units don't have particle in numerator
    c = c * ureg("1/particle")

    exp = 4 - alpha
    try:
        x_diff = xmax**exp - xmin**exp
    except DimensionalityError as exc:
        raise ValueError(
            "Cannot raise a Quantity to an array exponent. Use a single value for alpha"
            "or remove the units from xmin and xmax.",
        ) from exc
    mass = density * shape_factor * c / exp * x_diff

    return mass


def powerlaw_compute_number(
    xmin: Quantity, xmax: Quantity, alpha: float | ArrayLike, c: Quantity
) -> Quantity:
    """
    Compute number of particles in size range assuming power law size distribution

    Power law number size distribution:
        n(x) = C * x^-alpha

    Args:
        xmin, xmax: Size range bounds
        alpha: Power law parameter alpha, the exponent (must be > 1)
        c: Power law parameter C, the scaling factor

    All arguments must have size 1 or n.

    Returns:
        Number of particles in size range(s) [xmin, xmax]

    Derivation:
                 n(x) = C * x^-alpha
        N(xmin, xmax) = Integral[C * x^-alpha] dx
                      = C / (1 - alpha) * [xmax^(1 - alpha) - xmin^(1 - alpha)]
    """

    alpha = np.asarray(alpha).astype(float)  # numpy disallows negative integer powers
    if np.any(alpha <= 1):
        raise ValueError(f"alpha must be > 1; got {alpha}")

    exp = 1 - alpha
    try:
        x_diff = xmax**exp - xmin**exp
    except DimensionalityError as exc:
        raise ValueError(
            "Cannot raise a Quantity to an array exponent. Use a single value for alpha"
            "or remove the units from xmin and xmax.",
        ) from exc
    number = c / exp * x_diff

    return number
