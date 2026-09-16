# isort: off
import xarray as xr
import pint_xarray
import cf_xarray.units
# isort: on

import numpy as np
from numpy.typing import ArrayLike
from pint import Quantity
from pint.errors import DimensionalityError
from pint_xarray import unit_registry as ureg
from scipy import optimize


def powerlaw_compute_bin_edges(
    xmax: float, vol_mean_sizes: list[float], alpha: float, figs: int | None = None
) -> list[float]:
    """
    Compute size bin edges assuming a power law size distribution

    Args:
        xmax: Maximum size (upper edge of largest bin)
        vol_mean_sizes: Volume mean size of each bin
        alpha: Power law parameter alpha, the exponent (must be > 1)
        figs: Round bin edges to figs signficant figures (default: no rounding)

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

    def fun(x0: float) -> float:
        vol = exp1 / exp4 * (xmax**exp4 - x0**exp4) / (xmax**exp1 - x0**exp1)
        return xmean**3 - vol

    # Solve for bin upper edges using numerical root-finding
    edges = [xmax]
    exp1 = 1 - alpha
    exp4 = 4 - alpha
    xmin = xmax
    for i, xmean in enumerate(reversed(vol_mean_sizes)):
        xmax = xmin  # upper edge is lower edge of next largest bin
        xmin = optimize.brentq(fun, a=0.01, b=xmean, xtol=1e-6, maxiter=1000)
        if figs:
            # round to figs significant figures
            mags = 10 ** (figs - 1 - np.floor(np.log10(xmin)))
            xmin = np.round(xmin * mags) / mags
            xmin = xmin.item()
        edges = [xmin] + edges

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
