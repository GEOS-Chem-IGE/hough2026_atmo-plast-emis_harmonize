"""Helper functions for constraining simulation outputs"""

import numpy as np
from scipy import optimize
from numpy.typing import NDArray

# isort: off
import xarray as xr
import pint_xarray
import cf_xarray.units
# isort: on


def compute_error(obs: NDArray, sim: NDArray) -> float:
    """Compute an error score from observed and simulated values

    Args:
        obs: Observed values
        sim: Simulated values

    Returns:
        Sum of squared differences between log10-scaled observed and simulated values
    """

    # Clip to avoid log10(0) issues
    error = np.log10(obs.clip(1e-10)) - np.log10(sim.clip(1e-10))
    error = (error**2).sum()

    return error.item()


def compute_penalty(x: NDArray) -> float:
    """Compute a penalty to avoid negative scaling factors

    Args:
        x: Scaling factors

    Returns:
        -1e10 * sum of negative scaling factors
    """

    penalty = -1e10 * np.minimum(x, 0).sum()

    return penalty.item()


def compute_score(
    init: NDArray,
    sim_conc: xr.DataArray,
    sim_depo: xr.DataArray,
    obs_conc: xr.DataArray,
    obs_depo: xr.DataArray,
) -> float:
    """Compute a score for the current scales

    Args:
        init: Current scaling factors as 1D array
        sim_conc: Simulated concentration
        sim_depo: Simulated deposition
        obs_conc: Observed concentration
        obs_depo: Observed deposition

    Returns:
        Error score

    This function:
        1. Multiplies simulation outputs by the scaling factors
        2. Computes errors (log10-scaled squared difference from observations)
        3. Computes penalty to avoid negative scaling factors
        4. Returns score = errors + penalty

    Error is computed separately for concentration and deposition because they have
    different units and so cannot be summed.
    """

    # Scale simulated concentration and deposition and sum
    scales = reshape_scales(init, coords=sim_conc.coords)
    sum_dims = [d for d in sim_conc.dims if d != "index"]
    sim_conc = (scales * sim_conc).sum(dim=sum_dims)
    sim_depo = (scales * sim_depo).sum(dim=sum_dims)

    # Compute score
    # Compute separately for concentration + deposition because units differ
    err_conc = compute_error(obs=obs_conc.values, sim=sim_conc.values)
    err_depo = compute_error(obs=obs_depo.values, sim=sim_depo.values)

    # Add penalty to avoid negative scales
    penalty = compute_penalty(init)
    score = err_conc + err_depo + penalty

    return score


def fit_scales(
    init: xr.DataArray,
    sim_conc: xr.DataArray,
    sim_depo: xr.DataArray,
    obs_conc: xr.DataArray,
    obs_depo: xr.DataArray,
    method: str = "Nelder-Mead",
    bounds: optimize.Bounds | None = None,
    tol: float | None = 1e-6,
    maxiter: int | None = None,
    **kwargs,
) -> xr.DataArray:
    """
    Compute scaling factors that minimize difference between simulated and
    observed atmospheric microplastic concentration and deposition

    Args:
        init: Initial scales; dims must be [source, size] or [source]
        sim_conc: Simulated concentration; dims must be [{init dims}, lat, lon]
        sim_depo: Simulated deposition; dims must be [{init dims}, lat, lon]
        obs_conc: Observed concentration; dims must be [index] and coords [lat, lon]
        obs_depo: Observed deposition; dims must be [index] and coords [lat, lon]
        method: Solver to use (default: 'Nelder-Mead'). See scipy.optimize.minimize.
        bounds: Optional bounds for the scaling factors
        tol: Optional tolerance (default: 1e-6). See scipy.optimize.minimize for details.
            For the Nelder-Mead solver, the tolerance applies to both the scaling factors
            and the score returned by optimize_step().
        maxiter: Optional number of iterations to perform. The Nelder-Mead solver defaults
            to 200 times the number of scaling factors.
        **kwargs: Additional arguments passed to scipy.optimize.minimize()

    Returns:
        scales: xr.DataArray with dims [source] or [source, size]

    If init has dims [source] then this function optimizes the total magnitude of each
    source. If init has dims [source, size] then it optimizes both the magnitude and
    size distribution of each source.
    """

    # Check input dimensions
    def check_dims(obj, expected, obj_name):
        if isinstance(expected, str):
            expected = [expected]
        if set(obj.dims) != set(expected):
            raise ValueError(f"{obj_name} must have dims {expected}; got {obj.dims}")

    scale_dims = ["source", "size"] if "size" in init.dims else ["source"]
    check_dims(init, scale_dims, "init")
    check_dims(sim_conc, scale_dims + ["lat", "lon"], "sim_conc")
    check_dims(sim_depo, scale_dims + ["lat", "lon"], "sim_depo")
    check_dims(obs_conc, ["index"], "obs_conc")
    check_dims(obs_depo, ["index"], "obs_depo")

    # Confirm coordinates are lat -90 to 90, lon -180 to 180
    def check_coords(obj, obj_name):
        for coord in ["lat", "lon"]:
            if coord not in obj.coords:
                raise ValueError(f"{obj_name} must have '{coord}' coordinates")
        if obj["lat"].min() < -90 or obj["lat"].max() > 90:
            raise ValueError(f"{obj_name} lon must be -90 to 90")
        if obj["lon"].min() < -180 or obj["lon"].max() > 180:
            raise ValueError(f"{obj_name} lon must be -180 to 180")

    check_coords(sim_conc, "sim_conc")
    check_coords(sim_depo, "sim_depo")
    check_coords(obs_conc, "obs_conc")
    check_coords(obs_depo, "obs_depo")

    # Raise if any observations have lon > simulation data's max lon
    # Need this b/c "nearest" indexing does not wrap around dateline
    max_lon = sim_conc["lon"].max()
    if obs_conc["lon"].max() > max_lon:
        raise NotImplementedError(
            f"Can't handle observed concentration with lon > {max_lon}"
        )
    max_lon = sim_depo["lon"].max()
    if obs_depo["lon"].max() > max_lon:
        raise NotImplementedError(
            f"Can't handle observed deposition with lon > {max_lon}"
        )

    # Exclude any zero observations
    obs_conc = obs_conc.loc[obs_conc > 0]
    obs_depo = obs_depo.loc[obs_depo > 0]

    # Get simulated data corresponding to observations
    # Use nearest coordinate
    sim_conc = sim_conc.sel(lat=obs_conc["lat"], lon=obs_conc["lon"], method="nearest")
    sim_depo = sim_depo.sel(lat=obs_depo["lat"], lon=obs_depo["lon"], method="nearest")

    # Standardize units
    conc_units = "µg/m3"
    sim_conc = sim_conc.pint.quantify().pint.to(conc_units).pint.dequantify()
    obs_conc = obs_conc.pint.quantify().pint.to(conc_units).pint.dequantify()
    depo_units = "Mg/km2/yr"
    sim_depo = sim_depo.pint.quantify().pint.to(depo_units).pint.dequantify()
    obs_depo = obs_depo.pint.quantify().pint.to(depo_units).pint.dequantify()

    # Ensure all values are in memory
    sim_conc = sim_conc.compute()
    sim_depo = sim_depo.compute()
    obs_conc = obs_conc.compute()
    obs_depo = obs_depo.compute()

    # Optimize
    options = {"disp": True}
    if maxiter is not None:
        options["maxiter"] = maxiter
    result = optimize.minimize(
        fun=compute_score,
        x0=init.values.flatten(),
        args=(sim_conc, sim_depo, obs_conc, obs_depo),
        method=method,
        bounds=bounds,
        tol=tol,
        options=options,
        **kwargs,
    )

    # Unflatten fitted scales
    scales = reshape_scales(result.x, coords=sim_conc.coords).rename("scale_factor")

    return scales


def reshape_scales(x: NDArray, coords: xr.Coordinates) -> xr.DataArray:
    """Reshape 1D vector of scaling factors to match simulation outputs

    Needed because scipy.optimize.minimize requires initial conditions as 1D array

    Args:
        x: scaling factors (1D)
        coords: Coordinates of simulation outputs

    Returns:
        DataArray of scaling factors with same coordinates as simulation outputs
    """

    dims = [d for d in coords.dims if d != "index"]
    scales = xr.DataArray(
        x.reshape([coords.sizes[d] for d in dims]),
        coords={d: coords[d] for d in dims},
    )

    return scales
