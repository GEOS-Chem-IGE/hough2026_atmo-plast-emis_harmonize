# isort: off
import xarray as xr
import pint_xarray
import cf_xarray.units
# isort: on

import numpy as np
from pint_xarray import unit_registry as ureg

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
