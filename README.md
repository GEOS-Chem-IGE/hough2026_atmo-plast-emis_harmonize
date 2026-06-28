Size-harmonize atmospheric microplastics
========================================

[![DOI](https://zenodo.org/badge/1277827810.svg)](https://zenodo.org/badge/latestdoi/1277827810)

This repository contains code for the paper "Reduced global atmospheric microplastic emissions from size-harmonized observations" (Hough et al., 2026). The code:

1. Size-harmonizes observations of atmospheric microplastic concentration and deposition
2. Constrains modeled atmospheric microplastics with the size-harmonized observations

The size-harmonzied observations are in the [results/](results/) directory (`obs_size-harmonized_*.csv`)

The constrained simulation outputs are available at https://doi.org/10.5281/zenodo.20922804


To reproduce
------------

### Create and activate the environment

You will need [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) and [conda-lock](https://conda.github.io/conda-lock/) to create the environment. You may also use [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) or [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) instead of micromamba.

```bash
conda-lock install --micromamba -n hough2026_atmo-plast-emis_harmonize
micromamba activate hough2026_atmo-plast-emis_harmonize
```

### Run the size-harmonization notebooks

These notebooks size-harmonize observed atmospheric microplastic concentration and deposition to the 0.1-100 µm size range used in the simulations. The size-harmonized data are already in [results](results/) (`obs_size-harmonized_*.csv`); running these notebooks will overwrite them.

> [!TIP]
> If you get a `No such file or directory` error when running the notebooks, try updating your configuration to set the root directory for Jupyter Notebooks to the root of this repository. Alternatively, you can create symlinks in the [notebooks](notebooks/) directory to the [data](data/), [results](results/), [simulations](simulations/), and [utils](utils/) directories.

1. Run [prep-fu2023-obs.ipynb](notebooks/prep-fu2023-obs.ipynb) to postprocess Fu et al. (2023)'s collected observations of atmospheric microplastic concentration and deposition, adding source DOI and observed particle size range

2. Run [prep-mpsizebase.ipynb](notebooks/prep-mpsizebase.ipynb) to extract power law parameters for atmospheric microplastic number particle size distributions from the MPsizeBase database (Sonke et al., 2025)

3. Run [size-harmonize-obs.ipynb](notebooks/size-harmonize-obs.ipynb) to size-harmonize the postprocessed observations to the 0.1-100 µm size range used in the simulations using the median power law slope (parameter $\alpha$) of atmospheric microplastics from MPsizeBase

### Get the raw simulation outputs

The raw (unconstrained) simulation outputs are available at https://doi.org/10.5281/zenodo.20847720. Extract the `main.zip` and `alt.zip` files and place the resulting `main/` and `alt/` directories in `simulations/` like this:

```
simulations/
├── alt/
│   ├── 1-ocen/
│   ├── 2-mmpw/
│   ├── 3-agri/
│   ├── 4-resi/
│   ├── 5-road/
│   └── README.md
└── main/
    ├── 1-ocen/
    ├── 2-mmpw/
    ├── 3-agri/
    ├── 4-resi/
    ├── 5-road/
    └── README.md
```

Alternatively, you can recreate the raw outputs by running the simulation code at https://doi.org/10.5281/zenodo.21068954.

### Run the optimization notebooks

These notebooks constrain the raw simulation outputs to match as closely as possible the size-harmonized observations.

1. Run [process-sim-outputs.ipynb](notebooks/process-sim-outputs.ipynb): to postprocess the raw simulation outputs.

2. Run [constrain-simulations.ipynb](notebooks/constrain-simulations.ipynb): to compute scaling factors that adjust the postprocessed raw simulation outputs to minimize the difference from the size-harmonized observations.

### Extra

You do not need to run the [compute-emission-factors.ipynb](notebooks/compute-emission-factors.ipynb) notebook, which computes the initial base emissions factors and size distribution factors used in the *main* simulation. The file is included here for reference.


Contents
--------

### [data/](data/)

Original and postprocessed datasets

### [notebooks/](notebooks/)

Notebooks to size-harmonize observations and constrain raw simulation outputs

### [results/](results/)

Size-harmonized observations and constrained simulation outputs produced by the notebooks

### [simulations/](simulations/)

Raw simulation outputs

### [utils/](utils/)

Functions used in the notebooks


References
----------

Fu, Y., Pang, Q., Ga, S. L. Z., Wu, P., Wang, Y., Mao, M., Yuan, Z., Xu, X., Liu, K., Wang, X., Li, D., & Zhang, Y. (2023). Modeling atmospheric microplastic cycle by GEOS-Chem: An optimized estimation by a global dataset suggests likely 50 times lower ocean emissions. *One Earth*, *6*(6), 705–714. https://doi.org/10.1016/j.oneear.2023.05.012

Hough, I., Angot, H., Price, R., Dobiasova, N., Segur, T., Jahangir, E., Zhang, Y., Voisin, D., Sonke, J.E., & Thomas, J.L. (2026) Reduced global atmospheric microplastic emissions from size-harmonized observations.

Sonke, J., Segur, T., Hough, I., Dobiasova, N., Voisin, D., Yakovenko, N., Margenat, H., Hagelskjaer, O., Abbasi, S., Bucci, S., Richon, C., Angot, H., Thomas, J. L., & Roux, G. L. (2025). MPsizeBase: A database for particle size distributed environmental microplastic data. EarthArXiv. https://doi.org/10.31223/X5XX7R
