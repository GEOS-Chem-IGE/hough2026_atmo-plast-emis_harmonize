Size-harmonized atmospheric microplastics
=========================================

This directory contains size-harmonized observations and constrained simulation outputs from "Reduced global atmospheric microplastic emissions from size-harmonized observations" (Hough et al., 2026).


Contents
--------

[sensitivity/](sensitivity/): artifacts of analyses evaluating how sensitive our results are to the slope of the atmospheric microplastic particle size distribution and to the observations used to constrain simulated emissions. See the "Sensitivity" sections of [size-harmonize-obs.ipynb](/notebooks/size-harmonize-obs.ipynb) and [constrain-simulations.ipynb](/notebooks/constrain-simulations.ipynb) for details.

`obs_fu2023-size-harmonzied.csv`: atmospheric microplastic concentration and deposition observations, size-harmonized to the 0.1-100 µm particle size range. All observations were assumed to have a power law particle number size distribution with $\alpha$ = -2.7 (i.e. the log-log slope of the particle number size distribution is -2.7). Mass was computed assuming ellipsoidal particles with length = size, width = 0.68 length, height = 0.4 width, and density 1 g/cm3. See [size-harmonize-obs.ipynb](/notebooks/size-harmonize-obs.ipynb) for details.

`scales_{main|alt}.nc`: optimized scaling factors that adjust the postprocessed raw outputs of the *main* and *alternate* simulations to match observations. The *main* simulation is constrained by size-harmonized observations while the *alternate* simulation is constrained by unharmonized observations. See [constrain-simulations.ipynb](/notebooks/constrain-simulations.ipynb) for details.

`sim_{main|alt}.nc`: postprocessed raw outputs of the *main* and *alternate* simulations (time-averaged over 2018-2020). See [process-sim-outputs.ipynb](/notebooks/process-sim-outputs.ipynb) for details.

`sim_{main|alt}_constrained.nc`: constrained outputs of the *main* and *alternate* simulations obtained by multiplying the postprocessed raw outputs by the optimized scaling factors. See [constrain-simulations.ipynb](/notebooks/constrain-simulations.ipynb) for details.

`trace_{main|alt}.pkl`: trace of the scale optimization loop for the *main* and *alternate* simulations. See [constrain-simulations.ipynb](/notebooks/constrain-simulations.ipynb) for details.


References
----------

Hough, I., Angot, H., Price, R., Dobiasova, N., Segur, T., Jahangir, E., Zhang, Y., Voisin, D., Sonke, J.E., & Thomas, J.L. (2026) Reduced global atmospheric microplastic emissions from size-harmonized observations.
