Size-harmonized atmospheric microplastics
=========================================

This directory contains the size-harmonized observations and constrained simulation outputs described in "Reduced global atmospheric microplastic emissions from size-harmonized observations" (Hough et al., 2026).

The data are archived at DOI:[10.5281/zenodo.22865719](https://doi.org/10.5281/zenodo.22865719)


Contents
--------

[sensitivity/](sensitivity/): artifacts of analyses evaluating how sensitive our results are to the slope of the atmospheric microplastic particle size distribution and to the observations used to constrain simulated emissions. See the "Sensitivity" sections of [size-harmonize-obs.ipynb](/notebooks/size-harmonize-obs.ipynb) and [constrain-simulations.ipynb](/notebooks/constrain-simulations.ipynb) for details.

`obs_fu2023-revised.csv`: observational dataset of Fu et al. (2023), revised to ensure lat/lon and microplastic particle counts (`particle/m3` or `particle/m2/d`) correspond to the study-reported values. Microplastic mass (`ug/m3` or `t/km2/yr`) was recomputed from the updated particle counts assuming a fixed microplastic particle mass of 57 ng/particle over land and 100 ng/particle over oceans. See [prep-obs-fu2023.ipynb](/notebooks/prep-obs-fu2023.ipynb) for details.

`obs_fu2023-size-harmonzied.csv`: observational dataset of Fu et al. (2023) size-harmonized to the 0.1-100 µm particle size range. All observations were assumed to have a power law particle number size distribution with $\alpha$ = -2.7 (i.e. the log-log slope of the particle number size distribution is -2.7). Mass was computed assuming ellipsoidal particles with length = size, width = 0.68 length, height = 0.4 width, and density 1 g/cm3. See [size-harmonize-obs.ipynb](/notebooks/size-harmonize-obs.ipynb) for details.

`scales_{main|alt}.nc`: optimized scaling factors that adjust the postprocessed raw outputs of the *main* and *alternate* simulations to match observations. The *main* simulation is constrained by size-harmonized observations while the *alternate* simulation is constrained by unharmonized observations. See [constrain-simulations.ipynb](/notebooks/constrain-simulations.ipynb) for details.

`sim_{main|alt}.nc`: postprocessed raw outputs of the *main* and *alternate* simulations (time-averaged over 2018-2020). See [process-sim-outputs.ipynb](/notebooks/process-sim-outputs.ipynb) for details.

`sim_{main|alt}_constrained.nc`: constrained outputs of the *main* and *alternate* simulations obtained by multiplying the postprocessed raw outputs by the optimized scaling factors. See [constrain-simulations.ipynb](/notebooks/constrain-simulations.ipynb) for details.

`trace_{main|alt}.pkl`: trace of the scale optimization loop for the *main* and *alternate* simulations. See [constrain-simulations.ipynb](/notebooks/constrain-simulations.ipynb) for details.


References
----------

Fu, Y., Pang, Q., Ga, S. L. Z., Wu, P., Wang, Y., Mao, M., Yuan, Z., Xu, X., Liu, K., Wang, X., Li, D., & Zhang, Y. (2023). Modeling atmospheric microplastic cycle by GEOS-Chem: An optimized estimation by a global dataset suggests likely 50 times lower ocean emissions. *One Earth*, *6*(6), 705–714. https://doi.org/10.1016/j.oneear.2023.05.012

Hough, I., Angot, H., Price, R., Dobiasova, N., Segur, T., Jahangir, E., Zhang, Y., Voisin, D., Sonke, J.E., & Thomas, J.L. (2026) Reduced global atmospheric microplastic emissions from size-harmonized observations.
