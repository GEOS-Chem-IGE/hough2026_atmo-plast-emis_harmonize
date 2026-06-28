Data
====

Original and postprocessed observations of atmospheric microplastic concentration and deposition.


Contents
--------

### [original/](original/)

Original data files (before postprocessing).

`concentration.txt`, `deposition.txt`: observations of atmospheric microplastic concentration and deposition collected from the literature by Fu et al. (2023). These files were provided by Yanxu Zhang. The `number` column contains the number concentration reported by the study. The `mass` contains the corresponding mass; Fu et al. (2023) computed this from the number concentration assuming a fixed microplastic particle mass of 57 ng/particle over land and 100 ng/particle over oceans.

`MPsizeBase v9-1-2026.xlsx`: MPsizeBase v2, a database of environmental microplastic particle size distributions and the corresponding power law parameters (Sonke et al., 2025). Downloaded from [doi:10.5281/zenodo.18201301](https://doi.org/10.5281/zenodo.18201301).

### [processed](processed/)

Postprocessed data derived from the [original](original/) files.

`mpsizebase_atmo_.csv`: power law particle size distribution parameters for atmospheric microplastics extracted from MPsizeBase. See [prep-mpsizebase.ipynb](/notebooks/prep-mpsizebase.ipynb) for details.

`obs_fu2023_*.csv`: Fu et al. (2023)'s collected observations, annotated with study DOI and observed particle size range. See [prep-fu2023-obs.ipynb](/notebooks/prep-fu2023-obs.ipynb) for details.

`obs_revised_*.csv`: the annotated observations, updated to ensure lat/lon and microplastic particle counts (`particle/m3` or `particle/m2/d`) correspond to the study-reported values. Microplastic mass (`ug/m3` or `t/km2/yr`) was recomputed from the updated particle counts assuming a fixed microplastic particle mass of 57 ng/particle over land and 100 ng/particle over oceans. See [prep-fu2023-obs.ipynb](/notebooks/prep-fu2023-obs.ipynb) for details.


References
----------

Fu, Y., Pang, Q., Ga, S. L. Z., Wu, P., Wang, Y., Mao, M., Yuan, Z., Xu, X., Liu, K., Wang, X., Li, D., & Zhang, Y. (2023). Modeling atmospheric microplastic cycle by GEOS-Chem: An optimized estimation by a global dataset suggests likely 50 times lower ocean emissions. *One Earth*, *6*(6), 705–714. https://doi.org/10.1016/j.oneear.2023.05.012

Sonke, J., Segur, T., Hough, I., Dobiasova, N., Voisin, D., Yakovenko, N., Margenat, H., Hagelskjaer, O., Abbasi, S., Bucci, S., Richon, C., Angot, H., Thomas, J. L., & Roux, G. L. (2025). MPsizeBase: A database for particle size distributed environmental microplastic data. EarthArXiv. https://doi.org/10.31223/X5XX7R
