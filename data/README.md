Data
====

Original and postprocessed observations of atmospheric microplastic concentration and deposition.


Contents
--------

### [original/](original/)

Original data files (before postprocessing).

`evangelou2026.csv`: observations of atmospheric microplastic concentration and deposition collected from the literature by Evangelou et al. (2026). This file was downloaded from [doi:10.5281/zenodo.17660114](https://doi.org/10.5281/zenodo.17660114) (original filename: atmo_MP_dataset.csv).

`fu2023_concentration.txt`, `fu2023_deposition.txt`: observations of atmospheric microplastic concentration and deposition collected from the literature by Fu et al. (2023). These files were provided by Yanxu Zhang. The `number` column contains the number concentration reported by the study. The `mass` contains the corresponding mass; Fu et al. (2023) computed this from the number concentration assuming a fixed microplastic particle mass of 57 ng/particle over land and 100 ng/particle over oceans.

`MPsizeBase v19-5-2026.xlsx`: MPsizeBase v5, a database of environmental microplastic particle size distributions and the corresponding power law parameters (Sonke et al., 2025). Downloaded from [doi:10.5281/zenodo.20832480](https://doi.org/10.5281/zenodo.20832480).

### [processed](processed/)

Postprocessed data derived from the [original](original/) files.

`evangelou2026-revised.xlsx`: Evangelou et al. (2026)'s observation dataset, revised to correct some study authors, DOIs, lat/lon, and microplastic particle counts. See [prep-obs-evangelou2026.ipynb](/notebooks/prep-obs-evangelou2026.ipynb) for details.

`mpsizebase_atmo_2026-05-19.csv`: power law particle size distribution parameters for atmospheric microplastics extracted from MPsizeBase. See [prep-mpsizebase.ipynb](/notebooks/prep-mpsizebase.ipynb) for details.

`obs_evangelou2026-revised.csv`: revised observational dataset of Evangelou et al. (2026), aggregated to a single mean concentration or deposition value per sample location and microplastic shape. See [prep-obs-evangelou2026.ipynb](/notebooks/prep-obs-evangelou2026.ipynb) for details.

`obs_fu2023.csv`: observational dataset of Fu et al. (2023), annotated with study DOI and observed particle size range. See [prep-obs-fu2023.ipynb](/notebooks/prep-obs-fu2023.ipynb) for details.


References
----------

Evangelou, I., Bucci, S., & Stohl, A. (2026). Atmospheric microplastic emissions from land and ocean. *Nature*, *649*(8099), 1186–1189. https://doi.org/10.1038/s41586-025-09998-6

Fu, Y., Pang, Q., Ga, S. L. Z., Wu, P., Wang, Y., Mao, M., Yuan, Z., Xu, X., Liu, K., Wang, X., Li, D., & Zhang, Y. (2023). Modeling atmospheric microplastic cycle by GEOS-Chem: An optimized estimation by a global dataset suggests likely 50 times lower ocean emissions. *One Earth*, *6*(6), 705–714. https://doi.org/10.1016/j.oneear.2023.05.012

Sonke, J., Segur, T., Hough, I., Dobiasova, N., Voisin, D., Yakovenko, N., Margenat, H., Hagelskjaer, O., Abbasi, S., Bucci, S., Richon, C., Angot, H., Thomas, J. L., & Roux, G. L. (2025). MPsizeBase: A database for particle size distributed environmental microplastic data. *EarthArXiv*. https://doi.org/10.31223/X5XX7R
