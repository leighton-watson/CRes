# JVGR2019

This folder contains the script files needed to recreate figures from Watson, L. M., Dunham, E. M., and Johnson, J. B. (2019) Simulation and inversion of harmonic infrasound from open-vent volcanoes using an efficient quasi-1D crater model, Journal of Volcanology and Geothermal Research, [https://doi.org:10.1016/j.jvolgeores.2019.05.007](https://doi.org/10.1016/j.jvolgeores.2019.05.007).

### Structure of repository ###
* **FigXX_YYYY.m** - script files that generate Figures 3 to 12. Figures 1 and 2 are non-reproducible. 
* **Data** - contains data that are required to generates figures. The data includes text files of crater geometries or outputs from infraFDTD (Kim and Less, 2011 [https://doi.org/10.1029/2010GL046615](https://doi.org/10.1029/2010GL046615)) This folder contains subfolders containing the data required for each individual figure. Note that not all figures require additional data. 
* **Figures** - high resolution .pdf versions of all figures.
* **Functions** - functions that are required to generate figures. These are mostly slightly edited versions of the functions that are contained in `source`.
* **Inversion** - code and data related to the inversion results discussed in Section 7 and shown in Figure 12.

