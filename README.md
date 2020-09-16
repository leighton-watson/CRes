# CRes

CRes (Crater Resonance) is a one-dimensional (1D) numerical method for solving the linear acoustic wave equation within a volcanic crater. For a specified crater geometry and excitation source at the base of the crater, CRes computes the velocity and pressure at the crater outlet and can propagate the signal to an infrasound receiver some distance from the outlet. The linear acoustic wave equation is written in terms of two first-order differential equations for pressure and acoustic flow and is solved by CRes using a finite-difference frequency-domain method. CRes is written in MATLAB  and runs efficiently on a standard desktop/laptop computer. 

For more details see: 
* Watson, L. M., Dunham, E. M., and Johnson, J. B. (2019) Simulation and inversion of harmonic infrasound from open-vent volcanoes using an efficient quasi-1D crater model, Journal of Volcanology and Geothermal Research, [https://doi.org:10.1016/j.jvolgeores.2019.05.007](https://doi.org/10.1016/j.jvolgeores.2019.05.007).
* Watson, L. M., Johnson, J. B., Sciotto, M., and Cannata, A. (2020) Changes in crater geometry revealed by inversion of harmonic infrasound observations: 24 December 2018 eruption of Mount Etna, Italy, Geophysical Research Letters.

The first release of this repository is archived at Zenodo: [https://doi.org/10.5281/zenodo.3235683](https://doi.org/10.5281/zenodo.3235683) and contains the codes associated with [Watson et al. (2019)](https://doi.org/10.1016/j.jvolgeores.2019.05.007). 

<a href="https://doi.org/10.5281/zenodo.3235683"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.3235683.svg" alt="DOI"></a>

CRes can be used to inverted harmonic infrasound observations for crater geometry. An example of this procedure applied at Villarrica Volcano, Chile, is featured in the first release with details described in [Watson et al. (2019)](https://doi.org/10.1016/j.jvolgeores.2019.05.007).

The second release of CRes contains the codes associated with [Watson et al. (2020)](https://doi.org/10.1016/j.jvolgeores.2019.05.007). The second release contains more information about inverting harmonic infrasound observations for crater geometry including (1) a user guide that describes the inversion methodology, (2) a synthetic example, and (3) examples of the inversion procedure applied at Mount Etna. The Mount Etna data set can be downloaded from [https://doi.org/10.13127/etna_infra/raw_20181223_25](https://doi.org/10.13127/etna_infra/raw_20181223_25). Note that the data is provided in SAC format and can be converted to loaded into Matlab using rdsac: [https://github.com/IPGP/sac-matlab](https://github.com/IPGP/sac-matlab).

The second release of this repository will be archived at Zenodo and contains the codes associated with Watson et al. (2020). 

### How do I get set up? ###
* Clone this respository to your local directory
* demo contains example script files.
* inversion_demo contains example files for the inversion procedure.
* doc contains the user guide.
* JVGR2019 contains files associated with [Watson et al. (2019)](https://doi.org/10.1016/j.jvolgeores.2019.05.007).
* GRL2020 contains files associated with Watson et al. (2020).

### Who do I talk to? ###
* Leighton Watson: leightonmwatson@gmail.com




