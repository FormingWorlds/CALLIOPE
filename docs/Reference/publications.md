# Publications

## Methods papers (cite when using CALLIOPE)

If you use CALLIOPE in published work, please cite the following three methods papers, which together describe (i) the original mass-balance + Henry's law framework, (ii) the multi-species redox-coupled extension, and (iii) the magma-ocean evolution context that defined the present species set.

- **Bower, D.J., Kitzmann, D., Wolf, A.S., Sanan, P., Dorn, C., & Oza, A.V. (2019).** Linking the evolution of terrestrial interiors and an early outgassed atmosphere to astrophysical observations. *Astronomy & Astrophysics, 631*, A103. [https://doi.org/10.1051/0004-6361/201935710](https://doi.org/10.1051/0004-6361/201935710)
- **Bower, D.J., Hakim, K., Sossi, P.A., & Sanan, P. (2022).** Retention of water in terrestrial magma oceans and carbon-rich early atmospheres. *The Planetary Science Journal, 3*(4), 93. [https://doi.org/10.3847/PSJ/ac5fb1](https://doi.org/10.3847/PSJ/ac5fb1)
- **Nicholls, H., Lichtenberg, T., Bower, D.J., & Pierrehumbert, R. (2024).** Magma ocean evolution at arbitrary redox state. *Journal of Geophysical Research: Planets, 129*, e2024JE008576. [https://doi.org/10.1029/2024JE008576](https://doi.org/10.1029/2024JE008576)

## Underlying chemistry and solubility-law sources

CALLIOPE inherits its calibration from the following experimental and thermochemical-fit papers. Cite as appropriate to the species and conditions you exercise.

### Equilibrium constants

- **Chase, M.W. (1998).** *NIST-JANAF Thermochemical Tables*, 4th edition, Journal of Physical and Chemical Reference Data Monograph 9. (Source for the JANAF fits used in `janaf_H2`, `janaf_CO`, `janaf_SO2`, `janaf_H2S`, `janaf_NH3`.)
- **Schaefer, L., & Fegley, B. (2017).** Redox states of initial atmospheres outgassed on rocky planets and planetesimals. *The Astrophysical Journal, 843*(2), 120. [https://doi.org/10.3847/1538-4357/aa784f](https://doi.org/10.3847/1538-4357/aa784f) (IVTHANTHERMO source for `schaefer_H`, `schaefer_C`, `schaefer_CH4`.)

### Oxygen-fugacity buffers

- **O'Neill, H.St.C., & Eggins, S.M. (2002).** The effect of melt composition on trace element partitioning: an experimental investigation of the activity coefficients of FeO, NiO, CoO, MoO$_2$ and MoO$_3$ in silicate melts. *Chemical Geology, 186*, 151-181. [https://doi.org/10.1016/S0009-2541(01)00414-4](https://doi.org/10.1016/S0009-2541(01)00414-4) (Source for the IW buffer parameterisation `oneill`.)
- **Fischer, R.A., Campbell, A.J., Reaman, D.M., Miller, N.A., Heinz, D.L., Dera, P., & Prakapenka, V.B. (2011).** Phase relations in the Fe-FeSi system at high pressures and temperatures. *Earth and Planetary Science Letters, 373*, 54-64. [https://doi.org/10.1016/j.epsl.2013.04.035](https://doi.org/10.1016/j.epsl.2013.04.035) (Source for the alternative IW buffer `fischer`.)
- **Sossi, P.A., Burnham, A.D., Badro, J., Lanzirotti, A., Newville, M., & O'Neill, H.St.C. (2020).** Redox state of Earth's magma ocean and its Venus-like early atmosphere. *Science Advances, 6*, eabd1387. [https://doi.org/10.1126/sciadv.abd1387](https://doi.org/10.1126/sciadv.abd1387) (Reference for the modern Earth $\Delta\mathrm{IW} \approx +3.5$ used as a default.)

### Solubility laws

- **Sossi, P.A., Tollan, P.M.E., Badro, J., Bower, D.J. (2022).** Solubility of water in peridotite liquids and the prevalence of steam atmospheres on rocky planets. *Earth and Planetary Science Letters, 601*, 117894. [https://doi.org/10.1016/j.epsl.2022.117894](https://doi.org/10.1016/j.epsl.2022.117894) (CALLIOPE H$_2$O default `peridotite`.)
- **Newcombe, M.E., Brett, A., Beckett, J.R., Baker, M.B., Newman, S., Guan, Y., Eiler, J.M., & Stolper, E.M. (2017).** Solubility of water in lunar basalt at low pH$_2$O. *Geochimica et Cosmochimica Acta, 200*, 330-352. [https://doi.org/10.1016/j.gca.2016.12.026](https://doi.org/10.1016/j.gca.2016.12.026) (CALLIOPE H$_2$O `lunar_glass` and `anorthite_diopside`.)
- **Dixon, J.E., Stolper, E.M., & Holloway, J.R. (1995).** An experimental study of water and carbon dioxide solubilities in mid-ocean ridge basaltic liquids. Part I: Calibration and solubility models. *Journal of Petrology, 36*(6), 1607-1631. [https://doi.org/10.1093/oxfordjournals.petrology.a037267](https://doi.org/10.1093/oxfordjournals.petrology.a037267) (CALLIOPE CO$_2$ `basalt_dixon` and H$_2$O `basalt_dixon`.)
- **Hamilton, D.L. (1964).** The solubility of water in melts of basaltic compositions. (CALLIOPE H$_2$O `basalt_wilson` background.)
- **Wilson, L., & Head, J.W. (1981).** Ascent and eruption of basaltic magma on the Earth and Moon. *Journal of Geophysical Research, 86*(B4), 2971-3001. [https://doi.org/10.1029/JB086iB04p02971](https://doi.org/10.1029/JB086iB04p02971) (CALLIOPE H$_2$O `basalt_wilson`.)
- **Armstrong, K., Frost, D.J., McCammon, C.A., Rubie, D.C., & Boffa Ballaran, T. (2015).** Deep magma ocean formation set the oxidation state of Earth's mantle. *Science, 365*, 903-906. (CALLIOPE CO solubility `mafic_armstrong`.)
- **Ardia, P., Hirschmann, M.M., Withers, A.C., & Stanley, B.D. (2013).** Solubility of CH$_4$ in a synthetic basaltic melt, with applications to atmosphere-magma ocean-core partitioning of volatiles and to the evolution of the Martian atmosphere. *Geochimica et Cosmochimica Acta, 114*, 52-71. [https://doi.org/10.1016/j.gca.2013.03.028](https://doi.org/10.1016/j.gca.2013.03.028) (CALLIOPE CH$_4$ `basalt_ardia`.)
- **Libourel, G., Marty, B., & Humbert, F. (2003).** Nitrogen solubility in basaltic melt. Part I. Effect of oxygen fugacity. *Geochimica et Cosmochimica Acta, 67*(21), 4123-4135. [https://doi.org/10.1016/S0016-7037(03)00259-X](https://doi.org/10.1016/S0016-7037(03)00259-X) (CALLIOPE N$_2$ `libourel`.)
- **Dasgupta, R., Falksen, E., Pal, A., & Tsuno, K. (2022).** Nitrogen partitioning between silicate melt and ferroan- to magnesiowuestitic-spinel/oxide solids - implications for nitrogen storage in rocky planetary interiors. *Geochimica et Cosmochimica Acta, 336*, 291-307. [https://doi.org/10.1016/j.gca.2022.09.012](https://doi.org/10.1016/j.gca.2022.09.012) (CALLIOPE N$_2$ default `dasgupta`.)
- **Gaillard, F., Bouhifd, M.A., Furi, E., Malavergne, V., Marrocchi, Y., Noack, L., Ortenzi, G., Roskosz, M., & Vulpius, S. (2022).** The diverse planetary ingassing/outgassing paths produced over billions of years of magmatic activity. *Earth and Planetary Science Letters, 583*, 117440. [https://doi.org/10.1016/j.epsl.2021.117255](https://doi.org/10.1016/j.epsl.2021.117255) (CALLIOPE S$_2$ `gaillard`.)

## Applications using CALLIOPE within PROTEUS

These are publications that have applied CALLIOPE within coupled PROTEUS runs.

- **Nicholls, H., Pierrehumbert, R.T., Lichtenberg, T., Soucasse, L., & Smeets, S. (2025).** Convective shutdown in the atmospheres of lava worlds. *Monthly Notices of the Royal Astronomical Society, 536*(3), 2957-2971. [https://doi.org/10.1093/mnras/stae2772](https://doi.org/10.1093/mnras/stae2772)
- **Nicholls, H., Lichtenberg, T., Chatterjee, R.D., Guimond, C.M., Postolec, E., & Pierrehumbert, R.T. (2026).** Volatile-rich evolution of molten super-Earth L 98-59 d. *Nature Astronomy*. [https://doi.org/10.1038/s41550-026-02815-8](https://doi.org/10.1038/s41550-026-02815-8)
- **Hammond, M., Guimond, C.M., Lichtenberg, T., Nicholls, H., Fisher, C., Luque, R., Meier, T.G., Taylor, J., Changeat, Q., Dang, L., Herbort, O., & Teske, J. (2025).** Reliable detections of atmospheres on rocky exoplanets with photometric JWST phase curves. *Nature Astronomy*.

CALLIOPE is also part of the wider PROTEUS publication record, the live list for which is at [proteus-framework.org/publications](https://proteus-framework.org/publications).
