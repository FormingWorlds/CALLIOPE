# Publications

## Methods papers (cite when using CALLIOPE)

If you use CALLIOPE in published work, please cite the following three methods papers, which together describe (i) the original mass-balance + Henry's law framework, (ii) the multi-species redox-coupled extension, and (iii) the magma-ocean evolution context that defined the present species set.

- **Bower, D.J., Kitzmann, D., Wolf, A.S., Sanan, P., Dorn, C., & Oza, A.V. (2019).** Linking the evolution of terrestrial interiors and an early outgassed atmosphere to astrophysical observations. *Astronomy & Astrophysics, 631*, A103. \[[ADS](https://ui.adsabs.harvard.edu/abs/2019A%26A...631A.103B) | [DOI](https://doi.org/10.1051/0004-6361/201935710)\]
- **Bower, D.J., Hakim, K., Sossi, P.A., & Sanan, P. (2022).** Retention of water in terrestrial magma oceans and carbon-rich early atmospheres. *The Planetary Science Journal, 3*(4), 93. \[[ADS](https://ui.adsabs.harvard.edu/abs/2022PSJ.....3...93B) | [DOI](https://doi.org/10.3847/PSJ/ac5fb1)\]
- **Nicholls, H., Lichtenberg, T., Bower, D.J., & Pierrehumbert, R. (2024).** Magma ocean evolution at arbitrary redox state. *Journal of Geophysical Research: Planets, 129*, e2024JE008576. \[[ADS](https://ui.adsabs.harvard.edu/abs/2024JGRE..12908576N) | [DOI](https://doi.org/10.1029/2024JE008576) | [arXiv](https://arxiv.org/abs/2411.19137)\]

## Underlying chemistry and solubility-law sources

CALLIOPE inherits its calibration from the following experimental and thermochemical-fit papers. Cite as appropriate to the species and conditions you exercise.

### Equilibrium constants

- **Chase, M.W. (1998).** *NIST-JANAF Thermochemical Tables*, 4th edition, Journal of Physical and Chemical Reference Data Monograph 9. Source for the JANAF fits used in `janaf_H2`, `janaf_CO`, `janaf_SO2`, `janaf_H2S`, `janaf_NH3`. \[[NIST landing page](https://janaf.nist.gov/)\]
- **Schaefer, L., & Fegley, B. (2017).** Redox states of initial atmospheres outgassed on rocky planets and planetesimals. *The Astrophysical Journal, 843*(2), 120. \[[ADS](https://ui.adsabs.harvard.edu/abs/2017ApJ...843..120S) | [DOI](https://doi.org/10.3847/1538-4357/aa784f)\] (IVTHANTHERMO source for `schaefer_H`, `schaefer_C`, `schaefer_CH4`.)

### Oxygen-fugacity buffers

- **O'Neill, H.St.C., & Eggins, S.M. (2002).** The effect of melt composition on trace element partitioning: an experimental investigation of the activity coefficients of FeO, NiO, CoO, MoO$_2$ and MoO$_3$ in silicate melts. *Chemical Geology, 186*, 151-181. \[[ADS](https://ui.adsabs.harvard.edu/abs/2002ChGeo.186..151O) | [DOI](https://doi.org/10.1016/S0009-2541(01)00414-4)\] (Source for the IW buffer parameterisation `oneill`.)
- **Fischer, R.A., Campbell, A.J., Reaman, D.M., Miller, N.A., Heinz, D.L., Dera, P., & Prakapenka, V.B. (2013).** Phase relations in the Fe-FeSi system at high pressures and temperatures. *Earth and Planetary Science Letters, 373*, 54-64. \[[ADS](https://ui.adsabs.harvard.edu/abs/2013E%26PSL.373...54F) | [DOI](https://doi.org/10.1016/j.epsl.2013.04.035)\] (Source for the alternative IW buffer `fischer`.)
- **Sossi, P.A., Burnham, A.D., Badro, J., Lanzirotti, A., Newville, M., & O'Neill, H.St.C. (2020).** Redox state of Earth's magma ocean and its Venus-like early atmosphere. *Science Advances, 6*, eabd1387. \[[ADS](https://ui.adsabs.harvard.edu/abs/2020SciA....6.1387S) | [DOI](https://doi.org/10.1126/sciadv.abd1387)\] (Reference for the modern Earth $\Delta\mathrm{IW} \approx +3.5$ used as a default.)

### Solubility laws

- **Sossi, P.A., Tollan, P.M.E., Badro, J., & Bower, D.J. (2023).** Solubility of water in peridotite liquids and the prevalence of steam atmospheres on rocky planets. *Earth and Planetary Science Letters, 601*, 117894. \[[ADS](https://ui.adsabs.harvard.edu/abs/2023E%26PSL.60117894S) | [DOI](https://doi.org/10.1016/j.epsl.2022.117894) | [arXiv](https://arxiv.org/abs/2211.13344)\] (CALLIOPE H$_2$O default `peridotite`.)
- **Newcombe, M.E., Brett, A., Beckett, J.R., Baker, M.B., Newman, S., Guan, Y., Eiler, J.M., & Stolper, E.M. (2017).** Solubility of water in lunar basalt at low pH$_2$O. *Geochimica et Cosmochimica Acta, 200*, 330-352. \[[ADS](https://ui.adsabs.harvard.edu/abs/2017GeCoA.200..330N) | [DOI](https://doi.org/10.1016/j.gca.2016.12.026)\] (CALLIOPE H$_2$O `lunar_glass` and `anorthite_diopside`.)
- **Dixon, J.E., Stolper, E.M., & Holloway, J.R. (1995).** An experimental study of water and carbon dioxide solubilities in mid-ocean ridge basaltic liquids. Part I: Calibration and solubility models. *Journal of Petrology, 36*(6), 1607-1631. \[[ADS](https://ui.adsabs.harvard.edu/abs/1995JPet...36.1607D) | [DOI](https://doi.org/10.1093/oxfordjournals.petrology.a037267)\] (CALLIOPE CO$_2$ `basalt_dixon` and H$_2$O `basalt_dixon`.)
- **Hamilton, D.L., Burnham, C.W., & Osborn, E.F. (1964).** The solubility of water and effects of oxygen fugacity and water content on crystallization in mafic magmas. *Journal of Petrology, 5*(1), 21-39. \[[DOI](https://doi.org/10.1093/petrology/5.1.21)\] (Underlying H$_2$O solubility data behind `basalt_wilson`.)
- **Wilson, L., & Head, J.W. (1981).** Ascent and eruption of basaltic magma on the Earth and Moon. *Journal of Geophysical Research, 86*(B4), 2971-3001. \[[ADS](https://ui.adsabs.harvard.edu/abs/1981JGR....86.2971W) | [DOI](https://doi.org/10.1029/JB086iB04p02971)\] (CALLIOPE H$_2$O `basalt_wilson` parametrisation.)
- **Armstrong, L.S., Hirschmann, M.M., Stanley, B.D., Falksen, E.G., & Jacobsen, S.D. (2015).** Speciation and solubility of reduced C-O-H-N volatiles in mafic melt: implications for volcanism, atmospheric evolution, and deep volatile cycles in the terrestrial planets. *Geochimica et Cosmochimica Acta, 171*, 283-302. \[[ADS](https://ui.adsabs.harvard.edu/abs/2015GeCoA.171..283A) | [DOI](https://doi.org/10.1016/j.gca.2015.07.007)\] (CALLIOPE CO solubility `mafic_armstrong`.)
- **Ardia, P., Hirschmann, M.M., Withers, A.C., & Stanley, B.D. (2013).** Solubility of CH$_4$ in a synthetic basaltic melt, with applications to atmosphere-magma ocean-core partitioning of volatiles and to the evolution of the Martian atmosphere. *Geochimica et Cosmochimica Acta, 114*, 52-71. \[[ADS](https://ui.adsabs.harvard.edu/abs/2013GeCoA.114...52A) | [DOI](https://doi.org/10.1016/j.gca.2013.03.028)\] (CALLIOPE CH$_4$ `basalt_ardia`.)
- **Libourel, G., Marty, B., & Humbert, F. (2003).** Nitrogen solubility in basaltic melt. Part I. Effect of oxygen fugacity. *Geochimica et Cosmochimica Acta, 67*(21), 4123-4135. \[[ADS](https://ui.adsabs.harvard.edu/abs/2003GeCoA..67.4123L) | [DOI](https://doi.org/10.1016/S0016-7037(03)00259-X)\] (CALLIOPE N$_2$ `libourel`.)
- **Dasgupta, R., Falksen, E., Pal, A., & Sun, C. (2022).** The fate of nitrogen during parent body partial melting and accretion of the inner Solar System bodies at reducing conditions. *Geochimica et Cosmochimica Acta, 336*, 291-307. \[[ADS](https://ui.adsabs.harvard.edu/abs/2022GeCoA.336..291D) | [DOI](https://doi.org/10.1016/j.gca.2022.09.012)\] (CALLIOPE N$_2$ default `dasgupta`.)
- **Gaillard, F., Bernadou, F., Roskosz, M., Bouhifd, M.A., Marrocchi, Y., Iacono-Marziano, G., Moreira, M., Scaillet, B., & Rogerie, G. (2022).** Redox controls during magma ocean degassing. *Earth and Planetary Science Letters, 577*, 117255. \[[ADS](https://ui.adsabs.harvard.edu/abs/2022E%26PSL.57717255G) | [DOI](https://doi.org/10.1016/j.epsl.2021.117255)\] (CALLIOPE S$_2$ `gaillard`.)

## Applications using CALLIOPE within PROTEUS

These are publications that have applied CALLIOPE within coupled PROTEUS runs.

- **Nicholls, H., Pierrehumbert, R.T., Lichtenberg, T., Soucasse, L., & Smeets, S. (2025).** Convective shutdown in the atmospheres of lava worlds. *Monthly Notices of the Royal Astronomical Society, 536*(3), 2957-2971. \[[ADS](https://ui.adsabs.harvard.edu/abs/2025MNRAS.536.2957N) | [DOI](https://doi.org/10.1093/mnras/stae2772) | [arXiv](https://arxiv.org/abs/2412.11987)\]
- **Hammond, M., Guimond, C.M., Lichtenberg, T., Nicholls, H., Fisher, C., Luque, R., Meier, T.G., Taylor, J., Changeat, Q., Dang, L., Hay, H.C.F.C., Herbort, O., & Teske, J. (2025).** Reliable detections of atmospheres on rocky exoplanets with photometric JWST phase curves. *The Astrophysical Journal Letters, 978*, L40. \[[ADS](https://ui.adsabs.harvard.edu/abs/2025ApJ...978L..40H) | [arXiv](https://arxiv.org/abs/2409.04386)\]
- **Nicholls, H., Lichtenberg, T., Chatterjee, R.D., Guimond, C.M., Postolec, E., & Pierrehumbert, R.T. (2026).** Volatile-rich evolution of molten super-Earth L 98-59 d. *Nature Astronomy*. \[[ADS](https://ui.adsabs.harvard.edu/abs/2026NatAs.tmp...61N) | [DOI](https://doi.org/10.1038/s41550-026-02815-8)\]

## Related software

- **Bower, D.J., Thompson, M.A., Hakim, K., Tian, M., & Sossi, P.A. (2025).** Diversity of low-mass planet atmospheres in the C-H-O-N-S-Cl system with interior dissolution, nonideality, and condensation: application to TRAPPIST-1e and sub-Neptunes. *The Astrophysical Journal, 995*, 59. \[[ADS](https://ui.adsabs.harvard.edu/abs/2025ApJ...995...59B) | [arXiv](https://arxiv.org/abs/2507.00499)\] The atmodeller successor framework, used as the alternative outgassing module within PROTEUS.

CALLIOPE is also part of the wider PROTEUS publication record, the live list for which is at [proteus-framework.org/publications](https://proteus-framework.org/publications).
