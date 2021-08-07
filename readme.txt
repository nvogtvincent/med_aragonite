################################################################################
#                                 README                                       #
################################################################################

Directories:
- data: data obtained directly from CMIP/ERA5/Observation compilations
- processed_data: processed netcdf files containing omega/sst/rate variables
  calculated from the files in the data directory
- figures: figures generated from scripts

Scripts:
- aragonite_parameters.py: calculates the inorganic aragonite precipitation rate
  parameterisation used in supplement 3
- calculate_aragonite.py: calculates omega and precipitation rate from CMIP data
- calculate_aragonite_monclim.py: generates monthly climatology from the output of
  calculate_aragonite
- model_validation.py: CMIP6-ERA5-Proxy-Obs comparisons for validation
- precipitation_rate_time_series.py: Plots a time series of estimated inorganic
  aragonite precipitation rate for a particular region of the ocean under full,
  sst-only and omega-only scenarios

Prerequisite python packages (conda recommended):
- numpy
- scipy
- matplotlib
- netcdf4
- PyCO2SYS
- seawater
- progress (only for calculate_aragonite.py, dependency can be easily removed)

Terms of use:
- CMIP6 datasets from UKESM1-0-LL, GFDL-CM4 and CanESM5 are included under the
  Creative Commons Attribution-ShareAlike 4.0 International License. Consult
  https://pcmdi.llnl.gov/CMIP6/TermsOfUse for terms of use governing CMIP6
  output. For citations of the CMIP6 datasets used in this study, please see the
  accompanying paper or the bottom of this readme.
- ERA5 data is included under the terms of use for Copernicus Products, please
  see https://cds.climate.copernicus.eu/api/v2/terms/static/licence-to-use-copernicus-products.pdf.

References:
Bialik, O. M., & Sisma-Ventura, G. (2016). Proxy-based reconstruction of surface water acidification and carbonate saturation of the Levant Sea during the Anthropocene. Anthropocene, 16, 42–53. https://doi.org/10.1016/j.ancene.2016.08.001
Burton, E. A., & Walter, L. M. (1987). Relative precipitation rates of aragonite and Mg calcite from seawater: temperature or carbonate ion control? Geology, 15(2), 111–114. https://doi.org/10.1130/0091-7613(1987)15<111:RPROAA>2.0.CO;2
Eyring, V., Bony, S., Meehl, G. A., Senior, C. A., Stevens, B., Stouffer, R. J., & Taylor, K. E. (2016). Overview of the Coupled Model Intercomparison Project Phase 6 (CMIP6) experimental design and organization. Geoscientific Model Development, 9(5), 1937–1958. https://doi.org/10.5194/gmd-9-1937-2016
Gidden, M. J., Riahi, K., Smith, S. J., Fujimori, S., Luderer, G., Kriegler, E., Van Vuuren, D. P., Van Den Berg, M., Feng, L., Klein, D., Calvin, K., Doelman, J. C., Frank, S., Fricko, O., Harmsen, M., Hasegawa, T., Havlik, P., Hilaire, J., Hoesly, R., … Takahashi, K. (2019). Global emissions pathways under different socioeconomic scenarios for use in CMIP6: A dataset of harmonized emissions trajectories through the end of the century. Geoscientific Model Development, 12(4), 1443–1475. https://doi.org/10.5194/gmd-12-1443-2019
Held, I. M., Guo, H., Adcroft, A., Dunne, J. P., Horowitz, L. W., Krasting, J., Shevliakova, E., Winton, M., Zhao, M., Bushuk, M., Wittenberg, A. T., Wyman, B., Xiang, B., Zhang, R., Anderson, W., Balaji, V., Donner, L., Dunne, K., Durachta, J., … Zadeh, N. (2019). Structure and Performance of GFDL’s CM4.0 Climate Model. Journal of Advances in Modeling Earth Systems, 11(11), 3691–3727. https://doi.org/10.1029/2019MS001829
Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz-Sabater, J., Nicolas, J., Peubey, C., Radu, R., Schepers, D., Simmons, A., Soci, C., Abdalla, S., Abellan, X., Balsamo, G., Bechtold, P., Biavati, G., Bidlot, J., Bonavita, M., … Thépaut, J. N. (2020). The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society, 146(730), 1999–2049. https://doi.org/10.1002/qj.3803
Humphreys, M. P., Gregor, L., Pierrot, D., van Heuven, S., Lewis, E. R., & Wallace, D. W. R. (2020). PyCO2SYS: marine carbonate system calculations in Python. Zenodo. https://doi.org/10.5281/ZENODO.3967359
Lawrence, B. N., Bennett, V. L., Churchill, J., Juckes, M., Kershaw, P., Pascoe, S., Pepler, S., Pritchard, M., & Stephens, A. (2013). Storing and manipulating environmental big data with JASMIN. Proceedings - 2013 IEEE International Conference on Big Data, Big Data 2013, 68–75. https://doi.org/10.1109/BigData.2013.6691556
Lewis, E., & Wallace, D. (1998). Program developed for CO2 system calculations. ORNL/CDIAC-105, Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, U.S. Department of Energy. http://cdiac.esd.ornl.gov/oceans/co2rprtnbk.html
O’Neill, B. C., Kriegler, E., Ebi, K. L., Kemp-Benedict, E., Riahi, K., Rothman, D. S., van Ruijven, B. J., van Vuuren, D. P., Birkmann, J., Kok, K., Levy, M., & Solecki, W. (2017). The roads ahead: Narratives for shared socioeconomic pathways describing world futures in the 21st century. Global Environmental Change, 42, 169–180. https://doi.org/10.1016/j.gloenvcha.2015.01.004
O’Neill, B. C., Tebaldi, C., Van Vuuren, D. P., Eyring, V., Friedlingstein, P., Hurtt, G., Knutti, R., Kriegler, E., Lamarque, J. F., Lowe, J., Meehl, G. A., Moss, R., Riahi, K., & Sanderson, B. M. (2016). The Scenario Model Intercomparison Project (ScenarioMIP) for CMIP6. Geoscientific Model Development, 9(9), 3461–3482. https://doi.org/10.5194/gmd-9-3461-2016
Orr, J. C., Najjar, R. G., Aumont, O., Bopp, L., Bullister, J. L., Danabasoglu, G., Doney, S. C., Dunne, J. P., Dutay, J. C., Graven, H., Griffies, S. M., John, J. G., Joos, F., Levin, I., Lindsay, K., Matear, R. J., McKinley, G. A., Mouchet, A., Oschlies, A., … Yool, A. (2017). Biogeochemical protocols and diagnostics for the CMIP6 Ocean Model Intercomparison Project (OMIP). Geoscientific Model Development, 10(6), 2169–2199. https://doi.org/10.5194/gmd-10-2169-2017
Sellar, A. A., Jones, C. G., Mulcahy, J. P., Tang, Y., Yool, A., Wiltshire, A., O’Connor, F. M., Stringer, M., Hill, R., Palmieri, J., Woodward, S., de Mora, L., Kuhlbrodt, T., Rumbold, S. T., Kelley, D. I., Ellis, R., Johnson, C. E., Walton, J., Abraham, N. L., … Zerroukat, M. (2019). UKESM1: Description and Evaluation of the U.K. Earth System Model. Journal of Advances in Modeling Earth Systems, 11(12), 4513–4558. https://doi.org/10.1029/2019MS001739
Sisma-Ventura, G., Bialik, O. M., Yam, R., Herut, B., & Silverman, J. (2017). pCO2 variability in the surface waters of the ultra-oligotrophic Levantine Sea: Exploring the air–sea CO2 fluxes in a fast warming region. Marine Chemistry, 196(February), 13–23. https://doi.org/10.1016/j.marchem.2017.06.006
Sisma-Ventura, G., Yam, R., Kress, N., & Shemesh, A. (2016). Water column distribution of stable isotopes and carbonate properties in the South-eastern Levantine basin (Eastern Mediterranean): Vertical and temporal change. Journal of Marine Systems, 158, 13–25. https://doi.org/10.1016/j.jmarsys.2016.01.012
Sulpis, O., Lauvset, S. K., & Hagens, M. (2020). Current estimates of K1∗and K2∗appear inconsistent with measured CO2 system parameters in cold oceanic regions. Ocean Science, 16(4), 847–862. https://doi.org/10.5194/os-16-847-2020
Swart, N. C., Cole, J. N. S., Kharin, V. V., Lazare, M., Scinocca, J. F., Gillett, N. P., Anstey, J., Arora, V., Christian, J. R., Hanna, S., Jiao, Y., Lee, W. G., Majaess, F., Saenko, O. A., Seiler, C., Seinen, C., Shao, A., Sigmond, M., Solheim, L., … Winter, B. (2019). The Canadian Earth System Model version 5 (CanESM5.0.3). Geoscientific Model Development, 12(11), 4823–4873. https://doi.org/10.5194/gmd-12-4823-2019
