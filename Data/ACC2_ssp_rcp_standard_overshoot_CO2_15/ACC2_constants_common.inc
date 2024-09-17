********** ACC2 version 4.2: ACC2_constants_common.inc *****************************************************************
********** Scaler and parameter declarations and assignments (past and future modes) ===================================

* Note: Values of some scalars and parameters are assigned at the end of this file.


*** Scalars definitions and assignments ================================================================================

SCALARS

*** Flag to control model and inversion setup --------------------------------------------------------------------------

MODEFLAG       To control the program mode (past or future) (see ACC2_settings_past.inc and ACC_settings_future.inc) [1]

* For the usages for xxxSWCH, see ACC2_switchboard_common.inc.
CO2FFSWCH

CO2LUSWCH

ANTCH4N2OSWCH

RADMISSWCH


TMIX1750SWCH

DMPTMIXSWCH

BETASWCH

Q10BASSWCH

Q10TMPSWCH

TNPPSENSSWCH

CO2OCNUPPRESWCH

CO2LNDUPPRESWCH

NATCH4N2OSWCH

TAUCH4N2OSWCH

RAD2XCO2SWCH

SCLRADSWCH

SENSSWCH

KAPPASWCH

RADMISAR1PROSWCH

T2MAR1PROSWCH


CO2UPANTSWCH

CONCO2SWCH

CONCH4N2OSWCH

T2MSWCH


DOECLIMSWCH


DYNCARBSWCH

OCNCO2TFSWCH

LNDCO2TFNPPSWCH

LNDCO2TFRESSWCH

LNDCO2TFQ10SWCH

CO2UPPRESWCH


CO2LUSGMSWCH

CONCO2SGMSWCH

ANTCH4N2OSGMSWCH

CONCH4SGMSWCH

CONN2OSGMSWCH

ENSOCO2SWCH

ENSOT2MSWCH

VOLSWCH

RADMISAR1SWCH

T2MAR1SWCH


MISTRNSWCH

SOLSWCH


DELTAT          Time step of model run                                                   [year] /1/
* The model is written only for the time step of 1 year. Do not alter the time step assignment above.


*** CO2 concentration --------------------------------------------------------------------------------------------------

$offdigit
CNVCO2          Conversion factor for atmospheric CO2 (GtC) -> (ppm)            [ppm per GtC]   /0.471380471380471380/
* The value above is calculated from the assumed preindustrial CO2 atmospheric load (594 GtC)
* divided by the assumed preindustrial CO2 concentration (280 ppm).
$ondigit


*** Ocean CO2 uptake ---------------------------------------------------------------------------------------------------
*
* An unusual unit [Giga L] is used in order to establish a well-scaled problem.

ALK             Average alkalinity                                              [eq per Giga L] /2.435e6/
* Hoffert et al. (1981) pp.290-291. Also refer to Sundquist and Plummer (1981) p.267

BTS             Total borate concentration                                      [mol per Giga L] /0.409e6/
* Hoffert et al. (1981) pp.290-291. Also refer to Sundquist and Plummer (1981) p.267

ALPHAC          Conversion factor for surface layer DIC from mol per Giga L to GtC [GtC per (mol per Giga L)]
                                                                                                /2.42214338e-4/

HS              Depth of oceanic mixed layer                                    [m]             /73/

SAL             Average salinity in the mixed layer                             [1]             /34.76/
* Hoffert et al. (1981) pp.290-291. Also refer to Sundquist and Plummer (1981) p.267

DENS            Average seawater density in the mixed layer                     [kg per L]      /1.025/
* Hoffert et al. (1981) pp.290-291. Also refer to Sundquist and Plummer (1981) p.267

REV1750         Preindustrial Revelle factor                                    [1]             /9/
* Revelle and Munk (1977) cited in Mackenzie and Lerman (2006) p.265

REVCOF          Coefficient for Revelle factor parameterization                 [1]             /4/
* Revelle and Munk (1977) cited in Mackenzie and Lerman (2006) p.265


*** Anthropogenic land CO2 uptake --------------------------------------------------------------------------------------

FNPP            Preindustrial net primary production of land biosphere          [GTC per year]  /60/


*** Concentrations of non-CO2 gases ------------------------------------------------------------------------------------

* Methane

CNVCH4          Conversion factor from mass (TgCH4) to concentration (ppb)      [TgCH4 per ppb] /2.746/

TAUCH4SS        Methane lifetime with respect to stratospheric depletion and soil uptake [year] /68.2/


* N2O

CNVN2O          Conversion factor from mass (TgN) to concentration (ppb)        [TgN per ppb]   /4.8/

CONN2O2000R     Reference concentration of nitrous oxide in 2000                [ppb]           /330/
* Note that the value here is different from the value computed from the past run.
* The reference value is used only in the N2O concentration equation CONN2OQ to scale N2O lifetime.
* This treatment is for consistency to the sensitivity experiment used for calibration (see IPCC (2001) Table 4.5).

SENTAUN2O       Sensitivity coefficient of N2O lifetime                         [ln(year) per ln(ppb)] /-0.046/


* O3TRP

CONO3TRPPRE     Reference concentration of preindustrial tropospheric ozone     [DU]            /25/


* SF6

CNVSF6          Conversion factor from mass (Gg) to concentration (ppt)         [Gg per ppt]    /25.1/

TAUSF6          SF6 lifetime                                                    [year]          /3200/


* SO2

NATSO2          Natural SO2 emission                                            [TgS per year]  /42/


* EESC

AFRCFC11        Absolute fractional release of CFC-11                           [1]             /0.8/

ALPEESC         Relative impact of Br compared to Cl                            [1]             /45/


*** Radiative Forcings -------------------------------------------------------------------------------------------------

OVERLFAC1       CH4-N2O overlap modeling factor                                 [W per m^2]     /0.47/

OVERLFAC2       CH4-N2O overlap modeling factor                                 [1]             /2.01e-5/

OVERLFAC3       CH4-N2O overlap modeling factor                                 [1]             /5.31e-15/

OVERLEXP1       CH4-N2O overlap exponent                                        [1]             /0.75/

OVERLEXP2       CH4-N2O overlap exponent                                        [1]             /1.52/

CONCH4PRE       Preindustrial CH4 concentration                                 [ppb]           /700/
* (IPCC (2001) Table 6.1)

RHOCH4          Coefficient for CH4 radiative forcing calculation               [W per (m^2*ppb^0.5)] /0.036/
* Note to the unit (see the Joos et al. (2001) equation A8).

RHON2O          Coefficient for N2O radiative forcing calculation               [W per (m^2*ppb^0.5)] /0.12/
* Note to the unit (see the Joos et al. (2001) equation A10).

RHOO3TRP        Radiative efficiency of O3TRP                                   [W per (m^2*DU)] /0.042/

ALPHCH4H2O      Factor for stratospheric H2O production from CH4                [1]             /0.035/

CONN2OPRE       Preindustrial N2O concentration                                 [ppb]           /270/
* (IPCC (2001) Table 6.1)

RHOO3STR        Radiative efficiency of O3STR                                   [W per (m^2*ppt)] /-7.32e-5/

CONEESC1970     EESC in 1970                                                    [ppt]           /1250/

RHOSF6          Radiative efficiency of SF6                                     [W per (m^2*ppt)] /5.2e-4/

RADSO2DR2000    Direct effect of sulfur aerosol forcing in 2000                 [W per m^2]     /-0.4/

RADSO2ID2000    Indirect effect of sulfur aerosol forcing in 2000               [W per m^2]     /-0.8/

RADCRBAER2000   Carbonaceous aerosol forcing in 2000                            [W per m^2]     /-0.1/


*** Sea level rise -----------------------------------------------------------------------------------------------------

PCSIPCC         IPCC Central rate of sea level rise due to thermal expansion between 1910 and 1990 [m]   /0.0005/
* (see IPCC TAR Table 11.10)

SLRGLWT1990     Sea level rise (glaciers and small ice caps) without glacier area contraction in 1990 [m] /0.027/

SLRGR1990       Sea level rise due to Greenland icesheet up to 1990                              [m]     /0.00400/
*               (between 1910 and 1990 in IPCC TAR Table 11.10)

SLRAN1990       Sea level rise due to Antarctica icesheet up to 1990                             [m]     /-0.0080/
*               (between 1910 and 1990 in IPCC TAR Table 11.10)

T2M1990AOGCM    Temperature in 1990 relative to the late 19th century                            [degree in C] /0.78/

T2M19ST         Temperature of the late 19th century relative to glacier steady state temperature [degree in C] /0.15/
* These two values were taken from IPCC TAR p.682 and 683.

SENSGL          Sensitivity of global glacier mass balance to global-mean temperature change             /0.000637/
*               [m per year per degree in C]

COEFGL1         Coefficient for sea level equation with glacier area contraction                 [1]     /0.934/

COEFGL2         Coefficient for sea level equation with glacier area contraction                 [per m] /-1.165/

SENSGR          Sensitivity of global Greenland mass balance to global-mean temperature change           /0.000029/
*               [m per year per degree in C]

SENSAN          Sensitivity of global Antarctica mass balance to global-mean temperature change          /-0.000478/
*               [m per year per degree in C]

SLRI            Sea level rise due to ongoing icesheet adjustment to past climate change         [m per year] /0.00025/

SLRP            Sea level rise due to runoff from thawing of permafrost                          [m per year] /0.00011/

SLRS            Sea level rise due to deposition of sediment on the ocean floor                  [m per year] /0.000025/


*** DOECLIM ------------------------------------------------------------------------------------------------------------

PI              Pi                                                              [1]             /3.1416/

* Ocean model parameters

SHC             Specific heat capacity of 1 m^3 seawater                        [W per m^3 K]   /0.13/

CKD             Conversion factor from [cm^2 per s] to [m^2 per a]                              /3155.8/

DEPTH           Depth of Mixed layer                                            [m]             /60/

BOT             Depth of interior ocean                                         [m]             /4000/

RML             Ratio of mixed layer area at surface to that at bottom          [1]             /0.95/

BSI             Ocean surface air warming enhancement                           [1]             /1.3/

* Heat capacities of land and ocean box

CAL             Heat capacity of land atmosphere                                [W per m^2 K]   /0.52/

CML             Heat capacity of mixed layer atmosphere                         [W per m^2 K]   /7.8/

* Land - sea heat exchange parameterization

AK              Parameter for land-sea heat exchange                            [1]             /0.31/

BK              Parameter for land-sea heat exchange                            [W per m^2 K]   /1.59/

* Parameterization of land-ocean differences

FLA             Global land fraction                                            [1]             /0.29/

RLAM            Enhancement of climate sensitivity over land                    [1]             /1.43/

ALBL            Earth albedo over the land                                      [1]             /0.35/

ALBO            Earth albedo over the ocean                                     [1]             /0.28/

DQL             Enhancement in insolation over the land                         [1]             /0.96/

DQO             Enhancement in insolation over the ocean                        [1]             /1.02/

ROZO            Enhancement in ozone forcing over the land                      [1]             /1.5/

RSUL            Enhancement in direct sulphate forcing over the land            [1]             /2.0/

RFF             Enhancement in carbonaceous aerosol forcing over the land       [1]             /3.0/


SWITCH          Counter for solver iteration                                    [1]


SOICOFL         Coefficient for SOI for land surface temperature                                /-0.148/

SOICOFO         Coefficient for SOI for ocean surface temperature                               /-0.063/

* Parameterization of land-ocean differences

GSA             Global surface area                                             [million km2]   /510.072/;
* Reference https://www.cia.gov/library/publications/the-world-factbook/geos/xx.html


*** Parameters definitions and assignments =============================================================================


PARAMETERS

*** Ocean CO2 uptake ---------------------------------------------------------------------------------------------------

HL(NOC)         Layer thicknesses of the box model                              [m]
                / 1 768
                  2 479
                  3 1299
                  4 2723/
* HL('1') is the effective thickness of the combined layer (= Atmosphere + oceanic mixed layer).
* R01 version of HAMOCC run (see Hooss, 2001, Figure 3.1 and Table 3.2) is adopted above.
* When you choose a different version, be sure to change the assignment of HS as well.
* For example, R93 version of HAMOCC run (see Hooss, 2001, Figure 3.1 and Table 3.2) are as follows:
* 1 663.73475
* 2 262.1235
* 3 1641.0377
* 4 6118.5473
* Note that R93 version would give too slow uptake.

ETA(NOC)        Newton diffusion transfer coefficient between layer NOC-1 and layer NOC [m per year]
                / 1 0
                  2 19.3
                  3 10.33
                  4 7.23/
* R01 version of HAMOCC run (see Hooss, 2001, Figure 3.1 and Table 3.2) is adopted above.


*** Land CO2 uptake ----------------------------------------------------------------------------------------------------

TAU(NTR)        Decay time scales                                               [year]
                / 1 2.18182
                  2 2.8575532
                  3 20.0
                  4 100.0/

BJ(NTR)         Coefficients for r_Decay                                        [1 per year]
                / 1  -0.71846
                  2  0.70211
                  3  0.013414
                  4  0.0029323/


*** Halocarbons --------------------------------------------------------------------------------------------------------

CNV(HAL)        Mass (Gg) to concentration (ppt) conversion factors for halocarbons [Gg per ppt]
                /CF4            15.1
                C2F6            23.7
                C4F10           40.8
                HFC23           12.0
                HFC32           8.9
                HFC4310mee      43.2
                HFC125          20.6
                HFC134a         17.5
                HFC143a         14.4
                HFC152a         11.3
                HFC227ea        29.2
                HFC236fa        26.1
                HFC245ca        23.0
                CFC11           23.6
                CFC12           20.7
                CFC113          32.1
                CFC114          29.3
                CFC115          26.5
                CCl4            26.4
                CH3CCl3         22.9
                Halon1211       28.4
                Halon1301       25.6
                Halon2402       44.6
                HCFC22          14.8
                HCFC141b        20.1
                HCFC142b        17.2
                HCFC123         26.2
                CH3Cl           8.7
                CH3Br           16.3/

TAUATM(HAL)     Lifetimes of halocarbons                                        [year]
                /CF4            50000
                C2F6            10000
                C4F10           2600
                HFC23           270
                HFC32           4.9
                HFC4310mee      15.9
                HFC125          29
                HFC134a         14
                HFC143a         52
                HFC152a         1.4
                HFC227ea        34.2
                HFC236fa        240
                HFC245ca        6.2
                CFC11           45
                CFC12           100
                CFC113          85
                CFC114          300
                CFC115          1700
                CCl4            26
                CH3CCl3         5.0
                Halon1211       16
                Halon1301       65
                Halon2402       20
                HCFC22          12
                HCFC141b        9.3
                HCFC142b        17.9
                HCFC123         1.3
                CH3Cl           1.3
                CH3Br           0.7/

NAT(HAL)        Natural emissions of halocarbons                                [Gg]
                /CF4            0
                C2F6            0
                C4F10           0
                HFC23           0
                HFC32           0
                HFC4310mee      0
                HFC125          0
                HFC134a         0
                HFC143a         0
                HFC152a         0
                HFC227ea        0
                HFC236fa        0
                HFC245ca        0
                CFC11           0
                CFC12           0
                CFC113          0
                CFC114          0
                CFC115          0
                CCl4            0
                CH3CCl3         0
                Halon1211       0
                Halon1301       0
                Halon2402       0
                HCFC22          0
                HCFC141b        0
                HCFC142b        0
                HCFC123         0
                CH3Cl           3684.7
                CH3Br           192/

RHORAD(HALR)    Radiative efficiency of halocarbons                             [W per (m^2*ppt)]
                /CF4            0.00008
                C2F6            0.00026
                C4F10           0.00033
                HFC23           0.00019
                HFC32           0.00011
                HFC4310mee      0.0004
                HFC125          0.00023
                HFC134a         0.00016
                HFC143a         0.00013
                HFC152a         0.00009
                HFC227ea        0.00026
                HFC236fa        0.00028
                HFC245ca        0.00023
                CFC11           0.00025
                CFC12           0.00032
                CFC113          0.00030
                CFC114          0.00031
                CFC115          0.00018
                CCl4            0.00013
                CH3CCl3         0.00006
                Halon1211       0.0003
                Halon1301       0.00032
                Halon2402       0.00033
                HCFC22          0.0002
                HCFC141b        0.00014
                HCFC142b        0.0002
                HCFC123         0.00014/

NHYD(HAL)       Number of hydrogen atoms in halocarbon molecules                [1]
                /CF4            0
                C2F6            0
                C4F10           0
                HFC23           1
                HFC32           2
                HFC4310mee      2
                HFC125          1
                HFC134a         2
                HFC143a         3
                HFC152a         4
                HFC227ea        1
                HFC236fa        2
                HFC245ca        3
                CFC11           0
                CFC12           0
                CFC113          0
                CFC114          0
                CFC115          0
                CCl4            0
                CH3CCl3         3
                Halon1211       0
                Halon1301       0
                Halon2402       0
                HCFC22          1
                HCFC141b        3
                HCFC142b        3
                HCFC123         1
                CH3Cl           3
                CH3Br           3/

NCLO(ODS)       Number of chlorine atoms in molecules of ozone depleting substances [1]
                /CFC11          3
                CFC12           2
                CFC113          3
                CFC114          2
                CFC115          1
                CCl4            4
                CH3CCl3         3
                Halon1211       1
                Halon1301       0
                Halon2402       0
                HCFC22          1
                HCFC141b        2
                HCFC142b        1
                HCFC123         2
                CH3Cl           1
                CH3Br           0/

NBRO(ODS)       Number of bromine atoms in molecules of ozone depleting substances [1]
                /CFC11          0
                CFC12           0
                CFC113          0
                CFC114          0
                CFC115          0
                CCl4            0
                CH3CCl3         0
                Halon1211       1
                Halon1301       1
                Halon2402       2
                HCFC22          0
                HCFC141b        0
                HCFC142b        0
                HCFC123         0
                CH3Cl           0
                CH3Br           1/

RFRF(ODS)       Relative fractional release factors (to that of CFC11)          [1]
                /CFC11          1
                CFC12           0.6
                CFC113          0.75
                CFC114          0.28
                CFC115          0
                CCl4            1.06
                CH3CCl3         1.08
                Halon1211       1.18
                Halon1301       0.62
                Halon2402       1.22
                HCFC22          0.35
                HCFC141b        0.72
                HCFC142b        0.36
                HCFC123         1.11
                CH3Cl           0.80
                CH3Br           1.12/
* Estimate of relative fractional release factor of CFC115 cannot be found in the literature.


*** Radiative forcing --------------------------------------------------------------------------------------------------

RADVOL(T)       Volcanic forcing                                                [W per m^2]


*** Surface temperature change -----------------------------------------------------------------------------------------

RPCTAP(NPCR)    Coefficients in R_PCT of impulse response function              [degree in C]
                 /1              5.91013486
                  2              14.6530283/
* Two sets of parameters are needed here, the first for the computation of PCTRES (parameters add
* to unity), the second for PCSRES (RPCTAP enters the equations through the factor SENSSCALE)


*** Sea level rise (due to thermal expansion of oceans) ----------------------------------------------------------------

RPCSA(NPCR)     Coefficients in R_PCS of impulse response function              [m]
                /1              7.49982018
                 2              0.31360310/

RPCST(NPCR)     Time scales in  R_PCS of impulse response function              [year]
                /1              836.3651570
                 2              30.5081561/


*** DOECLIM ------------------------------------------------------------------------------------------------------------

SQAVG(I)   Convolution part of integrand                                        [1];
