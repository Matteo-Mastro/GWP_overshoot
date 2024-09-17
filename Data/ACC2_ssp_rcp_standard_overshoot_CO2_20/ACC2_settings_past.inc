********** ACC2 version 4.2: ACC2_settings_past.inc ********************************************************************
********** Parameters and equations declarations, assignments, and definitions (past mode) =============================
*
****** Program mode ====================================================================================================
*
* Never change the MODEFLAG assignment below.
* The assignment below (=1) indicates that the program runs on the past mode.
* To run the future mode, execute the file called 'ACC2_main_future.gms'.
*
         MODEFLAG =       1;
*
*
****** Constant declarations ===========================================================================================

SCALARS

TMIX1750PRI     Prior estimate of preindustrial mixed layer temperature                         [degree in C]   /19.59/
* Hoffert et al. (1981) pp.290-291. Also refer to Sundquist and Plummer (1981) p.267

TMIX1750SGM     One-sigma prior uncertainty range of preindustrial mixed layer temperature      [degree in C]   /3/

DMPTMIXPRI      Prior estimate of atmosphere-mixed layer temperature scaling factor             [1]             /0.5/

DMPTMIXSGM      One-sigma prior uncertainty range of atmosphere-mixed layer temperature scaling factor [1]      /0.25/

$if not setglobal BETA $setglobal BETA 0.4
BETAPRI         Prior estimate of beta factor for CO2 fertilization                             [1]             /%BETA%/

BETASGM         One-sigma prior uncertainty range of beta factor for CO2 fertilization          [1]             /0.15/

$if not setglobal Q10 $setglobal Q10 2.00
Q10BASPRI       Prior estimate of Q10 for heterotrophic respiration                             [1]             /%Q10%/
* The prior estimate and its range are based on Tyelker et al. (2001) Table 2.

Q10BASSGM       One-sigma prior uncertainty range of Q10 for heterotrophic respiration          [1]             /0.25/

Q10TMPPRI       Prior estimate of temperature sensitivity to Q10                                [1]             /0.046/
* The prior estimate and its range are based on Tyelker et al. (2001) Figure 2 and p.227.

Q10TMPSGM       One-sigma prior uncertainty range of temperature sensitivity to Q10             [1]             /0.023/

TNPPSENSPRI     Prior estimate of temperature sensitivity to NPP                          [GtC per year per T]  /0.016/
* Based on Meyer et al. (1999) Figure 1.

TNPPSENSSGM     One-sigma prior uncertainty range of temperature sensitivity to NPP       [GtC per year per T]  /0.15/
* Based on the Friedlingstein et al. (2006) Figure 3b.

CO2OCNUPPREPRI  Prior estimate of preindustrial ocean CO2 uptake                                [GtC per year]  /-0.24/
* The estimate of -0.48 was obtained from email communication with Fred Mackenzie (Oct. 26, 2005)
* -0.04 is the estimate calculated from Figure 1 in Ver et al. (1999)
*CO2OCNUPPREPRI  Prior estimate of preindustrial ocean CO2 uptake                                [GtC per year]
*                                                                                                       /-0.3454877590/
* Estimated by full inversion shown in Tanaka and ONeill (2018, Nature Climate Change)

CO2OCNUPPRESGM  One-sigma prior uncertainty range of preindustrial ocean CO2 uptake             [GtC per year]  /0.12/
* The estimate above is according to email communication with Fred Mackenzie (Oct. 26, 2005)

CO2LNDUPPREPRI  Prior estimate of preindustrial land CO2 uptake                                 [GtC per year]  /0.30/
* The estimate of 0.36 - 0.6 was obtained from email communication with Fred Mackenzie (Oct. 26, 2005)
*CO2LNDUPPREPRI  Prior estimate of preindustrial land CO2 uptake                                 [GtC per year]
*                                                                                                        /0.1351753758/
* Estimated by full inversion shown in Tanaka and ONeill (2018, Nature Climate Change)

CO2LNDUPPRESGM  One-sigma prior uncertainty range of preindustrial land CO2 uptake              [GtC per year]  /0.15/
* The estimate of 0.06 is according to email communication with Fred Mackenzie (Oct. 26, 2005)

NATCH4PRI       Prior estimate of natural CH4 emission                          [TgCH4 per year]  /210/
*NATCH4PRI       Prior estimate of natural CH4 emission                          [TgCH4 per year]  /226.3089521219/
* Estimated by inversion to produce stable CH4 concentrations without anthropogenic emissions.
*NATCH4PRI       Prior estimate of natural CH4 emission                          [TgCH4 per year]  /319.9010940839/
* Estimated by full inversion shown in Tanaka and ONeill (2018, Nature Climate Change)

NATCH4SGM       One-sigma prior uncertainty range of natural CH4 emission       [TgCH4 per year]  /120/

NATN2OPRI       Prior estimate of natural N2O emission                          [TgN per year]    /10.2/
*NATN2OPRI       Prior estimate of natural N2O emission                          [TgN per year]    /10.8837305766/
* Estimated by inversion to produce stable N2O concentrations without anthropogenic emissions.
*NATN2OPRI       Prior estimate of natural N2O emission                          [TgN per year]    /11.3419025019/
* Estimated by full inversion shown in Tanaka and ONeill (2018, Nature Climate Change)

NATN2OSGM       One-sigma prior uncertainty range of natural N2O emission       [TgN per year]    /1.2/

TAUCH4OHPRI     Prior estimate of CH4 lifetime with respect to OH depletion                     [year]          /9.6/
*TAUCH4OHPRI     Prior estimate of CH4 lifetime with respect to OH depletion                     [year]  /9.7019949431/
* Estimated by inversion to produce stable CH4 concentrations without anthropogenic emissions.
*TAUCH4OHPRI     Prior estimate of CH4 lifetime with respect to OH depletion                     [year]  /8.5348127733/
* Estimated by full inversion shown in Tanaka and ONeill (2018, Nature Climate Change)

TAUCH4OHSGM     One-sigma prior uncertainty range of CH4 lifetime with respect to OH depletion  [year]          /2.1/

TAUN2OPRI       Prior estimate of N2O lifetime                                  [year]          /110/
*TAUN2OPRI       Prior estimate of N2O lifetime                                  [year]          /117.9826946717/
* Estimated by inversion to produce stable N2O concentrations without anthropogenic emissions.
*TAUN2OPRI       Prior estimate of N2O lifetime                                  [year]          /114.3572563246/
* Estimated by full inversion shown in Tanaka and ONeill (2018, Nature Climate Change)

TAUN2OSGM       One-sigma prior uncertainty range of N2O lifetime               [year]          /13.5/

RAD2XCO2PRI     Prior estimate of CO2 doubling forcing                          [W per m^2]     /3.7/

RAD2XCO2SGM     One-sigma prior uncertainty range of CO2 doubling forcing       [W per m^2]     /0.2/

SCLRADPRI       Prior estimate of forcing scaling factor                        [1]             /1.0/

SCLRADSGM       One-sigma prior uncertainty range of forcing scaling factor     [1]             /0.5/

$if not setglobal ECS $setglobal ECS 3
SENSPRI         Prior estimate of climate sensitivity                           [degree in C]   /%ECS%/

SENSSGM         One-sigma prior uncertainty range of climate sensitivity        [degree in C]   /1.5/

KAPPAPRI        Prior estimate of ocean vertical diffusivity                    [cm2 per second]                 /0.55/
* The prior range of 0 and 5 and best estimate of 0.55 cm^2/s were recommended by Kriegler (personal communication).
* Note that KAPPA is fixed by default at 0.55 cm^2/s.

KAPPASGM        One-sigma prior uncertainty range of ocean vertical diffusivity [cm2 per second]                 /0.275/

RADMISAR1PROPRI    Prior estimate of autoregressive propagator for missing forcing                     [1]       /0.5/
* Based on Kriegler (2005)

RADMISAR1PROSGM    One-sigma prior uncertainty of autoregressive propagator for missing forcing        [1]       /0.25/

T2MAR1PROPRI    Prior estimate of autoregressive propagator for temperature change                     [1]       /0.45/
* Based on Kriegler (2005)

T2MAR1PROSGM    One-sigma prior uncertainty of autoregressive propagator for temperature change        [1]       /0.2/


CO2LUSGMFAC     Factor to expand the one-sigma prior uncertainty range of land use CO2 emission [1]              /2.0/

ANTCH4N2OSGMFAC Factor to expand the one-sigma prior uncertainty range of anthropogenic CH4 and N2O emissions [1] /0.5/

CONCO2SGMFAC    Factor to expand the one-sigma prior uncertainty range of atmospheric CO2 concentration [1]      /4.0/

CONCH4SGMFAC    Factor to expand the one-sigma prior uncertainty range of atmospheric CH4 concentration [1]      /4.0/

CONN2OSGMFAC    Factor to expand the one-sigma prior uncertainty range of atmospheric N2O concentration [1]      /4.0/

VOLSGMFAC       Factor to expand the one-sigma prior uncertainty range of various time series [1]                /4.0/
*               when the volcanic forcing is stronger than -0.5 W per m2

TENDINV         End year for the inverse estimation                                           [Year]
* The value is assigned below.


PARAMETERS

CO2LUPRI(T)     Prior estimates of CO2 emission due to land use change          [GtC per year]

CO2LUSGMORG(T)  One-sigma prior uncertainty ranges for CO2 emission due to land use change [GtC per year]

NINO3(T)        NINO3 to be used in model calculation (including the time lag)  [1]

SOI(T)          SOI to be used in model calculation (including the time lag)    [1];


****** Assignments =====================================================================================================

CO2LUPRI(T) = CO2LUPRIHOU(T);
*
* In the standard inversion, the land use CO2 emission estimate by Houghton (2003) is used.
* One could use alternatively the estimate by Raddatz (in prep)
* CO2LUPRIHOU(T): Houghton (2003)
* CO2LUPRIRAD(T): Raddatz (in prep)

CO2LUSGMORG(T) = CO2LUSGMORGHOU(T);
* In the standard inversion, the land use CO2 emission uncertainty by Houghton (2003) is used.
* One could use alternatively the uncertainty by Raddatz (in prep)
* CO2LUSGMORGHOU(T): Houghton (2003)
* CO2LUSGMORGRAD(T): Raddatz (in prep)

RADVOL(T) = RADVOLAMM(T);
*RADVOL(T) = 0;
*
* In the standard inversion, the volcanic forcing of Ammann et al. (2003) is used among the others.
* RADVOLAMM(T): Ammann et al. (2003)
* Crete, Greenland (Hammer et al., 1980) and Robock and Free (1995)
* RADVOLBRT(T): Bertrand et al. (2002)
* Greenland Ice Sheet Project 2 (GISP2) (Zielinski, 1995)
* RADVOLCRW(T): Crowley (2003)
* Crete, Greenland (Hammer et al., 1980) and Greenland Ice Sheet Project 2 (GISP2) (Zielinski, 1995)
*
* Assignment of the solar forcing is found in the equation section below.

NINO3(T)  = NINO3LAG(T,'7');
* Assigning the length of NINO3 leading the natural change in atmospheric CO2 concentration.

SOI(T)    = SOILAG(T,'4');
* Assigning the length of SOI leading the surface air temperature.

TENDINV   = 251;
* Assigning the end year of the inverse estimation. Note that 1 for TENDINV corresponds to year 1750 and 251 to 2000.
* Use a number above 151. Otherwise the temperature adjustment does not work (see eqs T2M69Q, T2M68Q, T2MINS69Q,
* T2MRCN68Q, T2MPRIQ(T)).


****** Variable declarations ===========================================================================================

VARIABLES

*** Emissions and concentrations ---------------------------------------------------------------------------------------

ANT(T,FOD)               Anthropogenic emissions (both energy and non-energy related emissions)
*                        Natural emissions are fixed and treated as parameters.
*                        Units are gas dependent. Units used in the model follows those in IPCC (2001) Appendix A.
*                        CO2             [GtC per year]
*                        CH4             [TgCH4 per year]
*                        N2O             [TgN per year]
*                        NOx             [TgN per year]
*                        CO              [TgCO per year]
*                        VOC             [TgVOC per year]
*                        SO2             [TgS per year]

CO2FF(T)                 CO2 emissions due to fossil fuel combustion     [GtC per year]

CO2LU(T)                 CO2 emissions due to land use change            [GtC per year]

CON(T,FOF)               Atmospheric concentrations of direct and indirect forcing agents
*                        CO2             [ppm]
*                        CH4             [ppb]
*                        N2O             [ppb]
*                        rOH             [1]
* Concentrations of pollutants (NOx, CO, VOC) are not explicitly specified.
* rOH is the relative concentration (ratio of the current concentration to the concentration in the year 2000).


*** Relationships between natural CO2 concentration change and ENSO ----------------------------------------------------

DCONCO2ANT(T)   Changes in the atmospheric CO2 concentration due to anthropogenic activities            [GtC per year]

DCONCO2NAT(T)   Changes in the atmospheric CO2 concentration due to natural variability                 [GtC per year]

DCONCO2TTL(T)   Changes in the atmospheric CO2 concentration                                            [GtC per year]

YINTRCO2ENSO    Y intercept of the regression line to explain DCONCO2NAT as a function of NINO3         [ppm]

SLOPECO2ENSO    Slope of the regression line to explain DCONCO2NAT as a function of NINO3               [1]

AVENINO3        Average of NINO3                                                                        [K]

AVEDCONCO2NAT   Average of DCONCO2NAT                                                                   [GtC per year]

SX2             Variance of NINO3                                                                       [K^2]

SXY             Covariance of NINO3 and DCONCO2NAT                                                  [(GtC per year)^2]

RSQ1            Square of the correlation coefficient between natural CO2 concentration change and ENSO [1]

RSQ2            Fitness of the regression line to explain natural CO2 concentration change as a function of NINO3 [1]

CONCO2RESWTENSO(T) Atmospheric CO2 concentration residuals without the effect of ENSO                   [ppm]


*** Radiative forcing --------------------------------------------------------------------------------------------------

RADO3(T)         Radiative forcing due to total O3                      [W per m^2]


*** Temperature change -------------------------------------------------------------------------------------------------

T2MPRI(T)        Measurements for T2M                                   [degree in C]

T2M69            Average of T2M between 1961 and 1990                   [degree in C]

T2M68            Average of T2M between 1961 and 1980                   [degree in C]

T2MINS69         Average of T2MINS between 1961 and 1990                [degree in C]

T2MRCN68         Average of T2MRCN between 1961 and 1980                [degree in C]


*** Inverse calculation ------------------------------------------------------------------------------------------------

CO2FFSGMMOD(T)   Modified one-sigma prior uncertainty ranges of fossil fuel CO2 emission       [GtC per year]

CO2LUSGMMOD(T)   Modified one-sigma prior uncertainty ranges of land use CO2 emission          [GtC per year]

ANTSGMMOD(T,FOG) Modified one-sigma prior uncertainty ranges of anthropogenic CH4 and N2O emissions
*                                                                                              [TgCH4 per year]
*                                                                                              [TgN per year]

RADMISSGMMOD(T)  Modified one-sigma prior uncertainty ranges of missing forcing                [W per m^2]


CONCO2SGMMOD(T)  Modified one-sigma prior uncertainty ranges of atmospheric CO2 concentration  [ppm]

CONCH4SGMMOD(T)  Modified one-sigma prior uncertainty ranges of atmospheric CH4 concentration  [ppb]

CONN2OSGMMOD(T)  Modified one-sigma prior uncertainty ranges of atmospheric N2O concentration  [ppb]

T2MSGMMOD(T)     Modified one-sigma prior uncertainty ranges of surface air temperature change [W per m^2]


CO2FFRES(T)      Residuals of fossil fuel CO2 emission                                   [GtC per year]

CO2LURES(T)      Residuals of land use CO2 emission                                      [GtC per year]

ANTRES(T,FOG)    Residuals of anthropogenic CH4 and N2O emissions                        [TgCH4 per year] [TgN per year]

RADMISRES(T)     Residuals of missing forcing                                            [W per m^2]


TMIX1750RES      Residual of preindustrial mixed layer temperature                       [degree in C]

DMPTMIXRES       Residual of atmosphere-mixed layer temperature scaling factor           [1]

BETARES(T)       Residual of beta factor for CO2 fertilization                           [1]

Q10BASRES        Residual of Q10 for heterotrophic respiration                           [1]

Q10TMPRES        Residual of temperature sensitivity to Q10                              [1]

TNPPSENSRES      Residual of temperature sensitivity to NPP                              [1]

CO2OCNUPPRERES(T)        Residual of preindustrial ocean CO2 uptake                      [GtC per year]

CO2LNDUPPRERES(T)        Residual of preindustrial land CO2 uptake                       [GtC per year]

NATCH4RES        Residual of natural CH4 emission                                        [TgCH4 per year]

NATN2ORES        Residual of natural N2O emission                                        [TgN per year]

TAUCH4OHRES      Residual of CH4 lifetime with respect to OH depletion                   [year]

TAUN2ORES        Residual of N2O lifetime                                                [year]

RAD2XCO2RES      Residual of CO2 doubling forcing                                        [W per m^2]

SCLRADRES        Residual of forcing scaling factor                                      [K]

SENSRES(T)       Residual of climate sensitivity                                         [K]

KAPPARES         Residual of ocean vertical diffusivity                                  [cm2 per s]

RADMISAR1PRORES  Residual of forcing autoregressive propagator                           [1]

T2MAR1PRORES     Residual of temperature autoregressive propagator                       [1]


CO2OCNUPANTRES(T)        Residuals of ocean CO2 uptake                                   [GtC per year]

CO2LNDUPANTRES(T)        Residuals of land CO2 uptake                                    [GtC per year]

CONCO2RES(T)     Residuals of atmospheric CO2 concentration                              [ppm]

CONCH4RES(T)     Residuals of atmospheric CH4 concentration                              [ppb]

CONN2ORES(T)     Residuals of atmospheric N2O concentration                              [ppb]

T2MRES(T)        Residuals of surface air temperature change                             [degree in C]


CO2FFWSQRES(T)   Squared weighted residuals of CO2FF                                     [1]

CO2LUWSQRES(T)   Squared weighted residuals of CO2LU                                     [1]

ANTWSQRES(T,FOG) Squared weighted residuals of ANT                                       [1]

RADMISWSQRES(T)  Squared weighted residuals of RADMIS                                    [1]


CO2OCNUPANTWSQRES(T)     Squared weighted residuals of CO2OCNUPANT                       [1]

CO2LNDUPANTWSQRES(T)     Squared weighted residuals of CO2LNDUPANT                       [1]

CONCO2WSQRES(T)  Squared weighted residuals of CONCO2                                    [1]

CONCH4WSQRES(T)  Squared weighted residuals of CONCH4                                    [1]

CONN2OWSQRES(T)  Squared weighted residuals of CONN2O                                    [1]

T2MWSQRES(T)     Squared weighted residuals of T2M                                       [1]


CO2FFRESSUM      Summation of squared weighted residuals of CO2FF                        [1]

CO2LURESSUM      Summation of squared weighted residuals of CO2LU                        [1]

ANTRESSUMCH4     Summation of squared weighted residuals of ANTCH4                       [1]

ANTRESSUMN2O     Summation of squared weighted residuals of ANTN2O                       [1]

RADMISRESSUM     Summation of squared weighted residuals of RADMIS                       [1]


TMIX1750RESSUM   Summation of squared weighted residuals of TMIX1750                     [1]

DMPTMIXRESSUM    Summation of squared weighted residuals of DMPTMIX                      [1]

BETARESSUM       Summation of squared weighted residuals of BETA                         [1]

Q10BASRESSUM     Summation of squared weighted residuals of Q10BAS                       [1]

Q10TMPRESSUM     Summation of squared weighted residuals of Q10TMP                       [1]

TNPPSENSRESSUM   Summation of squared weighted residuals of TNPPSENS                     [1]

CO2OCNUPPRERESSUM        Summation of squared weighted residuals of CO2OCNUPPRE          [1]

CO2LNDUPPRERESSUM        Summation of squared weighted residuals of CO2LNDUPPRE          [1]

NATCH4RESSUM     Summation of squared weighted residuals of NATCH4                       [1]

NATN2ORESSUM     Summation of squared weighted residuals of NATN2O                       [1]

TAUCH4OHRESSUM   Summation of squared weighted residuals of TAUCH4OH                     [1]

TAUN2ORESSUM     Summation of squared weighted residuals of TAUN2O                       [1]

RAD2XCO2RESSUM   Summation of squared weighted residuals of RAD2XCO2                     [1]

SCLRADRESSUM     Summation of squared weighted residuals of SCLRAD                       [1]

SENSRESSUM       Summation of squared weighted residuals of SENS                         [1]

KAPPARESSUM      Summation of squared weighted residuals of KAPPA                        [1]

RADMISAR1PRORESSUM Summation of squared weighted residuals of RADMISAR1PRO               [1]

T2MAR1PRORESSUM  Summation of squared weighted residuals of T2MAR1PRO                    [1]


CO2OCNUPANTRESSUM        Summation of squared weighted residuals of CO2OCNUPANT          [1]

CO2LNDUPANTRESSUM        Summation of squared weighted residuals of CO2LNDUPANT          [1]

CONCO2RESSUM     Summation of squared weighted residuals of CONCO2                       [1]

CONCH4RESSUM     Summation of squared weighted residuals of CONCH4                       [1]

CONN2ORESSUM     Summation of squared weighted residuals of CONN2O                       [1]

T2MRESSUM        Summation of squared weighted residuals of T2M                          [1];


****** Equation declarations ===========================================================================================

EQUATIONS

COC1ATMBACQ             Assignment of the preindustrial CO2 load in the atmosphere

COC1MIXBACQ1            Assignment of the preindustrial CO2 load in the mixed layer

COC1MIXBACQ2            Assignment of the preindustrial CO2 load in the mixed layer

CTRBACQ(NTR)            Equilibrium condition for the land biosphere

CTR1750Q(NTR)           Assignment of the preindustrial CO2 load in the land biosphere


ANTCO2Q(T)              Assignment of the anthropogenic CO2 emissions

ANTASGNQ1(T,FOG)        Assignment of the anthropogenic emissions of CH4 and N2O

ANTASGNQ2(T,FOE)        Assignment of the anthropogenic emissions of NOx CO VOC and SO2

STDSTTCH4Q(T)           Steady state condition for the preindustrial CH4 concentration

STDSTTN2OQ(T)           Steady state condition for the preindustrial N2O concentration


DCONCO2ANTQA(T)         Changes in the atmospheric CO2 concentration due to anthropogenic activities

DCONCO2ANTQB(T)         Changes in the atmospheric CO2 concentration due to anthropogenic activities

DCONCO2ANTQC(T)         Changes in the atmospheric CO2 concentration due to anthropogenic activities

DCONCO2NATQA(T)         Changes in the atmospheric CO2 concentration due to natural variability

DCONCO2NATQB(T)         Changes in the atmospheric CO2 concentration due to natural variability

DCONCO2NATQC(T)         Changes in the atmospheric CO2 concentration due to natural variability

DCONCO2TTLQA(T)         Changes in the atmospheric CO2 concentration

DCONCO2TTLQB(T)         Changes in the atmospheric CO2 concentration

DCONCO2TTLQC(T)         Changes in the atmospheric CO2 concentration

DCONCO2ENSOQ1(T)        Changes in the atmospheric CO2 concentration due to ENSO

DCONCO2ENSOQ2(T)        Changes in the atmospheric CO2 concentration due to ENSO

DCONCO2ENSOQ3(T)        Changes in the atmospheric CO2 concentration due to ENSO

DCONCO2ENSOQ4(T)        Changes in the atmospheric CO2 concentration due to ENSO

AVENINO3Q               Average of NINO3

AVEDCONCO2NATQ          Average of DCONCO2NAT

SX2Q1                   Calculation of SX2

SX2Q2                   Null assignment of SX2

SXYQ1                   Calculation of SXY

SXYQ2                   Null assignment of SXY

SLOPECO2ENSOQ1          Calculation of SLOPECO2ENSO

SLOPECO2ENSOQ2          Null assignment of SLOPECO2ENSO

YINTRCO2ENSOQ1          Calculation of YINTRCO2ENSO

YINTRCO2ENSOQ2          Null assignment of YINTRCO2ENSO

RSQ11Q                  Correlation coefficient of natural CO2 concentration change and NINO3

RSQ12Q                  Null assignment for RSQ1

RSQ21Q                  Fitness of the computed regression lines to explain natural CO2 concentration change
*                       as a function of NINO3

RSQ22Q                  Null assignment for RSQ2

CONCO2RESWTENSOQ(T)     Residuals for the atmospheric CO2 concentration without the influence from ENSO


RADHALQ(T,HALR)         Calculation of radiative forcing due to all the halocarbons

RADSF6Q(T)              Calculation of radiative forcing due to SF6

RADO3Q(T)               Calculation of radiative forcing due to total O3

RADSOLQ(T)              Assignment of solar forcing

RADSOL2000Q             Calculation of solar forcing in 2000

RADSOL20AVEQ            Calculation of solar forcing averaged between 1901 and 2000

RADAERQ(T)              Calculation of total aerosol forcing

RADANTQ(T)              Radiative forcing due to all the anthropogenic factors

RADANTLQ(T)             Radiative forcing due to all the anthropogenic factors over the land used for DOECLIM

RADANTOQ(T)             Radiative forcing due to all the anthropogenic factors over the ocean used for DOECLIM

RADTTLPRIQ(T)           Calculation of total radiative forcing including the prior missing forcing

RADTTLQ(T)              Calculation of total radiative forcing including the posterior missing forcing

RADAERAVEQ              Calculation of total aerosol forcing averaged between 1950 and 2000

RADMISAVEQ              Calculation of total missing forcing averaged between 1950 and 2000


T2M1750Q                Assignment of T2M1750

T2M1990Q                Assignment of T2M1990

TAL1750Q                Assignment of TAL1750

TSS1750Q                Assignment of TSS1750

T2M69Q                  Average of T2M between 1960 and 1990

T2M68Q                  Average of T2M between 1960 and 1980

T2MINS69Q               Average of T2MINS between 1960 and 1990

T2MRCN68Q               Average of T2MRCN between 1960 and 1980

T2MPRIQ(T)              Prior estimate of temperature change


DT2MENSOQ1(T)           Assignment for DT2MENSO

DT2MENSOQ2(T)           Assignment for DT2MENSO


CO2FFQ1(T)              Equality for fossil fuel CO2 emission (uncertainty taken into account)

CO2FFQ2(T)              Equality for fossil fuel CO2 emission (uncertainty not taken into account)

CO2FFQ3(T)              Equality for fossil fuel CO2 emission (uncertainty not taken into account)

CO2LUQ1(T)              Equality for land use CO2 emission (uncertainty taken into account)

CO2LUQ2(T)              Equality for land use CO2 emission (uncertainty not taken into account)

CO2LUQ3(T)              Equality for land use CO2 emission (uncertainty not taken into account)

ANTSPQ1(T,FOG)          Equality for anthropogenic CH4 and N2O emission (uncertainty taken into account)

ANTSPQ2(T,FOG)          Equality for anthropogenic CH4 and N2O emission (uncertainty not taken into account)

ANTSPQ3(T,FOG)          Equality for anthropogenic CH4 and N2O emission (uncertainty not taken into account)

RADMISSPINQ1(T)         Equality for missing forcing (uncertainty and autocorrelations taken into account)

RADMISSPINQ2(T)         Equality for missing forcing (uncertainty and autocorrelations taken into account)

RADMISSPINQ3(T)         Equality for missing forcing (uncertainty taken into account)

RADMISSPINQ4(T)         Equality for missing forcing (uncertainty not taken into account)

RADMISSPINQ5(T)         Equality for missing forcing (uncertainty not taken into account)


TMIX1750Q1              Equality for preindustrial mixed layer temperature (uncertainty taken into account)

TMIX1750Q2              Equality for preindustrial mixed layer temperature (uncertainty not taken into account)

TMIX1750Q3              Equality for preindustrial mixed layer temperature (uncertainty not taken into account)

DMPTMIXQ1               Equality for atmosphere-mixed layer temperature scaling factor (uncertainty taken into account)

DMPTMIXQ2               Equality for atmosphere-mixed layer temperature scaling factor
*                       (uncertainty not taken into account)

DMPTMIXQ3               Equality for atmosphere-mixed layer temperature scaling factor
*                       (uncertainty not taken into account)

BETAQ1(T)               Equality for beta factor for CO2 fertilization (uncertainty taken into account)

BETAQ2(T)               Equality for beta factor for CO2 fertilization (uncertainty not taken into account)

BETAQ3(T)               Equality for beta factor for CO2 fertilization (uncertainty not taken into account)

Q10BASQ1                Equality for Q10 for heterotrophic respiration  (uncertainty taken into account)

Q10BASQ2                Equality for Q10 for heterotrophic respiration  (uncertainty not taken into account)

Q10BASQ3                Equality for Q10 for heterotrophic respiration  (uncertainty not taken into account)

Q10TMPQ1                Equality for temperature sensitivity to Q10 (uncertainty taken into account)

Q10TMPQ2                Equality for temperature sensitivity to Q10 (uncertainty not taken into account)

Q10TMPQ3                Equality for temperature sensitivity to Q10 (uncertainty not taken into account)

TNPPSENSQ1              Equality for temperature sensitivity to NPP (uncertainty taken into account)

TNPPSENSQ2              Equality for temperature sensitivity to NPP (uncertainty not taken into account)

TNPPSENSQ3              Equality for temperature sensitivity to NPP (uncertainty not taken into account)

CO2OCNUPPREQ1(T)        Equality for preindustrial ocean CO2 uptake (uncertainty taken into account)

CO2OCNUPPREQ2(T)        Equality for preindustrial ocean CO2 uptake (uncertainty not taken into account)

CO2OCNUPPREQ3(T)        Equality for preindustrial ocean CO2 uptake (uncertainty not taken into account)

CO2LNDUPPREQ1(T)        Equality for preindustrial land CO2 uptake (uncertainty taken into account)

CO2LNDUPPREQ2(T)        Equality for preindustrial land CO2 uptake (uncertainty not taken into account)

CO2LNDUPPREQ3(T)        Equality for preindustrial land CO2 uptake (uncertainty not taken into account)

NATCH4Q1                Equality for natural CH4 emission (uncertainty taken into account)

NATCH4Q2                Equality for natural CH4 emission (uncertainty not taken into account)

NATCH4Q3                Equality for natural CH4 emission (uncertainty not taken into account)

NATN2OQ1                Equality for natural N2O emission (uncertainty taken into account)

NATN2OQ2                Equality for natural N2O emission (uncertainty not taken into account)

NATN2OQ3                Equality for natural N2O emission (uncertainty not taken into account)

TAUCH4OHQ1              Equality for CH4 lifetime with respect to OH depletion (uncertainty taken into account)

TAUCH4OHQ2              Equality for CH4 lifetime with respect to OH depletion (uncertainty not taken into account)

TAUCH4OHQ3              Equality for CH4 lifetime with respect to OH depletion (uncertainty not taken into account)

TAUN2OQ1                Equality for N2O lifetime (uncertainty taken into account)

TAUN2OQ2                Equality for N2O lifetime (uncertainty not taken into account)

TAUN2OQ3                Equality for N2O lifetime (uncertainty not taken into account)

RAD2XCO2Q1              Equality for CO2 doubling forcing (uncertainty taken into account)

RAD2XCO2Q2              Equality for CO2 doubling forcing (uncertainty not taken into account)

RAD2XCO2Q3              Equality for CO2 doubling forcing (uncertainty not taken into account)

SCLRADQ1                Equality for forcing scaling factor (uncertainty taken into account)

SCLRADQ2                Equality for forcing scaling factor (uncertainty not taken into account)

SCLRADQ3                Equality for forcing scaling factor (uncertainty not taken into account)

SENSQ1(T)               Equality for climate sensitivity (uncertainty taken into account)

SENSQ2(T)               Equality for climate sensitivity (uncertainty not taken into account)

SENSQ3(T)               Equality for climate sensitivity (uncertainty not taken into account)

KAPPAQ1                 Equality for ocean vertical diffusivity (uncertainty taken into account)

KAPPAQ2                 Equality for ocean vertical diffusivity (uncertainty not taken into account)

KAPPAQ3                 Equality for ocean vertical diffusivity (uncertainty not taken into account)

RADMISAR1PROQ1          Equality for forcing autoregressive propagator (uncertainty taken into account)

RADMISAR1PROQ2          Equality for forcing autoregressive propagator (uncertainty not taken into account)

RADMISAR1PROQ3          Equality for forcing autoregressive propagator (uncertainty not taken into account)

T2MAR1PROQ1             Equality for temperature autoregressive propagator (uncertainty taken into account)

T2MAR1PROQ2             Equality for temperature autoregressive propagator (uncertainty not taken into account)

T2MAR1PROQ3             Equality for temperature autoregressive propagator (uncertainty not taken into account)


CO2OCNUPANTSPQ1(T)      Equality for ocean CO2 uptake (uncertainty taken into account)

CO2OCNUPANTSPQ3(T)      Equality for ocean CO2 uptake (uncertainty not taken into account)

CO2LNDUPANTSPQ1(T)      Equality for land CO2 uptake (uncertainty taken into account)

CO2LNDUPANTSPQ3(T)      Equality for land CO2 uptake (uncertainty not taken into account)

CONCO2SPINQ1(T)         Equality for atmospheric CO2 concentration (uncertainty taken into account)

CONCO2SPINQ3(T)         Equality for atmospheric CO2 concentration (uncertainty not taken into account)

CONCH4SPINQ1(T)         Equality for atmospheric CH4 concentration (uncertainty taken into account)

CONCH4SPINQ3(T)         Equality for atmospheric CH4 concentration (uncertainty not taken into account)

CONN2OSPINQ1(T)         Equality for atmospheric N2O concentration (uncertainty taken into account)

CONN2OSPINQ3(T)         Equality for atmospheric N2O concentration (uncertainty not taken into account)

T2MSPQ1(T)              Equality for surface air temperature change
*                       (uncertainty and autocorrelations taken into account)

T2MSPQ2(T)              Equality for surface air temperature change
*                       (uncertainty and autocorrelations taken into account)

T2MSPQ3(T)              Equality for surface air temperature change (uncertainty taken into account)

T2MSPQ4(T)              Equality for surface air temperature change (uncertainty not taken into account)


CO2FFSGMMODQ(T)         Modifications of the one-sigma prior uncertainty ranges of fossil fuel CO2 emission

CO2LUSGMMODQ1(T)        Modifications of the one-sigma prior uncertainty ranges of land use CO2 emission
* (with a larger uncertainty range of land use CO2 emission)

CO2LUSGMMODQ2(T)        Modifications of the one-sigma prior uncertainty ranges of land use CO2 emission
* (without a larger uncertainty range of land use CO2 emission)

ANTCH4SGMMOD1(T)        Modifications of the one-sigma prior uncertainty ranges of anthropogenic CH4 emission
* (with larger uncertainty ranges than estimated in the original literature)

ANTCH4SGMMOD2(T)        Modifications of the one-sigma prior uncertainty ranges of anthropogenic CH4 emission
* (without larger uncertainty ranges than estimated in the original literature)

ANTN2OSGMMOD1(T)        Modifications of the one-sigma prior uncertainty ranges of anthropogenic N2O emission
* (with larger uncertainty ranges than estimated in the original literature)

ANTN2OSGMMOD2(T)        Modifications of the one-sigma prior uncertainty ranges of anthropogenic N2O emission
* (without larger uncertainty ranges than estimated in the original literature)

RADMISSGMMODQ1(T)       Modifications of the one-sigma prior uncertainty ranges of missing forcing
* (extending the uncertainty range during the periods of large volcanic eruptions)

RADMISSGMMODQ2(T)       Modifications of the one-sigma prior uncertainty ranges of missing forcing
* (not extending the uncertainty range during the periods of large volcanic eruptions)


CONCO2SGMMODQ1(T)       Modifications of the one-sigma prior uncertainty ranges of atmospheric CO2 concentration
* (with a larger uncertainty range of CO2 concentration and larger uncertainty ranges during large volcanic eruptions)

CONCO2SGMMODQ2(T)       Modifications of the one-sigma prior uncertainty ranges of atmospheric CO2 concentration
* (with a larger uncertainty range of CO2 concentration)

CONCO2SGMMODQ3(T)       Modifications of the one-sigma prior uncertainty ranges of atmospheric CO2 concentration
* (without a larger uncertainty range of CO2 concentration
*  and larger uncertainty ranges during large volcanic eruptions)

CONCO2SGMMODQ4(T)       Modifications of the one-sigma prior uncertainty ranges of atmospheric CO2 concentration
* (without a larger uncertainty range of CO2 concentration)

CONCH4SGMMODQ1(T)       Modifications of the one-sigma prior uncertainty ranges of atmospheric CH4 concentration
* (with a larger uncertainty range of CH4 concentration and larger uncertainty ranges during large volcanic eruptions)

CONCH4SGMMODQ2(T)       Modifications of the one-sigma prior uncertainty ranges of atmospheric CH4 concentration
* (with a larger uncertainty range of CH4 concentration)

CONCH4SGMMODQ3(T)       Modifications of the one-sigma prior uncertainty ranges of atmospheric CH4 concentration
* (without a larger uncertainty range of CH4 concentration
*  and larger uncertainty ranges during large volcanic eruptions)

CONCH4SGMMODQ4(T)       Modifications of the one-sigma prior uncertainty ranges of atmospheric CH4 concentration
* (without a larger uncertainty range of CH4 concentration)

CONN2OSGMMODQ1(T)       Modifications of the one-sigma prior uncertainty ranges of atmospheric N2O concentration
* (with a larger uncertainty range of N2O concentration and larger uncertainty ranges during large volcanic eruptions)

CONN2OSGMMODQ2(T)       Modifications of the one-sigma prior uncertainty ranges of atmospheric N2O concentration
* (with a larger uncertainty range of N2O concentration)

CONN2OSGMMODQ3(T)       Modifications of the one-sigma prior uncertainty ranges of atmospheric N2O concentration
* (without a larger uncertainty range of N2O concentration
*  and larger uncertainty ranges during large volcanic eruptions)

CONN2OSGMMODQ4(T)       Modifications of the one-sigma prior uncertainty ranges of atmospheric N2O concentration
* (without a larger uncertainty range of N2O concentration)

T2MSGMMODQ1(T)          Modifications of the one-sigma prior uncertainty ranges of surface air temperature change
* (extending the uncertainty range during the periods with large volcanic eruptions)
* (with AR1 model for temperature change)

T2MSGMMODQ2(T)          Modifications of the one-sigma prior uncertainty ranges of surface air temperature change
* (not extending the uncertainty range during the periods with large volcanic eruptions)
* (with AR1 model for temperature change)

T2MSGMMODQ3(T)          Modifications of the one-sigma prior uncertainty ranges of surface air temperature change
* (extending the uncertainty range during the periods with large volcanic eruptions)
* (without AR1 model for temperature change)

T2MSGMMODQ4(T)          Modifications of the one-sigma prior uncertainty ranges of surface air temperature change
* (not extending the uncertainty range during the periods with large volcanic eruptions)
* (without AR1 model for temperature change)


CO2FFWSQRESQ1(T)        Squared weighted residuals of fossil fuel CO2 emission

CO2FFWSQRESQ2(T)        Squared weighted residuals of fossil fuel CO2 emission

CO2LUWSQRESQ1(T)        Squared weighted residuals of land use CO2 emission

CO2LUWSQRESQ2(T)        Squared weighted residuals of land use CO2 emission

ANTWSQRESQ1(T,FOG)      Squared weighted residuals of anthropogenic CH4 and N2O emission

ANTWSQRESQ2(T,FOG)      Squared weighted residuals of anthropogenic CH4 and N2O emission

RADMISWSQRESQ1(T)       Squared weighted residuals of missing forcing

RADMISWSQRESQ2(T)       Squared weighted residuals of missing forcing


CO2OCNUPANTWSQRESQ1(T)  Squared weighted residuals of ocean CO2 uptake

CO2OCNUPANTWSQRESQ2(T)  Squared weighted residuals of ocean CO2 uptake

CO2LNDUPANTWSQRESQ1(T)  Squared weighted residuals of land CO2 uptake

CO2LNDUPANTWSQRESQ2(T)  Squared weighted residuals of land CO2 uptake

CONCO2WSQRESQ1(T)       Squared weighted residuals of atmospheric CO2 concentration

CONCO2WSQRESQ2(T)       Squared weighted residuals of atmospheric CO2 concentration

CONCH4WSQRESQ1(T)       Squared weighted residuals of atmospheric CH4 concentration

CONCH4WSQRESQ2(T)       Squared weighted residuals of atmospheric CH4 concentration

CONN2OWSQRESQ1(T)       Squared weighted residuals of atmospheric N2O concentration

CONN2OWSQRESQ2(T)       Squared weighted residuals of atmospheric N2O concentration

T2MWSQRESQ1(T)          Squared weighted residuals of surface air temperature change

T2MWSQRESQ2(T)          Squared weighted residuals of surface air temperature change


CO2FFRESSUMQ1           Summation of squared weighted residuals of fossil fuel CO2 emission

CO2FFRESSUMQ2           Summation of squared weighted residuals of fossil fuel CO2 emission

CO2LURESSUMQ1           Summation of squared weighted residuals of land use CO2 emission

CO2LURESSUMQ2           Summation of squared weighted residuals of land use CO2 emission

ANTRESSUMCH4Q1          Summation of squared weighted residuals of anthropogenic CH4 emission

ANTRESSUMCH4Q2          Summation of squared weighted residuals of anthropogenic CH4 emission

ANTRESSUMN2OQ1          Summation of squared weighted residuals of anthropogenic N2O emission

ANTRESSUMN2OQ2          Summation of squared weighted residuals of anthropogenic N2O emission

RADMISRESSUMQ1          Summation of squared weighted residuals of missing forcing

RADMISRESSUMQ2          Summation of squared weighted residuals of missing forcing


TMIX1750RESSUMQ1        Summation of squared weighted residuals of preindustrial mixed layer temperature

TMIX1750RESSUMQ2        Summation of squared weighted residuals of preindustrial mixed layer temperature

DMPTMIXRESSUMQ1         Summation of squared weighted residuals of atmosphere-mixed layer temperature scaling factor

DMPTMIXRESSUMQ2         Summation of squared weighted residuals of atmosphere-mixed layer temperature scaling factor

BETARESSUMQ1            Summation of squared weighted residuals of beta factor for CO2 fertilization

BETARESSUMQ2            Summation of squared weighted residuals of beta factor for CO2 fertilization

Q10BASRESSUMQ1          Summation of squared weighted residuals of Q10 for heterotrophic respiration

Q10BASRESSUMQ2          Summation of squared weighted residuals of Q10 for heterotrophic respiration

Q10TMPRESSUMQ1          Summation of squared weighted residuals of temperature sensitivity to Q10

Q10TMPRESSUMQ2          Summation of squared weighted residuals of temperature sensitivity to Q10

TNPPSENSRESSUMQ1        Summation of squared weighted residuals of temperature sensitivity to NPP

TNPPSENSRESSUMQ2        Summation of squared weighted residuals of temperature sensitivity to NPP

CO2OCNUPPRERESSUMQ1     Summation of squared weighted residuals of preindustrial ocean CO2 uptake

CO2OCNUPPRERESSUMQ2     Summation of squared weighted residuals of preindustrial ocean CO2 uptake

CO2LNDUPPRERESSUMQ1     Summation of squared weighted residuals of preindustrial land CO2 uptake

CO2LNDUPPRERESSUMQ2     Summation of squared weighted residuals of preindustrial land CO2 uptake

NATCH4RESSUMQ1          Summation of squared weighted residuals of natural CH4 emission

NATCH4RESSUMQ2          Summation of squared weighted residuals of natural CH4 emission

NATN2ORESSUMQ1          Summation of squared weighted residuals of natural N2O emission

NATN2ORESSUMQ2          Summation of squared weighted residuals of natural N2O emission

TAUCH4OHRESSUMQ1        Summation of squared weighted residuals of CH4 lifetime with respect to OH depletion

TAUCH4OHRESSUMQ2        Summation of squared weighted residuals of CH4 lifetime with respect to OH depletion

TAUN2ORESSUMQ1          Summation of squared weighted residuals of N2O lifetime

TAUN2ORESSUMQ2          Summation of squared weighted residuals of N2O lifetime

RAD2XCO2RESSUMQ1        Summation of squared weighted residuals of CO2 doubling forcing

RAD2XCO2RESSUMQ2        Summation of squared weighted residuals of CO2 doubling forcing

SCLRADRESSUMQ1          Summation of squared weighted residuals of forcing scaling factor

SCLRADRESSUMQ2          Summation of squared weighted residuals of forcing scaling factor

SENSRESSUMQ1            Summation of squared weighted residuals of climate sensitivity

SENSRESSUMQ2            Summation of squared weighted residuals of climate sensitivity

KAPPARESSUMQ1           Summation of squared weighted residuals of ocean vertical diffusivity

KAPPARESSUMQ2           Summation of squared weighted residuals of ocean vertical diffusivity

RADMISAR1PRORESSUMQ1    Summation of squared weighted residuals of forcing autoregressive propagator

RADMISAR1PRORESSUMQ2    Summation of squared weighted residuals of forcing autoregressive propagator

T2MAR1PRORESSUMQ1       Summation of squared weighted residuals of temperature autoregressive propagator

T2MAR1PRORESSUMQ2       Summation of squared weighted residuals of temperature autoregressive propagator


CO2OCNUPANTRESSUMQ1     Summation of squared weighted residuals of ocean CO2 uptake

CO2OCNUPANTRESSUMQ2     Summation of squared weighted residuals of ocean CO2 uptake

CO2LNDUPANTRESSUMQ1     Summation of squared weighted residuals of land CO2 uptake

CO2LNDUPANTRESSUMQ2     Summation of squared weighted residuals of land CO2 uptake

CONCO2RESSUMQ1          Summation of squared weighted residuals of atmospheric CO2 concentration

CONCO2RESSUMQ2          Summation of squared weighted residuals of atmospheric CO2 concentration

CONCH4RESSUMQ1          Summation of squared weighted residuals of atmospheric CH4 concentration

CONCH4RESSUMQ2          Summation of squared weighted residuals of atmospheric CH4 concentration

CONN2ORESSUMQ1          Summation of squared weighted residuals of atmospheric N2O concentration

CONN2ORESSUMQ2          Summation of squared weighted residuals of atmospheric N2O concentration

T2MRESSUMQ1             Summation of squared weighted residuals of surface air temperature change

T2MRESSUMQ2             Summation of squared weighted residuals of surface air temperature change;


****** Equation Definitions ============================================================================================

*** Initial setting for carbon cycle -----------------------------------------------------------------------------------

COC1ATMBACQ                       .. COC1ATMBAC =E= CONCO2INT('1750')/CNVCO2;

COC1MIXBACQ1 $(DYNCARBSWCH > 0.5) .. COC1MIXBAC =E= CONCO2INT('1750')/ALPHAS('1750')*(1+AK1S('1750')/H('1750')
                                                    +AK1S('1750')*AK2S('1750')/H('1750')**2)*ALPHAC;
* Special case of equation COC1MIXQ(T) in ACC2_equations_common.inc

COC1MIXBACQ2 $(DYNCARBSWCH < 0.5) .. COC1MIXBAC =E= 504.025863;
* According to SIAM input used for ICM1.1

CTRBACQ(NTR)                      .. NPPBAC('1750',NTR) =E= HTRRES('1750',NTR);

CTR1750Q(NTR)                     .. CTR1750(NTR) =E= CTR('1750',NTR);


*** Emissions of non-CO2 gases -----------------------------------------------------------------------------------------

ANTCO2Q(T)                       .. ANT(T,'CO2')  =E= CO2FF(T)+CO2LU(T);
*ANTCO2Q(T)                       .. ANT(T,'CO2')  =E= 0;

ANTASGNQ1(T,FOG)                 .. ANTPRI(T,FOG) =E= ANTSCN(T,FOG);
*ANTASGNQ1(T,FOG)                 .. ANTPRI(T,FOG) =E= 0;

ANTASGNQ2(T,FOE) $(Ord(FOE) > 2) .. ANT(T,FOE)    =E= ANTSCN(T,FOE);
*ANTASGNQ2(T,FOE) $(Ord(FOE) > 2) .. ANT(T,FOE)    =E= 0;

* The following two equations provide steady state conditions for CH4 and N2O,
* ensuring that CH4 and N2O concentrations do not change in 1750.
* The CH4 and N2O concentration measurements in a first few time steps are nearly flat,
* providing also quasi-steady state conditions.

STDSTTCH4Q(T) $((NATCH4N2OSWCH > 0.5) or (TAUCH4N2OSWCH > 0.5))
                       .. 1/CNVCH4*NATCH4 =E= (1/(1/(CON('1750','rOH')/TAUCH4OH+1/TAUCH4SS)))*CON('1750','CH4');

STDSTTN2OQ(T) $((NATCH4N2OSWCH > 0.5) or (TAUCH4N2OSWCH > 0.5))
                       .. 1/CNVN2O*NATN2O =E= 1/TAUN2O*(CON('1750','N2O')/CONN2O2000R)**(-SENTAUN2O)*CON('1750','N2O');


*** Natural CO2 concentration change and ENSO --------------------------------------------------------------------------
*
* Below a Model I regression line is calculated and used to relate the natural CO2 concentration change to ENSO.
* A Model II regression line would be precise statistically, but it cannot be calculated
* because the nonlinear program (NLP) in GAMS does not accommodate ranking operations required
* for computing Model II regression lines. Nevertheless, choice of a regression model does not significantly alter
* the inversion results. For the discussion of the Model I and II regression lines, see Laws (1997), for example.

DCONCO2ANTQA(T) $(Ord(T) = 1)                      .. DCONCO2ANT(T) =E= 0;

DCONCO2ANTQB(T) $((Ord(T) > 1) and (Ord(T) < 251)) .. DCONCO2ANT(T) =E= ANT(T,'CO2')-CO2OCNUPANT(T)-CO2LNDUPANT(T);

DCONCO2ANTQC(T) $(Ord(T) = 251)                    .. DCONCO2ANT(T) =E= 0;

DCONCO2NATQA(T) $(Ord(T) = 1)                      .. DCONCO2NAT(T) =E= 0;

DCONCO2NATQB(T) $((Ord(T) > 1) and (Ord(T) < 251)) .. DCONCO2NAT(T) =E= DCONCO2TTL(T)-DCONCO2ANT(T);
* Preindustrial ocean CO2 degassing and the preindustrial land CO2 uptake are included in DCONCO2NAT.

DCONCO2NATQC(T) $(Ord(T) = 251)                    .. DCONCO2NAT(T) =E= 0;

DCONCO2TTLQA(T) $(Ord(T) = 1)                      .. DCONCO2TTL(T) =E= 0;

DCONCO2TTLQB(T) $((Ord(T) > 1) and (Ord(T) < 251)) .. DCONCO2TTL(T) =E= ((CONCO2PRI(T+1)-CONCO2PRI(T))
                                                                       +(CONCO2PRI(T)-CONCO2PRI(T-1)))/2/DELTAT/CNVCO2;

DCONCO2TTLQC(T) $(Ord(T) = 251)                    .. DCONCO2TTL(T) =E= 0;
* DCONCO2TTL has to be defined using the CO2 concentration measurement data, which contain the natural variability.

DCONCO2ENSOQ1(T) $(ENSOCO2SWCH < 0.5)                       .. DCONCO2ENSO(T) =E= 0;

DCONCO2ENSOQ2(T) $((ENSOCO2SWCH > 0.5) and (Ord(T) <  211)) .. DCONCO2ENSO(T) =E= 0;

DCONCO2ENSOQ3(T) $((ENSOCO2SWCH > 0.5) and (TENDINV < 196)) .. DCONCO2ENSO(T) =E= 0;

DCONCO2ENSOQ4(T) $((ENSOCO2SWCH > 0.5) and (Ord(T) >= 211)) .. DCONCO2ENSO(T) =E= YINTRCO2ENSO+SLOPECO2ENSO*NINO3(T);
* Strictly speaking, DCONCO2ENSO('2000') is an extrapolation from the regression line between 1960 and 1999.

AVENINO3Q $(ENSOCO2SWCH > 0.5)
                 .. AVENINO3 =E= SUM(T$((Ord(T) >= 211) and (Ord(T) < 251) and (RADVOL(T) > -0.5)),NINO3(T))
                                 /SUM(T$((Ord(T) >= 211) and (Ord(T) < 251) and (RADVOL(T) > -0.5)),1);

AVEDCONCO2NATQ $(ENSOCO2SWCH > 0.5)
                 .. AVEDCONCO2NAT =E= SUM(T$((Ord(T) >= 211) and (Ord(T) < 251) and (RADVOL(T) > -0.5)),DCONCO2NAT(T))
                                      /SUM(T$((Ord(T) >= 211) and (Ord(T) < 251) and (RADVOL(T) > -0.5)),1);

SX2Q1 $(ENSOCO2SWCH > 0.5)
                 .. SX2 =E= SUM(T$((Ord(T) >= 211) and (Ord(T) < 251) and (RADVOL(T) > -0.5)),
                            (NINO3(T)-AVENINO3)*(NINO3(T)-AVENINO3))
                            /(SUM(T$((Ord(T) >= 211) and (Ord(T) < 251) and (RADVOL(T) > -0.5)),1)-1);

SX2Q2 $(ENSOCO2SWCH < 0.5)
                 .. SX2 =E= 0.001;

SXYQ1 $(ENSOCO2SWCH > 0.5)
                 .. SXY =E= SUM(T$((Ord(T) >= 211) and (Ord(T) < 251) and (RADVOL(T) > -0.5)),
                            (NINO3(T)-AVENINO3)*(DCONCO2NAT(T)-AVEDCONCO2NAT))
                            /(SUM(T$((Ord(T) >= 211) and (Ord(T) < 251) and (RADVOL(T) > -0.5)),1)-1);

SXYQ2 $(ENSOCO2SWCH < 0.5)
                 .. SXY =E= 0.001;

SLOPECO2ENSOQ1 $(ENSOCO2SWCH > 0.5) .. SLOPECO2ENSO =E= SXY/SX2;

SLOPECO2ENSOQ2 $(ENSOCO2SWCH < 0.5) .. SLOPECO2ENSO =E= 0;

YINTRCO2ENSOQ1 $(ENSOCO2SWCH > 0.5) .. YINTRCO2ENSO =E= AVEDCONCO2NAT-SLOPECO2ENSO*AVENINO3;

YINTRCO2ENSOQ2 $(ENSOCO2SWCH < 0.5) .. YINTRCO2ENSO =E= 0;


RSQ11Q $(ENSOCO2SWCH > 0.5)
                 .. RSQ1 =E= SQR(SUM(T$((Ord(T) >= 211) and (Ord(T) < 251) and (RADVOL(T) > -0.5)),
                              (NINO3(T)-AVENINO3)*(DCONCO2NAT(T)-AVEDCONCO2NAT)))
                             /SUM(T$((Ord(T) >= 211) and (Ord(T) < 251) and (RADVOL(T) > -0.5)),SQR(NINO3(T)-AVENINO3))
                             /SUM(T$((Ord(T) >= 211) and (Ord(T) < 251) and (RADVOL(T) > -0.5)),
                              SQR(DCONCO2NAT(T)-AVEDCONCO2NAT));
* Excluding the year 2000 because DCONCO2NAT cannot be defined in 2000.

RSQ12Q $(ENSOCO2SWCH < 0.5) .. RSQ1 =E= 0;

RSQ21Q $(ENSOCO2SWCH > 0.5)
                 .. RSQ2 =E= 1-SUM(T$((Ord(T) >= 211) and (Ord(T) < 251) and (RADVOL(T) > -0.5)),
                             SQR(DCONCO2NAT(T)-DCONCO2ENSO(T)))
                             /SUM(T$((Ord(T) >= 211) and (Ord(T) < 251) and (RADVOL(T) > -0.5)),
                             SQR(DCONCO2NAT(T)-AVEDCONCO2NAT));
* Excluding the year 2000 because DCONCO2NAT cannot be defined in 2000.

RSQ22Q $(ENSOCO2SWCH < 0.5) .. RSQ2 =E= 0;


* Only for model output of CONCO2RESWTENSO(T).
CONCO2RESWTENSOQ(T)         .. CONCO2INT(T) =E= CONCO2PRI(T)+CONCO2RESWTENSO(T);


*** Radiative forcings -------------------------------------------------------------------------------------------------

RADHALQ(T,HALR)       .. RAD(T,HALR) =E= RHORAD(HALR)*CONSCN(T,HALR);
*RADHALQ(T,HALR)       .. RAD(T,HALR) =E= 0;
* CH3Cl and CH3Br are not included in the radiative forcing calculation.

RADSF6Q(T)            .. RAD(T,'SF6') =E= RHOSF6*CONSCN(T,'SF6');
*RADSF6Q(T)            .. RAD(T,'SF6') =E= 0;

RADO3Q(T)             .. RADO3(T) =E= RADO3SCN(T);
*RADO3Q(T)             .. RADO3(T) =E= 0;

RADSOLQ(T)            .. RADSOL(T) =E= RADSOLBLM(T);
*RADSOLQ(T)            .. RADSOL(T) =E= 0;

RADSOL2000Q           .. RADSOL2000 =E= RADSOL('2000');
*RADSOL2000Q           .. RADSOL2000 =E= 0;

RADSOL20AVEQ          .. RADSOL20AVE =E= SUM(T$((151 < Ord(T)) and (Ord(T) <= 251)),RADSOL(T))/100;
*RADSOL20AVEQ          .. RADSOL20AVE =E= 0;

RADAERQ(T)            .. RADAER(T) =E= RAD(T,'SO2')+RAD(T,'CRBAER')+RAD(T,'AERIND');

RADANTQ(T)            .. RADANT(T) =E=
                         RAD(T,'CO2')+RAD(T,'CH4')+RAD(T,'N2O')+RADHAL(T)+RAD(T,'SF6')+RAD(T,'H2OSTR')+RAD(T,'AERIND')+
                         RADO3(T)+
                         RAD(T,'SO2')+
                         RAD(T,'CRBAER');

RADANTLQ(T)           .. RADANTL(T) =E=
                         RAD(T,'CO2')+RAD(T,'CH4')+RAD(T,'N2O')+RADHAL(T)+RAD(T,'SF6')+RAD(T,'H2OSTR')+RAD(T,'AERIND')+
                         ROZO/(1-FLA+FLA*ROZO)*RADO3(T)+
                         RSUL/(1-FLA+FLA*RSUL)*RAD(T,'SO2')+
                         RFF/(1-FLA+FLA*RFF)*RAD(T,'CRBAER');

RADANTOQ(T)           .. RADANTO(T) =E=
                         RAD(T,'CO2')+RAD(T,'CH4')+RAD(T,'N2O')+RADHAL(T)+RAD(T,'SF6')+RAD(T,'H2OSTR')+RAD(T,'AERIND')+
                         1/(1-FLA+FLA*ROZO)*RADO3(T)+
                         1/(1-FLA+FLA*RSUL)*RAD(T,'SO2')+
                         1/(1-FLA+FLA*RFF)*RAD(T,'CRBAER');

RADTTLPRIQ(T)         .. RADTTLPRI(T) =E=
                         RAD(T,'CO2')+RAD(T,'CH4')+RAD(T,'N2O')+RADHAL(T)+RAD(T,'SF6')+RADO3(T)+RAD(T,'H2OSTR')
                         +[RAD(T,'SO2')+RAD(T,'CRBAER')+RAD(T,'AERIND')]/SCLRAD*SCLRADPRI+RADVOL(T)+RADSOL(T)
                         +RADMISPRI(T);
* It is calculated only for graphical presentation.
* RADTTLPRI is the total radiative forcing including the prior estimates of forcing scaling or missing forcing.
* RADTTLPRI is not exactly the prior information of the total radiative forcing
* as CO2, CH4, and N2O forcings are calculated based on the posterior parameter estimates.

RADTTLQ(T)            .. RADTTL(T) =E=
                         RAD(T,'CO2')+RAD(T,'CH4')+RAD(T,'N2O')+RADHAL(T)+RAD(T,'SF6')+RADO3(T)+RAD(T,'H2OSTR')
                         +RAD(T,'SO2')+RAD(T,'CRBAER')+RAD(T,'AERIND')+RADVOL(T)+RADSOL(T)+RADMIS(T);

RADAERAVEQ            .. RADAERAVE =E= SUM(T$(Ord(T) >= 201), RAD(T,'SO2')+RAD(T,'CRBAER')+RAD(T,'AERIND'))/51;

RADMISAVEQ            .. RADMISAVE =E= SUM(T$(Ord(T) >= 201), RADMIS(T))/51;


*** Temperature change -------------------------------------------------------------------------------------------------

T2M1750Q              .. T2M1750 =E= T2M('1750');

T2M1990Q              .. T2M1990 =E= T2M('1990');

TAL1750Q              .. TAL1750 =E= TAL('1750');

TSS1750Q              .. TSS1750 =E= TSS('1750');

T2M69Q                .. T2M69 =E= SUM(T$((Ord(T) GT (TENDINV-40)) and (Ord(T) LE (TENDINV-10))),T2M(T))/30;
* '30' means the number of time steps between 1990 and 1961 (if the inversion time horizon is till 2000).

T2M68Q                .. T2M68 =E= SUM(T$((Ord(T) GT (TENDINV-40)) and (Ord(T) LE (TENDINV-20))),T2M(T))/20;
* '20' means the number of time steps between 1980 and 1961 (if the inversion time horizon is till 2000).

T2MINS69Q             .. T2MINS69 =E= SUM(T$((Ord(T) GT (TENDINV-40)) and (Ord(T) LE (TENDINV-10))),T2MINS(T))/30;

T2MRCN68Q             .. T2MRCN68 =E= SUM(T$((Ord(T) GT (TENDINV-40)) and (Ord(T) LE (TENDINV-20))),T2MRCN(T))/20;

T2MPRIQ(T)            .. T2MPRI(T) =E= ((T2MINS(T)+T2M69-T2MINS69)$(Ord(T) GT 106))
                         +((T2MRCN(T)+(T2M68-T2MRCN68))$(Ord(T) LE 106));


*** Temperature change and ENSO ----------------------------------------------------------------------------------------

DT2MENSOQ1(T) $(Ord(T) < 181)  .. DT2MENSO(T) =E= 0;

DT2MENSOQ2(T) $(Ord(T) >= 181) .. DT2MENSO(T) =E= FLA*SOICOFL*SOI(T) + (1-FLA)*BSI*SOICOFO*SOI(T);


*** Equality for parameters and data -----------------------------------------------------------------------------------

CO2FFQ1(T) $(CO2FFSWCH > 0.5)             .. CO2FF(T) =E= CO2FFPRI(T)+CO2FFRES(T);

CO2FFQ2(T) $(CO2FFSWCH < 0.5)             .. CO2FF(T) =E= CO2FFPRI(T);

CO2FFQ3(T) $(CO2FFSWCH < 0.5)             .. CO2FFRES(T) =E= 0;

CO2LUQ1(T) $(CO2LUSWCH > 0.5)             .. CO2LU(T) =E= CO2LUPRI(T)+CO2LURES(T);

CO2LUQ2(T) $(CO2LUSWCH < 0.5)             .. CO2LU(T) =E= CO2LUPRI(T);

CO2LUQ3(T) $(CO2LUSWCH < 0.5)             .. CO2LURES(T) =E= 0;

ANTSPQ1(T,FOG) $(ANTCH4N2OSWCH > 0.5)     .. ANT(T,FOG) =E= ANTPRI(T,FOG)+ANTRES(T,FOG);

ANTSPQ2(T,FOG) $(ANTCH4N2OSWCH < 0.5)     .. ANT(T,FOG) =E= ANTPRI(T,FOG);

ANTSPQ3(T,FOG) $(ANTCH4N2OSWCH < 0.5)     .. ANTRES(T,FOG) =E= 0;

RADMISSPINQ1(T) $((RADMISSWCH > 0.5) and (RADMISAR1SWCH > 0.5) and (Ord(T) > 1))
                                 .. RADMIS(T) =E= RADMISPRI(T)+RADMISRES(T)+RADMISAR1PRO*(RADMIS(T-1)-RADMISPRI(T-1));

RADMISSPINQ2(T) $((RADMISSWCH > 0.5) and (RADMISAR1SWCH > 0.5) and (Ord(T) = 1))
                                          .. RADMIS(T) =E= RADMISPRI(T)+RADMISRES(T);

RADMISSPINQ3(T) $((RADMISSWCH > 0.5) and (RADMISAR1SWCH < 0.5))
                                          .. RADMIS(T) =E= RADMISPRI(T)+RADMISRES(T);

RADMISSPINQ4(T) $(RADMISSWCH < 0.5)       .. RADMIS(T) =E= RADMISPRI(T);

RADMISSPINQ5(T) $(RADMISSWCH < 0.5)       .. RADMISRES(T) =E= 0;


TMIX1750Q1 $(TMIX1750SWCH > 0.5)          .. TMIX1750 =E= TMIX1750PRI+TMIX1750RES;

TMIX1750Q2 $((TMIX1750SWCH < 0.5) or (OCNCO2TFSWCH < 0.5)) .. TMIX1750 =E= TMIX1750PRI;

TMIX1750Q3 $((TMIX1750SWCH < 0.5) or (OCNCO2TFSWCH < 0.5)) .. TMIX1750RES =E= 0;

DMPTMIXQ1 $(DMPTMIXSWCH > 0.5)            .. DMPTMIX =E= DMPTMIXPRI+DMPTMIXRES;

DMPTMIXQ2 $((DMPTMIXSWCH < 0.5) or (OCNCO2TFSWCH < 0.5)) .. DMPTMIX =E= DMPTMIXPRI;

DMPTMIXQ3 $((DMPTMIXSWCH < 0.5) or (OCNCO2TFSWCH < 0.5)) .. DMPTMIXRES =E= 0;

BETAQ1(T) $(BETASWCH > 0.5)               .. BETA(T) =E= BETAPRI+BETARES(T);

BETAQ2(T) $(BETASWCH < 0.5)               .. BETA(T) =E= BETAPRI;

BETAQ3(T) $(BETASWCH < 0.5)               .. BETARES(T) =E= 0;

Q10BASQ1 $(Q10BASSWCH > 0.5)              .. Q10BAS =E= Q10BASPRI+Q10BASRES;

Q10BASQ2 $((Q10BASSWCH < 0.5) or (LNDCO2TFRESSWCH < 0.5)) .. Q10BAS =E= Q10BASPRI;

Q10BASQ3 $((Q10BASSWCH < 0.5) or (LNDCO2TFRESSWCH < 0.5)) .. Q10BASRES =E= 0;

Q10TMPQ1 $(Q10TMPSWCH > 0.5)              .. Q10TMP =E= Q10TMPPRI+Q10TMPRES;

Q10TMPQ2 $((Q10TMPSWCH < 0.5) or (LNDCO2TFQ10SWCH < 0.5)) .. Q10TMP =E= Q10TMPPRI;

Q10TMPQ3 $((Q10TMPSWCH < 0.5) or (LNDCO2TFQ10SWCH < 0.5)) .. Q10TMPRES =E= 0;

TNPPSENSQ1 $(TNPPSENSSWCH > 0.5)          .. TNPPSENS =E= TNPPSENSPRI+TNPPSENSRES;

TNPPSENSQ2 $((TNPPSENSSWCH < 0.5) or (LNDCO2TFNPPSWCH < 0.5)) .. TNPPSENS =E= TNPPSENSPRI;

TNPPSENSQ3 $((TNPPSENSSWCH < 0.5) or (LNDCO2TFNPPSWCH < 0.5)) .. TNPPSENSRES =E= 0;

CO2OCNUPPREQ1(T) $(CO2OCNUPPRESWCH > 0.5) .. CO2OCNUPPRE(T) =E= CO2OCNUPPREPRI+CO2OCNUPPRERES(T);

CO2OCNUPPREQ2(T) $(CO2OCNUPPRESWCH < 0.5) .. CO2OCNUPPRE(T) =E= CO2OCNUPPREPRI;

CO2OCNUPPREQ3(T) $(CO2OCNUPPRESWCH < 0.5) .. CO2OCNUPPRERES(T) =E= 0;

CO2LNDUPPREQ1(T) $(CO2LNDUPPRESWCH > 0.5) .. CO2LNDUPPRE(T) =E= CO2LNDUPPREPRI+CO2LNDUPPRERES(T);

CO2LNDUPPREQ2(T) $(CO2LNDUPPRESWCH < 0.5) .. CO2LNDUPPRE(T) =E= CO2LNDUPPREPRI;

CO2LNDUPPREQ3(T) $(CO2LNDUPPRESWCH < 0.5) .. CO2LNDUPPRERES(T) =E= 0;

NATCH4Q1 $(NATCH4N2OSWCH > 0.5)           .. NATCH4 =E= NATCH4PRI+NATCH4RES;

NATCH4Q2 $(NATCH4N2OSWCH < 0.5)           .. NATCH4 =E= NATCH4PRI;

NATCH4Q3 $(NATCH4N2OSWCH < 0.5)           .. NATCH4RES =E= 0;

NATN2OQ1 $(NATCH4N2OSWCH > 0.5)           .. NATN2O =E= NATN2OPRI+NATN2ORES;

NATN2OQ2 $(NATCH4N2OSWCH < 0.5)           .. NATN2O =E= NATN2OPRI;

NATN2OQ3 $(NATCH4N2OSWCH < 0.5)           .. NATN2ORES =E= 0;

TAUCH4OHQ1 $(TAUCH4N2OSWCH > 0.5)         .. TAUCH4OH =E= TAUCH4OHPRI+TAUCH4OHRES;

TAUCH4OHQ2 $(TAUCH4N2OSWCH < 0.5)         .. TAUCH4OH =E= TAUCH4OHPRI;

TAUCH4OHQ3 $(TAUCH4N2OSWCH < 0.5)         .. TAUCH4OHRES =E= 0;

TAUN2OQ1 $(TAUCH4N2OSWCH > 0.5)           .. TAUN2O =E= TAUN2OPRI+TAUN2ORES;

TAUN2OQ2 $(TAUCH4N2OSWCH < 0.5)           .. TAUN2O =E= TAUN2OPRI;

TAUN2OQ3 $(TAUCH4N2OSWCH < 0.5)           .. TAUN2ORES =E= 0;

RAD2XCO2Q1 $(RAD2XCO2SWCH > 0.5)          .. RAD2XCO2 =E= RAD2XCO2PRI+RAD2XCO2RES;

RAD2XCO2Q2 $(RAD2XCO2SWCH < 0.5)          .. RAD2XCO2 =E= RAD2XCO2PRI;

RAD2XCO2Q3 $(RAD2XCO2SWCH < 0.5)          .. RAD2XCO2RES =E= 0;

SCLRADQ1 $(SCLRADSWCH > 0.5)              .. SCLRAD =E= SCLRADPRI+SCLRADRES;

SCLRADQ2 $(SCLRADSWCH < 0.5)              .. SCLRAD =E= SCLRADPRI;

SCLRADQ3 $(SCLRADSWCH < 0.5)              .. SCLRADRES =E= 0;

SENSQ1(T) $(SENSSWCH > 0.5)               .. SENS(T) =E= SENSPRI+SENSRES(T);

SENSQ2(T) $(SENSSWCH < 0.5)               .. SENS(T) =E= SENSPRI;

SENSQ3(T) $(SENSSWCH < 0.5)               .. SENSRES(T) =E= 0;

KAPPAQ1 $(KAPPASWCH > 0.5)                .. KAPPA =E= KAPPAPRI+KAPPARES;

KAPPAQ2 $(KAPPASWCH < 0.5)                .. KAPPA =E= KAPPAPRI;

KAPPAQ3 $(KAPPASWCH < 0.5)                .. KAPPARES =E= 0;

RADMISAR1PROQ1 $(RADMISAR1PROSWCH > 0.5)  .. RADMISAR1PRO =E= RADMISAR1PROPRI+RADMISAR1PRORES;

RADMISAR1PROQ2 $(RADMISAR1PROSWCH < 0.5)  .. RADMISAR1PRO =E= RADMISAR1PROPRI;

RADMISAR1PROQ3 $(RADMISAR1PROSWCH < 0.5)  .. RADMISAR1PRORES =E= 0;

T2MAR1PROQ1 $(T2MAR1PROSWCH > 0.5)        .. T2MAR1PRO =E= T2MAR1PROPRI+T2MAR1PRORES;

T2MAR1PROQ2 $(T2MAR1PROSWCH < 0.5)        .. T2MAR1PRO =E= T2MAR1PROPRI;

T2MAR1PROQ3 $(T2MAR1PROSWCH < 0.5)        .. T2MAR1PRORES =E= 0;


CO2OCNUPANTSPQ1(T) $(CO2UPANTSWCH > 0.5)  .. CO2OCNUPANT(T) =E= CO2OCNUPANTPRI(T)+CO2OCNUPANTRES(T);

CO2OCNUPANTSPQ3(T) $(CO2UPANTSWCH < 0.5)  .. CO2OCNUPANTRES(T) =E= 0;

CO2LNDUPANTSPQ1(T) $(CO2UPANTSWCH > 0.5)  .. CO2LNDUPANT(T) =E= CO2LNDUPANTPRI(T)+CO2LNDUPANTRES(T);

CO2LNDUPANTSPQ3(T) $(CO2UPANTSWCH < 0.5)  .. CO2LNDUPANTRES(T) =E= 0;

CONCO2SPINQ1(T) $(CONCO2SWCH > 0.5)       .. CON(T,'CO2') =E= CONCO2PRI(T)+CONCO2RES(T);

CONCO2SPINQ3(T) $(CONCO2SWCH < 0.5)       .. CONCO2RES(T) =E= 0;

CONCH4SPINQ1(T) $(CONCH4N2OSWCH > 0.5)    .. CON(T,'CH4') =E= CONCH4PRI(T)+CONCH4RES(T);

CONCH4SPINQ3(T) $(CONCH4N2OSWCH < 0.5)    .. CONCH4RES(T) =E= 0;

CONN2OSPINQ1(T) $(CONCH4N2OSWCH > 0.5)    .. CON(T,'N2O') =E= CONN2OPRI(T)+CONN2ORES(T);

CONN2OSPINQ3(T) $(CONCH4N2OSWCH < 0.5)    .. CONN2ORES(T) =E= 0;

T2MSPQ1(T) $((T2MSWCH > 0.5) and (T2MAR1SWCH > 0.5) and (Ord(T) > 1))
                                          .. T2M(T) =E= T2MPRI(T)+T2MRES(T)+T2MAR1PRO*(T2M(T-1)-T2MPRI(T-1));

T2MSPQ2(T) $((T2MSWCH > 0.5) and (T2MAR1SWCH > 0.5) and (Ord(T) = 1))
                                          .. T2M(T) =E= T2MPRI(T)+T2MRES(T);

T2MSPQ3(T) $((T2MSWCH > 0.5) and (T2MAR1SWCH < 0.5))
                                          .. T2M(T) =E= T2MPRI(T)+T2MRES(T);

T2MSPQ4(T) $(T2MSWCH < 0.5)               .. T2MRES(T) =E= 0;


*** Modifications of the prior uncertainty ranges ----------------------------------------------------------------------
*
* Original prior uncertainty ranges of parameters and data in literature are modified here for the ACC2 inversion.

CO2FFSGMMODQ(T)                           .. CO2FFSGMMOD(T) =E= [CO2FFSGMORG(T)];

CO2LUSGMMODQ1(T) $(CO2LUSGMSWCH > 0.5)    .. CO2LUSGMMOD(T) =E= [CO2LUSGMFAC*CO2LUSGMORG(T)];

CO2LUSGMMODQ2(T) $(CO2LUSGMSWCH < 0.5)    .. CO2LUSGMMOD(T) =E= [CO2LUSGMORG(T)];

ANTCH4SGMMOD1(T) $(ANTCH4N2OSGMSWCH > 0.5).. ANTSGMMOD(T,'CH4') =E= [ANTCH4N2OSGMFAC*ANTSGMORG(T,'CH4')];

ANTCH4SGMMOD2(T) $(ANTCH4N2OSGMSWCH < 0.5).. ANTSGMMOD(T,'CH4') =E= [ANTSGMORG(T,'CH4')];

ANTN2OSGMMOD1(T) $(ANTCH4N2OSGMSWCH > 0.5).. ANTSGMMOD(T,'N2O') =E= [ANTCH4N2OSGMFAC*ANTSGMORG(T,'N2O')];

ANTN2OSGMMOD2(T) $(ANTCH4N2OSGMSWCH < 0.5).. ANTSGMMOD(T,'N2O') =E= [ANTSGMORG(T,'N2O')];

RADMISSGMMODQ1(T) $(VOLSWCH > 0.5)
                                 .. RADMISSGMMOD(T) =E= [RADMISSGMORG(T)]$(RADVOL(T) > -0.5)+
                                                        [VOLSGMFAC*RADMISSGMORG(T)]$(RADVOL(T) <= -0.5);

RADMISSGMMODQ2(T) $(VOLSWCH < 0.5)
                                 .. RADMISSGMMOD(T) =E= RADMISSGMORG(T);


CONCO2SGMMODQ1(T) $((CONCO2SGMSWCH > 0.5) and (VOLSWCH > 0.5))
                                 .. CONCO2SGMMOD(T) =E=
                                    [CONCO2SGMORG(T)]$((RADVOL(T) >  -0.5) and (Ord(T) <  210))
                                    +[VOLSGMFAC*CONCO2SGMORG(T)]$((RADVOL(T) <= -0.5) and (Ord(T) <  210))
                                    +[CONCO2SGMFAC*CONCO2SGMORG(T)]$((RADVOL(T) >  -0.5) and (Ord(T) >= 210))
                                    +[CONCO2SGMFAC*VOLSGMFAC*CONCO2SGMORG(T)]$((RADVOL(T) <= -0.5) and (Ord(T) >= 210));

CONCO2SGMMODQ2(T) $((CONCO2SGMSWCH > 0.5) and (VOLSWCH < 0.5))
                                 .. CONCO2SGMMOD(T) =E= [CONCO2SGMORG(T)]$(Ord(T) <  210)
                                    +[CONCO2SGMFAC*CONCO2SGMORG(T)]$(Ord(T) >= 210);

CONCO2SGMMODQ3(T) $((CONCO2SGMSWCH < 0.5) and (VOLSWCH > 0.5))
                                 .. CONCO2SGMMOD(T) =E= [CONCO2SGMORG(T)]$(RADVOL(T) >  -0.5)
                                    +[VOLSGMFAC*CONCO2SGMORG(T)]$(RADVOL(T) <= -0.5);

CONCO2SGMMODQ4(T) $((CONCO2SGMSWCH < 0.5) and (VOLSWCH < 0.5))
                                 .. CONCO2SGMMOD(T) =E= CONCO2SGMORG(T);

CONCH4SGMMODQ1(T) $((CONCH4SGMSWCH > 0.5) and (VOLSWCH > 0.5))
                                 .. CONCH4SGMMOD(T) =E= [CONCH4SGMORG(T)]$((RADVOL(T) >  -0.5) and (Ord(T) <  235))
                                    +[VOLSGMFAC*CONCH4SGMORG(T)]$((RADVOL(T) <= -0.5) and (Ord(T) <  235))
                                    +[CONCH4SGMFAC*CONCH4SGMORG(T)]$((RADVOL(T) >  -0.5) and (Ord(T) >= 235))
                                    +[CONCH4SGMFAC*VOLSGMFAC*CONCH4SGMORG(T)]$((RADVOL(T) <= -0.5) and (Ord(T) >= 235));

CONCH4SGMMODQ2(T) $((CONCH4SGMSWCH > 0.5) and (VOLSWCH < 0.5))
                                 .. CONCH4SGMMOD(T) =E= [CONCH4SGMORG(T)]$(Ord(T) <  235)
                                    +[CONCH4SGMFAC*CONCH4SGMORG(T)]$(Ord(T) >= 235);

CONCH4SGMMODQ3(T) $((CONCH4SGMSWCH < 0.5) and (VOLSWCH > 0.5))
                                 .. CONCH4SGMMOD(T) =E= [CONCH4SGMORG(T)]$(RADVOL(T) >  -0.5)
                                    +[VOLSGMFAC*CONCH4SGMORG(T)]$(RADVOL(T) <= -0.5);

CONCH4SGMMODQ4(T) $((CONCH4SGMSWCH < 0.5) and (VOLSWCH < 0.5))
                                 .. CONCH4SGMMOD(T) =E= CONCH4SGMORG(T);

CONN2OSGMMODQ1(T) $((CONN2OSGMSWCH > 0.5) and (VOLSWCH > 0.5))
                                 .. CONN2OSGMMOD(T) =E= [CONN2OSGMORG(T)]$((RADVOL(T) >  -0.5) and (Ord(T) <  229))
                                    +[VOLSGMFAC*CONN2OSGMORG(T)]$((RADVOL(T) <= -0.5) and (Ord(T) <  229))
                                    +[CONN2OSGMFAC*CONN2OSGMORG(T)]$((RADVOL(T) >  -0.5) and (Ord(T) >= 229))
                                    +[CONN2OSGMFAC*VOLSGMFAC*CONN2OSGMORG(T)]$((RADVOL(T) <= -0.5) and (Ord(T) >= 229));

CONN2OSGMMODQ2(T) $((CONN2OSGMSWCH > 0.5) and (VOLSWCH < 0.5))
                                 .. CONN2OSGMMOD(T) =E= [CONN2OSGMORG(T)]$(Ord(T) <  229)
                                    +[CONN2OSGMFAC*CONN2OSGMORG(T)]$(Ord(T) >= 229);

CONN2OSGMMODQ3(T) $((CONN2OSGMSWCH < 0.5) and (VOLSWCH > 0.5))
                                 .. CONN2OSGMMOD(T) =E= [CONN2OSGMORG(T)]$(RADVOL(T) >  -0.5)
                                    +[VOLSGMFAC*CONN2OSGMORG(T)]$(RADVOL(T) <= -0.5);

CONN2OSGMMODQ4(T) $((CONN2OSGMSWCH < 0.5) and (VOLSWCH < 0.5))
                                 .. CONN2OSGMMOD(T) =E= CONN2OSGMORG(T);

T2MSGMMODQ1(T) $((VOLSWCH > 0.5) and (T2MAR1SWCH > 0.5))
                                 .. T2MSGMMOD(T) =E= [T2MSGMAR1(T)]$(RADVOL(T) > -0.5)
                                                     +[VOLSGMFAC*T2MSGMAR1(T)]$(RADVOL(T) <= -0.5);

T2MSGMMODQ2(T) $((VOLSWCH < 0.5) and (T2MAR1SWCH > 0.5))
                                 .. T2MSGMMOD(T) =E= T2MSGMAR1(T);

T2MSGMMODQ3(T) $((VOLSWCH > 0.5) and (T2MAR1SWCH < 0.5))
                                 .. T2MSGMMOD(T) =E= [T2MSGMORG(T)]$(RADVOL(T) > -0.5)
                                                     +[VOLSGMFAC*T2MSGMORG(T)]$(RADVOL(T) <= -0.5);

T2MSGMMODQ4(T) $((VOLSWCH < 0.5) and (T2MAR1SWCH < 0.5))
                                 .. T2MSGMMOD(T) =E= T2MSGMORG(T);


*** Residual terms in the cost function --------------------------------------------------------------------------------

CO2FFWSQRESQ1(T) $(CO2FFSWCH > 0.5)       .. CO2FFWSQRES(T) =E= 0.5*SQR(CO2FFRES(T)/CO2FFSGMMOD(T));

CO2FFWSQRESQ2(T) $(CO2FFSWCH < 0.5)       .. CO2FFWSQRES(T) =E= 0;

CO2LUWSQRESQ1(T) $(CO2LUSWCH > 0.5)       .. CO2LUWSQRES(T) =E= 0.5*SQR(CO2LURES(T)/CO2LUSGMMOD(T));

CO2LUWSQRESQ2(T) $(CO2LUSWCH < 0.5)       .. CO2LUWSQRES(T) =E= 0;

ANTWSQRESQ1(T,FOG) $(ANTCH4N2OSWCH > 0.5) .. ANTWSQRES(T,FOG) =E= 0.5*SQR(ANTRES(T,FOG)/ANTSGMMOD(T,FOG));

ANTWSQRESQ2(T,FOG) $(ANTCH4N2OSWCH < 0.5) .. ANTWSQRES(T,FOG) =E= 0;

RADMISWSQRESQ1(T) $(RADMISSWCH > 0.5)     .. RADMISWSQRES(T) =E= 0.5*SQR(RADMISRES(T)/RADMISSGMMOD(T));

RADMISWSQRESQ2(T) $(RADMISSWCH < 0.5)     .. RADMISWSQRES(T) =E= 0;


CO2OCNUPANTWSQRESQ1(T) $(CO2UPANTSWCH > 0.5) .. CO2OCNUPANTWSQRES(T) =E= 0.5*SQR(CO2OCNUPANTRES(T)/CO2OCNUPANTSGM(T));

CO2OCNUPANTWSQRESQ2(T) $(CO2UPANTSWCH < 0.5) .. CO2OCNUPANTWSQRES(T) =E= 0;

CO2LNDUPANTWSQRESQ1(T) $(CO2UPANTSWCH > 0.5) .. CO2LNDUPANTWSQRES(T) =E= 0.5*SQR(CO2LNDUPANTRES(T)/CO2LNDUPANTSGM(T));

CO2LNDUPANTWSQRESQ2(T) $(CO2UPANTSWCH < 0.5) .. CO2LNDUPANTWSQRES(T) =E= 0;

CONCO2WSQRESQ1(T) $(CONCO2SWCH > 0.5)     .. CONCO2WSQRES(T) =E= 0.5*SQR(CONCO2RES(T)/CONCO2SGMMOD(T));

CONCO2WSQRESQ2(T) $(CONCO2SWCH < 0.5)     .. CONCO2WSQRES(T) =E= 0;

CONCH4WSQRESQ1(T) $(CONCH4N2OSWCH > 0.5)  .. CONCH4WSQRES(T) =E= 0.5*SQR(CONCH4RES(T)/CONCH4SGMMOD(T));

CONCH4WSQRESQ2(T) $(CONCH4N2OSWCH < 0.5)  .. CONCH4WSQRES(T) =E= 0;

CONN2OWSQRESQ1(T) $(CONCH4N2OSWCH > 0.5)  .. CONN2OWSQRES(T) =E= 0.5*SQR(CONN2ORES(T)/CONN2OSGMMOD(T));

CONN2OWSQRESQ2(T) $(CONCH4N2OSWCH < 0.5)  .. CONN2OWSQRES(T) =E= 0;

T2MWSQRESQ1(T) $(T2MSWCH > 0.5)           .. T2MWSQRES(T) =E= 0.5*SQR(T2MRES(T)/T2MSGMMOD(T));

T2MWSQRESQ2(T) $(T2MSWCH < 0.5)           .. T2MWSQRES(T) =E= 0;


CO2FFRESSUMQ1 $(CO2FFSWCH > 0.5)          .. CO2FFRESSUM =E= SUM(T$(Ord(T) <= TENDINV),CO2FFWSQRES(T));

CO2FFRESSUMQ2 $(CO2FFSWCH < 0.5)          .. CO2FFRESSUM =E= 0;

CO2LURESSUMQ1 $(CO2LUSWCH > 0.5)          .. CO2LURESSUM =E= SUM(T$(Ord(T) <= TENDINV),CO2LUWSQRES(T));

CO2LURESSUMQ2 $(CO2LUSWCH < 0.5)          .. CO2LURESSUM =E= 0;

ANTRESSUMCH4Q1 $(ANTCH4N2OSWCH > 0.5)     .. ANTRESSUMCH4 =E= SUM(T$(Ord(T) <= TENDINV),ANTWSQRES(T,'CH4'));

ANTRESSUMCH4Q2 $(ANTCH4N2OSWCH < 0.5)     .. ANTRESSUMCH4 =E= 0;

ANTRESSUMN2OQ1 $(ANTCH4N2OSWCH > 0.5)     .. ANTRESSUMN2O =E= SUM(T$(Ord(T) <= TENDINV),ANTWSQRES(T,'N2O'));

ANTRESSUMN2OQ2 $(ANTCH4N2OSWCH < 0.5)     .. ANTRESSUMN2O =E= 0;

RADMISRESSUMQ1 $(RADMISSWCH > 0.5)        .. RADMISRESSUM =E= SUM(T$(Ord(T) <= TENDINV),RADMISWSQRES(T));

RADMISRESSUMQ2 $(RADMISSWCH < 0.5)        .. RADMISRESSUM =E= 0;


TMIX1750RESSUMQ1 $(TMIX1750SWCH > 0.5)    .. TMIX1750RESSUM =E= 0.5*SQR(TMIX1750RES/TMIX1750SGM);

TMIX1750RESSUMQ2 $(TMIX1750SWCH < 0.5)    .. TMIX1750RESSUM =E= 0;

DMPTMIXRESSUMQ1 $(DMPTMIXSWCH > 0.5)      .. DMPTMIXRESSUM =E= 0.5*SQR(DMPTMIXRES/DMPTMIXSGM);

DMPTMIXRESSUMQ2 $(DMPTMIXSWCH < 0.5)      .. DMPTMIXRESSUM =E= 0;

BETARESSUMQ1 $(BETASWCH > 0.5)            .. BETARESSUM =E= 0.5*SUM(T,SQR(BETARES(T)/BETASGM))/Card(T);
* Division of Card(T) (number of time steps in the past run) is needed
* because CO2 fertilization factor is in principle and by default time invariant
* although it can be experimentally varied in the model.

BETARESSUMQ2 $(BETASWCH < 0.5)            .. BETARESSUM =E= 0;

Q10BASRESSUMQ1 $(Q10BASSWCH > 0.5)        .. Q10BASRESSUM =E= 0.5*SQR(Q10BASRES/Q10BASSGM);

Q10BASRESSUMQ2 $(Q10BASSWCH < 0.5)        .. Q10BASRESSUM =E= 0;

Q10TMPRESSUMQ1 $(Q10TMPSWCH > 0.5)        .. Q10TMPRESSUM =E= 0.5*SQR(Q10TMPRES/Q10TMPSGM);

Q10TMPRESSUMQ2 $(Q10TMPSWCH < 0.5)        .. Q10TMPRESSUM =E= 0;

TNPPSENSRESSUMQ1 $(TNPPSENSSWCH > 0.5)    .. TNPPSENSRESSUM =E= 0.5*SQR(TNPPSENSRES/TNPPSENSSGM);

TNPPSENSRESSUMQ2 $(TNPPSENSSWCH < 0.5)    .. TNPPSENSRESSUM =E= 0;

CO2OCNUPPRERESSUMQ1 $(CO2OCNUPPRESWCH > 0.5)
                                 .. CO2OCNUPPRERESSUM =E= 0.5*SUM(T,SQR(CO2OCNUPPRERES(T)/CO2OCNUPPRESGM))/Card(T);
* Card(T) (number of time steps in the past run) is needed
* because preindustrial ocean CO2 uptake is assumed to be constant over time.

CO2OCNUPPRERESSUMQ2 $(CO2OCNUPPRESWCH < 0.5) .. CO2OCNUPPRERESSUM =E= 0;

CO2LNDUPPRERESSUMQ1 $(CO2LNDUPPRESWCH > 0.5)
                                 .. CO2LNDUPPRERESSUM =E= 0.5*SUM(T,SQR(CO2LNDUPPRERES(T)/CO2LNDUPPRESGM))/Card(T);
* Card(T) (number of time steps in the past run) is needed
* when preindustrial land CO2 uptake is assumed to be constant over time.

CO2LNDUPPRERESSUMQ2 $(CO2LNDUPPRESWCH < 0.5) .. CO2LNDUPPRERESSUM =E= 0;

NATCH4RESSUMQ1 $(NATCH4N2OSWCH > 0.5)     .. NATCH4RESSUM =E= 0.5*SQR(NATCH4RES/NATCH4SGM);

NATCH4RESSUMQ2 $(NATCH4N2OSWCH < 0.5)     .. NATCH4RESSUM =E= 0;

NATN2ORESSUMQ1 $(NATCH4N2OSWCH > 0.5)     .. NATN2ORESSUM =E= 0.5*SQR(NATN2ORES/NATN2OSGM);

NATN2ORESSUMQ2 $(NATCH4N2OSWCH < 0.5)     .. NATN2ORESSUM =E= 0;

TAUCH4OHRESSUMQ1 $(TAUCH4N2OSWCH > 0.5)   .. TAUCH4OHRESSUM =E= 0.5*SQR(TAUCH4OHRES/TAUCH4OHSGM);

TAUCH4OHRESSUMQ2 $(TAUCH4N2OSWCH < 0.5)   .. TAUCH4OHRESSUM =E= 0;

TAUN2ORESSUMQ1 $(TAUCH4N2OSWCH > 0.5)     .. TAUN2ORESSUM =E= 0.5*SQR(TAUN2ORES/TAUN2OSGM);

TAUN2ORESSUMQ2 $(TAUCH4N2OSWCH < 0.5)     .. TAUN2ORESSUM =E= 0;

RAD2XCO2RESSUMQ1 $(RAD2XCO2SWCH > 0.5)    .. RAD2XCO2RESSUM =E= 0.5*SQR(RAD2XCO2RES/RAD2XCO2SGM);

RAD2XCO2RESSUMQ2 $(RAD2XCO2SWCH < 0.5)    .. RAD2XCO2RESSUM =E= 0;

SCLRADRESSUMQ1 $(SCLRADSWCH > 0.5)        .. SCLRADRESSUM =E= 0.5*SQR(SCLRADRES/SCLRADSGM);

SCLRADRESSUMQ2 $(SCLRADSWCH < 0.5)        .. SCLRADRESSUM =E= 0;

SENSRESSUMQ1 $(SENSSWCH > 0.5)            .. SENSRESSUM =E= 0.5*SUM(T,SQR(SENSRES(T)/SENSSGM))/Card(T);
* Division of Card(T) (number of time steps in the past run) is needed
* because climate sensitivity is in principal time invariant although it can be experimentally varied in the model.

SENSRESSUMQ2 $(SENSSWCH < 0.5)            .. SENSRESSUM =E= 0;

KAPPARESSUMQ1 $(KAPPASWCH > 0.5)          .. KAPPARESSUM =E= 0.5*SQR(KAPPARES/KAPPASGM);

KAPPARESSUMQ2 $(KAPPASWCH < 0.5)          .. KAPPARESSUM =E= 0;

RADMISAR1PRORESSUMQ1 $(RADMISAR1PROSWCH > 0.5)  .. RADMISAR1PRORESSUM =E= 0.5*SQR(RADMISAR1PRORES/RADMISAR1PROSGM);

RADMISAR1PRORESSUMQ2 $(RADMISAR1PROSWCH < 0.5)  .. RADMISAR1PRORESSUM =E= 0;

T2MAR1PRORESSUMQ1 $(T2MAR1PROSWCH > 0.5)  .. T2MAR1PRORESSUM =E= 0.5*SQR(T2MAR1PRORES/T2MAR1PROSGM);

T2MAR1PRORESSUMQ2 $(T2MAR1PROSWCH < 0.5)  .. T2MAR1PRORESSUM =E= 0;


CO2OCNUPANTRESSUMQ1 $(CO2UPANTSWCH > 0.5) .. CO2OCNUPANTRESSUM =E= SUM(T$(Ord(T) <= TENDINV),CO2OCNUPANTWSQRES(T));

CO2OCNUPANTRESSUMQ2 $(CO2UPANTSWCH < 0.5) .. CO2OCNUPANTRESSUM =E= 0;

CO2LNDUPANTRESSUMQ1 $(CO2UPANTSWCH > 0.5) .. CO2LNDUPANTRESSUM =E= SUM(T$(Ord(T) <= TENDINV),CO2LNDUPANTWSQRES(T));

CO2LNDUPANTRESSUMQ2 $(CO2UPANTSWCH < 0.5) .. CO2LNDUPANTRESSUM =E= 0;

CONCO2RESSUMQ1 $(CONCO2SWCH > 0.5)        .. CONCO2RESSUM =E= SUM(T$(Ord(T) <= TENDINV),CONCO2WSQRES(T));

CONCO2RESSUMQ2 $(CONCO2SWCH < 0.5)        .. CONCO2RESSUM =E= 0;

CONCH4RESSUMQ1 $(CONCH4N2OSWCH > 0.5)     .. CONCH4RESSUM =E= SUM(T$(Ord(T) <= TENDINV),CONCH4WSQRES(T));

CONCH4RESSUMQ2 $(CONCH4N2OSWCH < 0.5)     .. CONCH4RESSUM =E= 0;

CONN2ORESSUMQ1 $(CONCH4N2OSWCH > 0.5)     .. CONN2ORESSUM =E= SUM(T$(Ord(T) <= TENDINV),CONN2OWSQRES(T));

CONN2ORESSUMQ2 $(CONCH4N2OSWCH < 0.5)     .. CONN2ORESSUM =E= 0;

T2MRESSUMQ1 $(T2MSWCH > 0.5)              .. T2MRESSUM =E= SUM(T$(Ord(T) <= TENDINV),T2MWSQRES(T));

T2MRESSUMQ2 $(T2MSWCH < 0.5)              .. T2MRESSUM =E= 0;


****** Variable scalings ===============================================================================================

ANT.SCALE(T,'CH4') = 100;

ANT.SCALE(T,'CO') = 100;

ANT.SCALE(T,'VOC') = 100;

CON.SCALE(T,'CO2') = 100;

CON.SCALE(T,'CH4') = 100;

CON.SCALE(T,'N2O') = 100;

DT2MENSO.SCALE(T) = 0.01;

CO2FFSGMMOD.SCALE(T) = 0.01;

CO2LUSGMMOD.SCALE(T) = 0.01;

CONCO2SGMMOD.SCALE(T) = 0.1;

CO2FFRES.SCALE(T) = 0.01;

CO2LURES.SCALE(T) = 0.01;

CO2OCNUPANTRES.SCALE(T) = 0.01;

CO2LNDUPANTRES.SCALE(T) = 0.01;

NATCH4RES.SCALE = 100;

CONCO2RES.SCALE(T) = 0.01;

CONCH4RES.SCALE(T) = 0.01;

CONN2ORES.SCALE(T) = 0.01;

RADO3.SCALE(T) = 0.01;

RADMISRES.SCALE(T) = 0.01;

T2MRES.SCALE(T) = 0.01;

CO2LURESSUM.SCALE = 100;

RADMISRESSUM.SCALE = 100;

T2MRESSUM.SCALE = 100;
