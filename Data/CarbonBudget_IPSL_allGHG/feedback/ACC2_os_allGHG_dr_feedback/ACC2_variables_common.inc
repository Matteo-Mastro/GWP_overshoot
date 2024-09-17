********** ACC2 version 4.2: ACC2_variables_common.inc *****************************************************************
********** Variable declarations (past and future modes) ===============================================================

VARIABLES

*** Ocean CO2 uptake ---------------------------------------------------------------------------------------------------

* The following dissociation constants are defined with molar units rather than molal units,
* which are conventionally used. See the relevant equations for unit conversion terms.

AK1S(T)         Dissociation constant of carbonate acid in actual seawater [mol per Giga L]
* To convert in the conventional unit [mol per kg], use AK1S(T)/(10**9*DENS)

AK1P(T)         Dissociation constant of carbonate acid in pure water   [mol per Giga L]

AK1A            Coefficient for AK1P to adjust for actual seawater

AK1B            Coefficient for AK1P to adjust for actual seawater

AK1C            Coefficient for AK1P to adjust for actual seawater


AK2S(T)         Dissociation constant of bicarbonate in actual seawater [mol per Giga L]
* To convert in the conventional unit [mol per kg], use AK2S(T)/(10**9*DENS)

AK2P(T)         Dissociation constant of bicarbonate in pure water      [mol per Giga L]

AK2A            Coefficient for AK2P to adjust for actual seawater

AK2B            Coefficient for AK2P to adjust for actual seawater

AK2C            Coefficient for AK2P to adjust for actual seawater


AKB(T)          Dissociation constant of boric acid                     [mol per Giga L]
* To convert in the conventional unit [mol per kg], use AKB(T)/(10**9*DENS)

AKW(T)          Dissociation constant of water                          [mol^2 per  (Giga L)^2]
* To convert in the conventional unit [(mol per kg)^2], use AKW(T)/(10**9*DENS)**2

ALPHAS(T)       Inverse of Henrys Law constant                          [ppm per (mol per Giga L)]
* Note that the definition of the Henrys Law constant: [CO2aq]=HCO2*ppCO2
* (e.g. equation (5) in Millero, 1995; equation (7.10) in Pilson, 1998)]
* To convert in the conventional unit [(mol per kg) per ppm], use 1/ALPHAS(T)/(10**9*DENS)*10^6

TMIX1750        Preindustrial mixed-layer temperature                   [degree in C]

H(T)            H+ concentration                                        [mol per Giga L]

CONCO2INT(T)    Atmospheric CO2 concentration (without ENSO effect)     [ppm]
* Internally calculated CONCO2INT(T) rather than final CON(T,'CO2') is used for box model calculations.
* (see the explanation under Concentration in ACC2_equations_common.inc)

PH(T)           PH of the ocean mixed layer                             [1]

DIC(T)          Concentration of DIC                                    [mol per kg]

CO2AQ(T)        Concentration of CO2aq                                  [mol per kg]

HCO3(T)         Concentration of HCO3-                                  [mol per kg]

CO32(T)         Concentration of CO32-                                  [mol per kg]

FRCO2AQ(T)      Fractional concentration of CO2aq                       [mol per kg]

FRHCO3(T)       Fractional concentration of HCO3-                       [mol per kg]

FRCO32(T)       Fractional concentration of CO32-                       [mol per kg]

COC(T,NOC)      Amount of anthropogenic carbon (NOC-layer of the analogue model)  [GtC]

COC1ALL(T)      Amount of preindustrial and anthropogenic carbon in the composite layer [GtC]

COC1ATM(T)      Amount of anthropogenic carbon in the atmosphere                  [GtC]

COC1ATMBAC      Amount of preindustrial carbon in the atmosphere                  [GtC]

COC1MIX(T)      Amount of anthropogenic carbon in the ocean mixed layer           [GtC]

COC1MIXBAC      Amount of preindustrial carbon in the ocean mixed layer           [GtC]
* This variable is estimated to be approximately 504.025863, according to SIAM input used for ICM1.1,
* In ACC2, the estimate of this variable is determined by the thermodynamic equilibrium
* in the preindustrial state of the carbon cycle system and with the assumed preindustrial CO2 concentration (= 277ppm).

REVELLE(T)      Revelle factor                                          [1]

DMPTMIX         Atmosphere-mixed layer temperature scaling factor       [1]

TMIX(T)         Average temperature in the mixed later                  [degree in C]

CO2OCNUPPRE(T)  Magnitude of preindustrial ocean CO2 uptake             [GtC per year]

CO2OCNUPANT(T)  Magnitude of ocean CO2 uptake                           [GtC per year]

CO2OCNUPTTL(T)  Magnitude of total ocean CO2 uptake                     [GtC per year]

SUMCOC(T)       Amount of carbon in the ocean relative to 1750          [GtC]

SUMCOC1999      Amount of carbon in the ocean in 1999 relative to 1750  [GtC]
* This variable is used to calculate the land CO2 uptake in 2000 in the future mode.


*** Land CO2 uptake ----------------------------------------------------------------------------------------------------

CTR(T,NTR)      Amount of preindustrial and anthropogenic carbon in the land biosphere [GtC]

CTR1750(NTR)    Amount of land carbon in 1750                           [GtC]

NPPBAC(T,NTR)   Magnitude of preindustrial NPP to the land biosphere    [GtC]

BETA(T)         Beta factor for CO2 fertilization                       [1]
* Time-invariant by default as controlled by equation BETACONT(T).

CO2FRT(T,NTR)   Magnitude of CO2 fertilization                          [GtC per year]

HTRRES(T,NTR)   Magnitude of heterotrophic respiration                  [GtC per year]

FLXCOC1(T)      Magnitude of total carbon flux for COC1                 [GtC per year]

FLXCOC2(T)      Magnitude of total carbon flux for COC2                 [GtC per year]

FLXCOC3(T)      Magnitude of total carbon flux for COC3                 [GtC per year]

FLXCOC4(T)      Magnitude of total carbon flux for COC4                 [GtC per year]

FLXCTR(T,NTR)   Magnitude of total carbon flux for CTR(NTR)             [GtC per year]

TNPPSENS        Sensitivity of terrestrial NPP to temperature change    [GtC per year per degree in C]

TNPPFAC(T)      Auxiliary variable related to temperature feedback to NPP [degree in C]

Q10EFF(T)       Effective Q10 factor for temperature feedback to land biosphere metabolism [1]

Q10BAS          Q10 for heterotrophic respiration                       [1]

Q10TMP          Temperature sensitivity to Q10                          [1]

TRESFAC(T)      Auxiliary variable related to temperature feedback to respiration [1]

CO2LNDUPPRE(T)  Magnitude of preindustrial land CO2 uptake              [GtC per year]

CO2LNDUPANT(T)  Magnitude of land CO2 uptake                            [GtC per year]

CO2LNDUPTTL(T)  Magnitude of total land CO2 uptake                      [GtC per year]

SUMCTR(T)       Amount of carbon in the land biosphere relative to 1750 [GtC]

SUMCTR1999      Amount of carbon in the land biosphere in 1999 relative to 1750 [GtC]
* This variable is used to calculate the land CO2 uptake in 2000 in the future mode.

CTRREL(T,NTR)   Changes in the land carbon amount relative to the amount in 1750 [GtC]


*** Carbon budget ------------------------------------------------------------------------------------------------------

CUMCO2(T)       Cummulative anthropogenic CO2 emissions                 [GtC]

CUMCO21999      Cummulative anthropogenic CO2 emissions by year 1999    [GtC]
* This variable is used to calculate the cumulative anthropogenic CO2 emissions in 2000 in the future mode.


*** Emissions (other than CO2) -----------------------------------------------------------------------------------------

ANTPRI(T,FOG)   Emissions of CH4 and N2O between 1750 and 2000                   [Units depending on the gases]
* For units, see the declaration of variable ANT(T,FOB) in ACC2_main_past.gms or ACC2_main_future.gms

NATCH4          Natural methane emission                                         [TgCH4 per year]

NATN2O          Natural nitrous oxide emission                                   [TgN per year]


*** Concentrations -----------------------------------------------------------------------------------------------------

DCONCO2ENSO(T)  Change in the atmospheric CO2 concentration due to ENSO          [GtC per year]

DCONCO2ENSO2000 Change in the atmospheric CO2 concentration due to ENSO in 2000  [GtC per year]

* The following concentrations prior 2000 are used only for the calculation of EESCl (in the future mode).
* They are calculated in the past mode and transferred to the future mode in ACC2_init_future.inc.
CON1999(ODS)    Atmospheric concentrations of ozone depleting substances in 1999

CON1998(ODS)    Atmospheric concentrations of ozone depleting substances in 1998

CON1997(ODS)    Atmospheric concentrations of ozone depleting substances in 1997


TAUCH4OH        Methane lifetime with respect to OH depletion            [year]

TAUN2O          Nitrous oxide lifetime                                   [year]


*** Radiative forcings -------------------------------------------------------------------------------------------------

RAD2XCO2        Radiative forcing of a doubling of atmospheric CO2 concentration         [W per m^2]

RAD(T,FOA)      Various components of total radiative forcing           [W per m^2]

RADAER(T)       Total aerosol forcing                                   [W per m^2]

RADTTL(T)       Combined radiative forcing of all components            [W per m^2]

RADANT(T)       Combined radiative forcing of all anthropogenic components               [W per m^2]
* Missing forcing is not included in RADANT(T).

RADHAL(T)       Combined radiative forcing of halocarbons               [W per m^2]

RADODS(T)       Combined radiative forcing of ozone depleting substances [W per m^2]

RADSOL(T)       Solar forcing                                           [W per m^2]

RADMIS(T)       Missing forcing                                         [W per m^2]

RADSOL2000      Solar forcing in 2000 used in the future mode           [W per m^2]

RADSOL20AVE     Solar forcing averaged between 1901 and 2000            [W per m^2]

RADTTLPRI(T)    Total radiative forcing using prior missing radiative forcing            [W per m^2]
* Note that this variable is not prior information of the total radiative forcing.
* This is, in other words, the posterior total radiative forcing except for the posterior missing radiative forcing.

RADAERAVE       Average of the aerosol forcing between 1950 and 2000    [W per m^2]

RADMISAVE       Average of the missing forcing between 1950 and 2000    [W per m^2]

SCLRAD          Scaling factor for radiative forcing                    [1]
* This variable is used to scale all the radiative forcing except for volcanic, solar, and missing forcings.
* See SCLRADSWCH in the switchboard file for further explanation.


*** Normalized radiative forcing and its time derivative ---------------------------------------------------------------

NORRAD(T)       Normalized radiative forcing                            [1]

SENS(T)         Climate sensitivity                                     [degree in C]
* Time-invariant by default as controlled by equation SENSCONT(T).


*** Principal components of Empirical Orthogonal Function --------------------------------------------------------------

SCALEMEANT      Factor translating first principal component of temperature in a global-mean value [1]

SENSSCALE(T)    Scale factor for the first past period                  [1]

SQROOT(T)       Factor for calculation of RPCTA and RPCTT for the first past period [1]

SCALEMEANS      Factor translating first principal component of sea level rise in a global-mean value [1]


*** Temperature change -------------------------------------------------------------------------------------------------

RPCTA(T,NPCR)    Coefficients in R_PCT of impulse response function for the first past period  [degree in C]

RPCTT(T,NPCR)    Time scales in R_PCT of impulse response function for the first past period   [year]

PCT(T)           Sum of principal components of temperature change      [degree in C]

PCTRES(T,NPCR)   Principal components of temperature change             [degree in C]

T2M(T)           Global-mean surface temperature change                 [degree in C]

T2M1750          Temperature in 1750 obtained during past phase         [degree in C]
* (needed for temperature feedback to land CO2 uptake)

T2M1990          Temperature in 1990 obtained during past phase         [degree in C]
* (needed for sea level calculation)

T2MRT1990(T)     Temperature change relative to 1990                    [degree in C]

* Note from Thomas Bruckner in ICM1.1:
* A comparison of the 3 box model and 2 box model approaches discussed in
* Hooss et al. (1999) pp. 38-40 and pp. 43-47, respectively, revealed
* that the 2 box model results are closer to IPCC Second Assessment Report
* results than those determined by the 3 box model.
* Consequently, ICM adopts the 2 box variant of NICCS.

*** Autoregressive model -----------------------------------------------------------------------------------------------

RADMISAR1PRO    Autoregressive propagator for missing forcing residuals                         [1]

T2MAR1PRO       Autoregressive propagator for global-mean surface air temperature residuals     [1]


*** Sea level rise (due to ocean thermal expansion)---------------------------------------------------------------------

PCS(T)          Sea level rise due to ocean thermal expansion                                   [m]

PCSRES(T,NPCR)  Principal components of sea level rise                                          [m]

PCS1990         Sea level rise due to ocean thermal expansion in 1990                           [m]

PCSRT1990(T)    Sea level rise due to ocean thermal expansion relative to the level in 1990     [m]


*** Sea level rise (due to factors other than ocean thermal expansion) -------------------------------------------------

SLRGLWT(T)      Sea level rise due to loss of mass of glaciers and small ice caps without glacier area contraction [m]

SLRGLWC(T)      Sea level rise due to loss of mass of glaciers and small ice caps with glacier area contraction [m]

SLRGLWCRT1990(T) Sea level rise due to loss of mass of glaciers and small ice caps with glacier area contraction
*               relative to 1990 [m]

SLRGLWT2000     Sea level rise in 2000 due to loss of mass of glaciers and small ice caps without glacier area
*               contraction [m]

SLRGLWC2000     Sea level rise in 2000 due to loss of mass of glaciers and small ice caps with glacier area contraction
*               [m]

SLRGLWC1990     Sea level rise in 1990 due to loss of mass of glaciers and small ice caps with glacier area contraction
*               [m]

SLRGR(T)        Greenland ice sheet contribution to sea level rise                              [m]

SLRGRRT1990(T)  Greenland ice sheet contribution to sea level rise relative to 1990             [m]

SLRAN(T)        Antarctic ice sheet contribution to sea level rise                              [m]

SLRANRT1990(T)  Antarctic ice sheet contribution to sea level rise relative to 1990             [m]

SLRABS(T)       Sea level rise due to four major factors since 1750                             [m]

SLRRT1750(T)    Sea level rise due to four major factors relative to 1750                       [m]

SLR1990         Sea level rise due to four major factors in 1990 relative to 1750               [m]

SLRRT1990(T)    Sea level rise due to four major and three minor factors relative to 1990       [m]


*** DOECLIM ------------------------------------------------------------------------------------------------------------

* Heat transfer coefficients connecting climate feedback strength over land (LP1, LP2)
* and ocean (OP1, OP2) and land-sea heat exchange (KP1, KP2) with global
* climate sensitivity (see Equations A.19, A.20, A.21 in Kriegler, 2005)

KP1             Heat transfer coefficients                       [W per m^2 K]

KP2             Heat transfer coefficients                       [W per m^2]

LP1             Heat transfer coefficients                       [W per m^2 K]

LP2             Heat transfer coefficients                       [W per m^2]

OP1             Heat transfer coefficients                       [W per m^2 K]

OP2             Heat transfer coefficients                       [W per m^2]

* Model parameters dependent on tco. See equations A.19-A.21 in Kriegler, 2005.

KSL             Land sea heat exchange coefficient               [W per m^2 K]

OLAM            Climate feedback strength over ocean             [W per m^2 K]

LLAM            Climate feedback strength over land              [W per m^2 K]

* Matrix entries for implicit numerical integration of DOECLIM   [W per m^2 K]
* Dependent on tco, kappa. See Equation A.27 in Kriegler, 2005.

A11

A12

A21

A22

B11

B12

B21

B22

* The following correction terms for 4th order implicit Runge-Kutta added by Elmar Kriegler (February 02, 2006)

C11

C12

C21

C22

* Radiative forcings over land and ocean

RADANTL(T)       Radiative forcing due to all the anthropogenic factors over the land    [W per m^2]

RADANTO(T)       Radiative forcing due to all the anthropogenic factors over the ocean   [W per m^2]

RADSOLL(T)       Solar forcing over the land                                             [W per m^2]

RADSOLO(T)       Solar forcing over the ocean                                            [W per m^2]

RADTTLL(T)       Total radiative forcing over the land                                   [W per m^2]

RADTTLO(T)       Total radiative forcing over the ocean                                  [W per m^2]

* The followings are associated with numerical implementation of forcing (added by Elmar Kriegler, 02.02.2006)

RADNUML(t)

RADNUMO(t)

* Temperature-related variables

TSS(T)           Sea surface temperature anomaly                                         [K]

TAL(T)           Land surface air temperature                                            [K]

GMT(T)           Global mean temperature anomaly                                         [K]

TSUM(T)          Temperature memory                                                      [K]

TSSMM(I)         TSS data used for temperature memory calculation                        [K]
* Note that the index is used in an opposite way of TSS(T).


* Convolution vectors

KER(I)

KER1(I)

*KER2(I)
*
*KER3(I)

KDIF             Ocean vertical diffusivity                                              [m2 per a]

KAPPA            Ocean vertical diffusivity                                              [cm2 per s]

TAL1750          Land surface air temperature in 1750                                    [K]

TSS1750          Sea surface air temperature in 1750                                     [K]

DT2MENSO(T)      El Nino related effect in global mean surface temperature               [K]


* Heat flux (for output)

FO(T)           Heat flux to the interior ocean from the mixed layer                     [W per m2]

HCIO(T)          Global total heat flux to the interior ocean from the mixed layer       [ZJ]
* ZJ is zetta joules i.e. 10^21J.

HCML(T)          Change in the heat content in the mixed-layer ocean-atmosphere system   [ZJ]

HCAL(T)          Change in the heat content in the land-atmosphere system                [ZJ]

HCTL(T)          Change in the heat content in the earth system                          [ZJ]

HCRD(T)          Heat input from the radiative forcing to the system                     [ZJ];


****** Variable scalings ===============================================================================================

AK1S.SCALE(T) = 1000;

AK1P.SCALE(T) = 100;

AK1B.SCALE = 1000;

AK2P.SCALE(T) = 0.01;

AK2A.SCALE = 100;

AK2B.SCALE = 1000;

AK2C.SCALE = 10;

AKW.SCALE(T) = 10000;

ALPHAS.SCALE(T) = 0.01;

CONCO2INT.SCALE(T) = 100;

DIC.SCALE(T) = 0.001;

CO2AQ.SCALE(T) = 0.000001;

HCO3.SCALE(T) = 0.001;

CO32.SCALE(T) = 0.0001;

FRCO2AQ.SCALE(T) = 0.001;

COC1ALL.SCALE(T) = 1000;

COC1ATMBAC.SCALE = 100;

COC1MIXBAC.SCALE = 100;

CTR.SCALE(T,NTR) = 100;

CTR1750.SCALE(NTR) = 100;

NPPBAC.SCALE(T,NTR) = 10;

CO2FRT.SCALE(T,NTR) = 0.01;

HTRRES.SCALE(T,NTR) = 100;

FLXCOC1.SCALE(T) = 0.1;

FLXCOC2.SCALE(T) = 0.1;

FLXCOC3.SCALE(T) = 0.1;

FLXCOC4.SCALE(T) = 0.1;

FLXCTR.SCALE(T,NTR) = 0.01;

SUMCTR.SCALE(T) = 1000;

NATCH4.SCALE = 100;

TAUN2O.SCALE = 100;

RAD.SCALE(T,FOA) = 0.01;

RADTTL.SCALE(T) = 0.01;

RADHAL.SCALE(T) = 0.01;

RADODS.SCALE(T) = 0.01;

RADMIS.SCALE(T) = 0.01;

RADTTLPRI.SCALE(T) = 0.01;

NORRAD.SCALE(T) = 0.01;

RPCTT.SCALE(T,NPCR) = 100;

PCT.SCALE(T) = 0.1;

PCTRES.SCALE(T,NPCR) = 0.01;

PCS.SCALE(T) = 0.01;

PCSRES.SCALE(T,NPCR) = 0.01;

RADANTL.SCALE(T) = 0.01;

RADANTO.SCALE(T) = 0.01;

RADSOLL.SCALE(T) = 0.01;

RADSOLO.SCALE(T) = 0.01;

RADTTLL.SCALE(T) = 0.01;

RADTTLO.SCALE(T) = 0.01;

RADNUML.SCALE(T) = 0.01;

RADNUMO.SCALE(T) = 0.01;

TSS.SCALE(T) = 0.01;

TAL.SCALE(T) = 0.01;

GMT.SCALE(T) = 0.01;

T2M.SCALE(T) = 0.01;

TSUM.SCALE(T) = 0.01;

KER.SCALE(I) = 0.01;

KER1.SCALE(I) = 0.01;

KDIF.SCALE = 1000;
