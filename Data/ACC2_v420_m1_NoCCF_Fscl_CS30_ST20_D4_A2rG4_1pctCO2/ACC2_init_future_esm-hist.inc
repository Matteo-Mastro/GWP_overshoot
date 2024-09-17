********** ACC2 version 4.2: ACC2_init_future.inc **********************************************************
********** Initial conditions (future mode) ================================================================

* In order to allow a warm start of the future run starting in the year 2000,
* initial values in 2000 and best estimates of the parameters (contained in ACC2_init_future.inc)
* must be calculated by a preceding past run.
* Therefore, when parameter values or equations are modified, the past mode have to be run again.
* A recalibrated ACC2_init_future.inc file is produced automatically by invoking the past mode


*** Best estimates of the parameters -----------------------------------------------------------------------

* Preindustrial mixed-layer temperature
 TMIX1750.FX                  =   19.5900000000;

* Atmosphere-mixed layer temperature scaling factor
 DMPTMIX.FX                   =    0.5000000000;

* Beta factor for CO2 fertilization
 BETA.FX('2000')              =    0.5239643167;

* Q10 for heterotrophic respiration
 Q10BAS.FX                    =    2.0000000000;

* Temperature sensitivity to Q10
* (not used in the standard inversion.)
 Q10TMP.FX                    =    0.0460000000;

* Temperature sensitivity to NPP
* (not used in the standard inversion.)
 TNPPSENS.FX                  =    0.0160000000;

* Preindustrial ocean CO2 uptake
* CO2OCNUPPRE.FX('2000')       =   -0.3452340437;

* Preindustrial land CO2 uptake
* CO2LNDUPPRE.FX('2000')       =    0.1355717675;

* Natural methane emission
* NATCH4.FX                    =  319.8896521064;

* Natural nitrous oxide emission
* NATN2O.FX                    =   11.3359325131;

* CH4 lifetime with respect to OH depletion
* TAUCH4OH.FX                  =    8.5351401533;

* N2O lifetime
* TAUN2O.FX                    =  114.4174817663;

* Radiative forcing of a doubling of atmospheric CO2 concentration
* (fixed in the standard inversion)
 RAD2XCO2.FX                  =    3.7000000000;

* Forcing scaling factor
* (fixed in the standard inversion)
 SCLRAD.FX                    =    0.9990283568;

* Climate sensitivity
 SENS.FX('2000')              =    3.0000000000;

* Ocean vertical diffusivity
* (fixed in the standard inversion)
 KAPPA.FX                     =    0.5500000000;

* Propagator for forcing residual AR(1) model
* (not used in the standard inversion)
* RADMISAR1PRO.FX              =    0.5000000000;

* Propagator for temperature residual AR(1) model
* (not used in the standard inversion)
* T2MAR1PRO.FX                 =    0.4500000000;

