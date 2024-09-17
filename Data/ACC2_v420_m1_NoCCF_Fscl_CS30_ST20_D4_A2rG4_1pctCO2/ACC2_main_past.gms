$offlisting
$ontext

     Aggregated Carbon Cycle, Atmospheric Chemistry, and Climate Model (ACC2)

                               Version 4.2

                 Past Mode (between year 1750 and 2000)


Release date:
23 March 2017

Developers:
Katsumasa Tanaka1, Elmar Kriegler2

1:National Institute for Environmental Studies (NIES), Japan
2:Potsdam Institute for Climate Impact Research (PIK), Germany

For any questions and suggestions, please contact Katsumasa Tanaka (tanaka.katsumasa@nies.go.jp).


                               (Copyright 2017)


                               History of model development

ACC2 was developed from Integrated Assessment of Climate Protection Strategies (ICLIPS)
Climate Model (ICM) version 1.1 (Bruckner et al., 2003) as a starting point.
The origin of ACC2 and ICM is Nonlinear Impulse-response representation of the coupled
Carbon cycle-Climate System (NICCS) (Hooss, 2001; Hooss et al., 2001).
Diffusion Ocean Energy balance CLIMate model (DOECLIM) (Kriegler, 2005) is implemented to
ACC2 to calculate temperature from radiative forcing.


                               Note for use

ACC2 version 4.1 was developed on GAMS 24.8.3 Windows 64-bit version (released in 28 January 2017).
The model may behave differently in other versions of GAMS or platforms.

Except for changes in the future emission scenarios, upon any changes to the model codes,
one must run the past mode again by executing ACC2_main_past.gms.
The past mode overwrites ACC2_init_future.inc,
which contains the estimates of time-dependent variables in 2000 and other variables
(including uncertain parameters for the inverse calculation).
Then ACC2_main_future.gms can be executed to run the simulation from 2000.

Note that the sea level component is not described in Tanaka (2008).
The sea level component is developed based on IPCC (2001, Appendix 11.1).


                               File structure

See Section D.3 in Tanaka (2008) (with some changes).


                               Acknowledgments

Financial supports:
Max Planck Institute for Biogeochemistry, Jena
Max Planck Institute for Meteorology, Hamburg
University of Hamburg
International Institute for Applied Systems Analysis (IIASA)


                               Main reference

Tanaka, K. (2008) Inverse estimation for the simple Earth system model ACC2 and its applications. Ph.D. dissertation.
Hamburg Universität. International Max Planck Research School on Earth System Modelling, Hamburg, Germany. 296 pp.
http://www.sub.uni-hamburg.de/opus/volltexte/2008/3654/
(all the references cited in the model code are found)

$offtext

*** Set declarations and definitions -----------------------------------------------------------------------------------

SET T          Time steps for past mode
               /1750*2000/;

SET SCN        Scenarios
               /1pctCO2,1pctCO2-4xext,1pctCO2-cdr,abrupt-4xCO2,abrupt-2xCO2,abrupt-0p5xCO2/;


*** Variables declarations and definitions -----------------------------------------------------------------------------

VARIABLES

CF             Cost function;


*** INCLUDE section (pre-run) ------------------------------------------------------------------------------------------

$INCLUDE ACC2_set_common.inc
$INCLUDE ACC2_constants_common.inc
$INCLUDE ACC2_variables_common.inc
$INCLUDE ACC2_scenarios_past.inc
$INCLUDE ACC2_settings_past.inc
*$INCLUDE ACC2_bounds_past.inc
$INCLUDE ACC2_switchboard_common.inc
$INCLUDE ACC2_init_past.inc
$INCLUDE ACC2_equations_common.inc
$INCLUDE ACC2_init_future_esm-hist.inc


*** Import from the interface ------------------------------------------------------------------------------------------

*$CALL 'GDXXRW INPUT=ACC2_switchboard_common.xls INDEX=BackgroundPre!a1'
*$GDXIN ACC2_switchboard_common.gdx
*$LOAD CO2FFSWCH CO2LUSWCH ANTCH4N2OSWCH RADMISSWCH TMIX1750SWCH DMPTMIXSWCH BETASWCH Q10BASSWCH Q10TMPSWCH TNPPSENSSWCH
*$LOAD CO2OCNUPPRESWCH CO2LNDUPPRESWCH NATCH4N2OSWCH TAUCH4N2OSWCH RAD2XCO2SWCH SCLRADSWCH SENSSWCH KAPPASWCH
*$LOAD RADMISAR1PROSWCH T2MAR1PROSWCH CO2UPANTSWCH CONCO2SWCH CONCH4N2OSWCH T2MSWCH DOECLIMSWCH DYNCARBSWCH OCNCO2TFSWCH
*$LOAD LNDCO2TFNPPSWCH LNDCO2TFRESSWCH LNDCO2TFQ10SWCH CO2UPPRESWCH CO2LUSGMSWCH CONCO2SGMSWCH ANTCH4N2OSGMSWCH
*$LOAD CONCH4SGMSWCH CONN2OSGMSWCH ENSOCO2SWCH ENSOT2MSWCH RADMISAR1SWCH T2MAR1SWCH VOLSWCH MISTRNSWCH SOLSWCH
*$GDXIN


*** Equation declarations and definitions ------------------------------------------------------------------------------

EQUATIONS

CFQ;

* Cost function
CFQ              .. CF =E=
                         CO2FFRESSUM
                         +CO2LURESSUM
                         +ANTRESSUMCH4
                         +ANTRESSUMN2O
                         +RADMISRESSUM

                         +TMIX1750RESSUM
                         +DMPTMIXRESSUM
                         +BETARESSUM
                         +Q10BASRESSUM
                         +Q10TMPRESSUM
                         +TNPPSENSRESSUM
                         +CO2OCNUPPRERESSUM
                         +CO2LNDUPPRERESSUM
                         +NATCH4RESSUM
                         +NATN2ORESSUM
                         +TAUCH4OHRESSUM
                         +TAUN2ORESSUM
                         +RAD2XCO2RESSUM
                         +SCLRADRESSUM
                         +SENSRESSUM
                         +KAPPARESSUM
                         +RADMISAR1PRORESSUM
                         +T2MAR1PRORESSUM

                         +CO2OCNUPANTRESSUM
                         +CO2LNDUPANTRESSUM
                         +CONCO2RESSUM
                         +CONCH4RESSUM
                         +CONN2ORESSUM
                         +T2MRESSUM;

* Q10TMPRESSUM, TNPPSENSRESSUM, RAD2XCO2RESSUM, SCLRADRESSUM, KAPPARESSUM, RADMISAR1PRORESSUM and T2MAR1PRORESSUM are
* fixed at zero in the standard inversion (coupled; missing forcing approach). See ACC2_switchboard_common.xls.


*** Model definition ---------------------------------------------------------------------------------------------------

MODEL PAST "ACC2 Past Mode" /All/;


*** Execution and associated options -----------------------------------------------------------------------------------

PAST.OPTFILE = 1;
* To let the solver use the option file (conopt.opt and conopt3.opt)
* The name of the options file should be the same with the name of the solver.
* If you have used option NLP=CONOPT, then the options file should be
* conopt.opt. (Arne Stolbjerg Drud, personal communication, Oct. 25, 2006)
* In case the option file is not recognized (no explicit error is shown!),
* two identical option files (conopt.opt and conopt3.opt) are still prepared.

PAST.SCALEOPT = 1;
* To activate both automatic and manual scalings for variables and equations

OPTION NLP = CONOPT3;
* To solve the inverse calculation by using CONOPT3
* The convention is that NLP=CONOPT will give you the latest version.
* If you would like to use another version, you must explicitly use
* NLP=CONOPT1 / CONOPT2 / CONOPT3.
* (Arne Stolbjerg Drud, personal communication, Oct. 25, 2006)

*OPTION NLP = GAMSCHK;
* To perform a detailed check on the model structure
* Note that the conopt solver has to be replaced with GAMSCHK to perform that.

OPTION SOLSLACK = 1;
* To generate an equation output that contains slack variables.

OPTION PROFILE = 0;
* To generate more information on program execution profiles

OPTION LIMCOL = 0;
* To control the number of columns that are listed for each variable
* in the COLUMN LISTING section of the listing file.
* The assignment of zero means suppressions of the equation listing.

OPTION LIMROW = 0;
* To control the number of rows that are listed for each equation
* in the EQUATION LISTING section of the listing file.
* The assignment of zero means suppressions of the variable listing.

OPTION SOLPRINT = off;
* Not to print the detailed solution information for equations and variables.
* Recommended for series of optimizations.

OPTION DMPSYM;
* To diagnose memory problems

OPTION RESLIM = 500000;
* To specify the time limit for model run (in seconds)
* (e.g. 500000 seconds is approximately 138 hours)

OPTION ITERLIM = 1000000;
* To specify the limit for the number of iterations

* GAMSBAS is no longer recommended as it does not work with all the variables and equations
* (Arne Stolbjerg Drud, personal communication, February 15, 2008). Instead, we use SAVEPOINT/LOADPOINT.
* OPTION NLP = GAMSBAS;
* To save a file that specifies bases
* $include ACC2_main_past.bas
* PAST.bratio=0;
* As suggested by Franz Nelissen (personal communication, October 7, 2005).

PAST.SAVEPOINT = 1;
* To save the current solution in the GDX file (PAST_p.gdx).

EXECUTE_LOADPOINT 'PAST_p.gdx';
* To load the previous solution from the GDX file (PAST_p.dgx).

* To solve a single optimization, use the following solve statement:
SOLVE PAST USING NLP MINIMIZING CF;
* This model is classified into NLP (Franz Nelissen at GAMS software GmbH, personal communication on Oct. 5, 2005).
* Use of DNLP is generally not recommended according to the GAMS User Guide.

$ontext
* To solve a series of optimizations, use the solve statement below:
$INCLUDE ACC2_series_past.inc
$offtext


*** Export to the interface --------------------------------------------------------------------------------------------
*
*EXECUTE_UNLOAD 'ACC2_switchboard_common.gdx', CF, CO2FFRESSUM, CO2LURESSUM, ANTRESSUMCH4, ANTRESSUMN2O,
*RADMISRESSUM, TMIX1750RESSUM, DMPTMIXRESSUM, BETARESSUM, Q10BASRESSUM, Q10TMPRESSUM, TNPPSENSRESSUM, CO2OCNUPPRERESSUM,
*CO2LNDUPPRERESSUM, NATCH4RESSUM, NATN2ORESSUM, TAUCH4OHRESSUM, TAUN2ORESSUM, RAD2XCO2RESSUM, SCLRADRESSUM, SENSRESSUM,
*KAPPARESSUM, RADMISAR1PRORESSUM, T2MAR1PRORESSUM, CO2OCNUPANTRESSUM, CO2LNDUPANTRESSUM, CONCO2RESSUM, CONCH4RESSUM,
*CONN2ORESSUM, T2MRESSUM, TMIX1750, DMPTMIX, BETA, Q10BAS, Q10TMP, TNPPSENS, CO2OCNUPPRE, CO2LNDUPPRE, NATCH4, NATN2O,
*TAUCH4OH, TAUN2O, RAD2XCO2, SCLRAD, SENS, KAPPA, RADMISAR1PRO, T2MAR1PRO, RADAERAVE, RADMISAVE, YINTRCO2ENSO,
*SLOPECO2ENSO, RSQ1, RSQ2;
*EXECUTE 'GDXXRW.EXE INPUT=ACC2_switchboard_common.gdx OUTPUT=ACC2_switchboard_common.xls INDEX=BackgroundPost!a1'


*** INCLUDE section (post-run) -----------------------------------------------------------------------------------------

$INCLUDE ACC2_alert_past.inc
$INCLUDE ACC2_transfer_past.inc
$INCLUDE ACC2_outputs_past.inc
