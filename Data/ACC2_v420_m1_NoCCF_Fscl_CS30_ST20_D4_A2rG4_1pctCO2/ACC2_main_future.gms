*$offlisting
$ontext

     Aggregated Carbon Cycle, Atmospheric Chemistry, and Climate Model (ACC2)

                               Version 4.2

                 Future Mode (between year 2000 and 2100 or beyond)


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
Hamburg Universit�t. International Max Planck Research School on Earth System Modelling, Hamburg, Germany. 296 pp.
http://www.sub.uni-hamburg.de/opus/volltexte/2008/3654/
(all the references cited in the model code are found)

$offtext

*** Set declarations and definitions -----------------------------------------------------------------------------------

SET T          Time steps for future mode
               /2000*2500/;
* If you run the model beyond the year 3000, the member for set I (in file ACC2_set_common.inc) needs to be expanded.

SET TOUT(T)    Time steps for future mode (used for table output)
               /2000,2005,2010,2015,2020,2025,2030,2035,2040,2045,2050
                2055,2060,2065,2070,2075,2080,2085,2090,2095,2100/;

SET SCN        Emission scenarios
               /A1B,A1T,A1FI,A2,B1,B2,IS92a,A2CBL,A2C45,B1CBL,B1C45,B2CBL,B2C45,B2CBC,B2C40,B2SBL,B2SCL,B2SMF
                A2rGB,A2rG9,A2rG6,A2rG4,1pctCO2,1pctCO2-4xext,1pctCO2-cdr,abrupt-4xCO2,abrupt-2xCO2,abrupt-0p5xCO2/;
                
SET SCN2       Emission scenarios
                /CO2-750PgC, CO2-1000PgC, CO2-2000PgC/;

* A1B:   SRES A1B
* A1T:   SRES A1T
* A1FI:  SRES A1FI
* A2:    SRES A2
* B1:    SRES B1
* B2:    SRES B2
*
* A2CBL: MESSAGE A2 baseline
* A2C45: MESSAGE A2 450 ppm in CO2-eq stabilization
* B1CBL: MESSAGE B1 baseline
* B1C45: MESSAGE B1 450 ppm in CO2-eq stabilization
* B2CBL: MESSAGE B2 baseline
* B2C45: MESSAGE B2 450 ppm in CO2-eq stabilization
* B2CBC: MESSAGE B2 450 ppm in CO2-eq stabilization without BioCCS (corresponds to 50% CO2-eq emission cut by 2050)
* B2C40: MESSAGE B2 400 ppm in CO2-eq stabilization
*
* B2SBL: MESSAGE-GAINS B2r baseline
* B2SCL: MESSAGE-GAINS B2r Current Legistration (CLE)
* B2SMF: MESSAGE-GAINS B2r Maximum Feasible Reduction (MFR)
*
* A2rGB: GGI A2r baseline
* A2rG9: GGI A2r 970 ppm in CO2-eq stabilization
* A2rG6: GGI A2r 670 ppm in CO2-eq stabilization
* A2rG4: GGI A2r 480 ppm in CO2-eq stabilization


*** Variables declarations and definitions -----------------------------------------------------------------------------

VARIABLES

CF             Cost function;


*** INCLUDE section (pre-run) ------------------------------------------------------------------------------------------

$INCLUDE ACC2_set_common.inc
$INCLUDE ACC2_constants_common.inc
$INCLUDE ACC2_variables_common.inc
$INCLUDE ACC2_scenarios_future.inc
$INCLUDE ACC2_settings_future.inc
*$INCLUDE ACC2_bounds_future.inc
$INCLUDE ACC2_switchboard_common.inc
$INCLUDE ACC2_init_future.inc
$INCLUDE ACC2_equations_common.inc
$INCLUDE ACC2_abatement_future.inc


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
* Use below when the prescribed scenario are used as it is. CF is dummy in this case but formally required in GAMS.
*CFQ                .. CF =E= 0;

* Use below when the prescribed scenario is modified to allow emission abatement
CFQ                .. CF =E= ABCTTL;

*** Assignment of emission scenario ------------------------------------------------------------------------------------

* For a single scenario run (for scenarios given till 2100 (e.g. SRES))
*ENG.FX(T,FOB) = (ENGSCN(T,FOB,'A2rGB'))$(ORD(T) LT 101)+(ENGSCN('2100',FOB,'A2rGB'))$(ORD(T) GE 101);

*NEG.FX(T,FOB) = (NEGSCN(T,FOB,'A2rGB'))$(ORD(T) LT 101)+(NEGSCN('2100',FOB,'A2rGB'))$(ORD(T) GE 101);

* For a single scenario run (for scenarios given till 2200 (e.g. A222N))
*ENG.FX(T,FOB)$(ORD(FOB) LT 4) = (ENGSCN(T,FOB,'A222N'))$((ORD(T) LT 201) and (ORD(FOB) LT 4))
*                                +(ENGSCN('2200',FOB,'A222N'))$((ORD(T) GE 201) and (ORD(FOB) LT 4));
*
*NEG.FX(T,FOB)$(ORD(FOB) LT 4) = (NEGSCN(T,FOB,'A222N'))$((ORD(T) LT 201) and (ORD(FOB) LT 4))
*                                +(NEGSCN('2200',FOB,'A222N'))$((ORD(T) GE 201) and (ORD(FOB) LT 4));
*
*ENG.FX(T,FOB)$(ORD(FOB) GE 4) = (ENGSCN(T,FOB,'A222N'))$((ORD(T) LT 101) and (ORD(FOB) GE 4))
*                                +(ENGSCN('2100',FOB,'A222N'))$((ORD(T) GE 101) and (ORD(FOB) GE 4));
*
*NEG.FX(T,FOB)$(ORD(FOB) GE 4) = (NEGSCN(T,FOB,'A222N'))$((ORD(T) LT 101) and (ORD(FOB) GE 4))
*                                +(NEGSCN('2100',FOB,'A222N'))$((ORD(T) GE 101) and (ORD(FOB) GE 4));


*** Model definition ---------------------------------------------------------------------------------------------------

MODEL FUTURE  "ACC2 Future Mode" /All/;


*** Execution and associated options -----------------------------------------------------------------------------------

FUTURE.OPTFILE = 1;

FUTURE.SCALEOPT = 1;

OPTION NLP = CONOPT3;

OPTION PROFILE = 1;

OPTION PROFILETOL = 1;

OPTION LIMCOL = 0;

OPTION LIMROW = 0;

*OPTION SOLPRINT = off;

OPTION DMPSYM;

OPTION RESLIM = 500000;
* To specify the time limit for model run (in seconds)
* (e.g. 500000 seconds is approximately 138 hours)

OPTION ITERLIM = 1000000;
* To specify the limit for the number of iterations

FUTURE.SAVEPOINT = 1;

EXECUTE_LOADPOINT 'FUTURE_p.gdx';

SOLVE FUTURE USING NLP MINIMIZING CF;


* Marginal abatement costs alternatively obtained from the marginals of the emission equations
MACALTREC(T,'CO2') = ANTQ1.M(T)*DSCCUR.L(T);
MACALTREC(T,'CH4') = ANTQ2.M(T)*DSCCUR.L(T)*1000;
MACALTREC(T,'N2O') = ANTQ3.M(T)*DSCCUR.L(T)*1000*(28/44);

* GCP (price ratio) for CH4 and N2O
GCPCH4REC(T)$(MACALTREC(T,'CO2') > 0) = MACALTREC(T,'CH4')/MACALTREC(T,'CO2')/(12/44);
GCPCH4REC(T)$(MACALTREC(T,'CO2') = 0) = 0;
GCPN2OREC(T)$(MACALTREC(T,'CO2') > 0) = MACALTREC(T,'N2O')/MACALTREC(T,'CO2')/(12/44);
GCPN2OREC(T)$(MACALTREC(T,'CO2') = 0) = 0;
* See the clarification about the emission units in the file "ACC2_abatement_future.inc"

ABECO2EQREC(T) = ENG.L(T,'CO2')*(ABL.L(T,'CO2')/100)
                 +(ENG.L(T,'CH4')+NEG.L(T,'CH4'))*(ABL.L(T,'CH4')/100)/1000*(12/44)*GCPCH4REC(T)
                 +(ENG.L(T,'N2O')+NEG.L(T,'N2O'))*(ABL.L(T,'N2O')/100)/1000/(28/44)*(12/44)*GCPN2OREC(T);

* For multiple scenario runs
*LOOP (SCN,
*
*ENG.FX(T,FOB)  = (ENGSCN(T,FOB,SCN))$(ORD(T) LT 101)+(ENGSCN('2100',FOB,SCN))$(ORD(T) GE 101);
*
*NEG.FX(T,FOB)  = (NEGSCN(T,FOB,SCN))$(ORD(T) LT 101)+(NEGSCN('2100',FOB,SCN))$(ORD(T) GE 101);
*
*SOLVE CLIMATE USING DNLP MINIMIZING GF;
*
*LOOP (T,
*SAVSUMR(T,SCN)=RADTTL.L(T);
*SAVSUMT(T,SCN)=T2MRT1990.L(T););
*
*Display SAVSUMR;
*Display SAVSUMT;);



* Export the results in gdx file to csv

*$ontext
File results_future_scn / results_future_ZECMIP-750PgC_standard.csv /;
results_future_scn.AP=0;
results_future_scn.PW=1000;
put results_future_scn;
PUT
'YEAR ',
'ANTCO2 ',
'ANTCH4 ',
'ANTN2O ',
'ANTCO2EQ ',
'ABLCO2 ',
'ABLCH4 ',
'ABLN2O ',
'MACCO2 ',
'MACCH4 ',
'MACN2O ',
'ABCCO2 ',
'ABCCH4 ',
'ABCN2O ',
'CONCO2 ',
'CONCH4 ',
'CONN2O ',
'CONrOH ',
'RADCO2 ',
'RADCH4 ',
'RADN2O ',
'RADO3TRP ',
'RADHALO ',
'RADSF6 ',
'RADSO2 ',
'RADCRBAER ',
'RADAERIND ',
'RADTTL ',
'T2M ',
'pH ',
'OCNCO2UP ',
'LNDCO2UP ',
'DSCCUR ',
'HCIO ',
'HCML ',
'HCAL '/;
LOOP(T, PUT T.TL:4,
ANT.L(T,'CO2'):10:4,
ANT.L(T,'CH4'):10:2,
ANT.L(T,'N2O'):8:2,
ANTCO2EQ.L(T):8:2,
ABL.L(T,'CO2'):8:2,
ABL.L(T,'CH4'):8:2,
ABL.L(T,'N2O'):8:2,
MAC.L(T,'CO2'):15:2,
MAC.L(T,'CH4'):15:2,
MAC.L(T,'N2O'):15:2,
ABC.L(T,'CO2'):15:2,
ABC.L(T,'CH4'):15:2,
ABC.L(T,'N2O'):15:2,
CON.L(T,'CO2'):9:2,
CON.L(T,'CH4'):9:2,
CON.L(T,'N2O'):9:2,
CON.L(T,'rOH'):9:4,
RAD.L(T,'CO2'):9:4,
RAD.L(T,'CH4'):9:4,
RAD.L(T,'N2O'):9:4,
RAD.L(T,'O3TRP'):9:4,
RADHAL.L(T):9:4,
RAD.L(T,'SF6'):9:4,
RAD.L(T,'SO2'):9:4,
RAD.L(T,'CRBAER'):9:4,
RAD.L(T,'AERIND'):9:4,
RADTTL.L(T):9:4,
T2M.L(T):9:4,
PH.L(T):9:4,
CO2OCNUPTTL.L(T):12:6,
CO2LNDUPTTL.L(T):12:6,
DSCCUR.L(T):13:2,
HCIO.L(T):20:6,
HCML.L(T):20:6,
HCAL.L(T):20:6/;);
*$offtext
*** INCLUDE section (post-run) -----------------------------------------------------------------------------------------

*$INCLUDE ACC2_outputs_future.inc
$INCLUDE ACC2_alert_future.inc
*$INCLUDE ACC2_transfer_metric_future.inc
