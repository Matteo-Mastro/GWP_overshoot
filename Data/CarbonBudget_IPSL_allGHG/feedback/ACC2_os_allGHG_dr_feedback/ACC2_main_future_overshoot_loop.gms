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
                A2rGB,A2rG9,A2rG6,A2rG4
                rcp60,ssp370,ssp370-lowNTCF,ssp434,ssp460,rcp26,ssp119,ssp126,rcp85,ssp245,rcp45,ssp534-over,ssp585/;

SET T2MDIM /15os/;
*SET T2MDIM /15os,2os/;

SET ENDDIM /2120,2130,2140,2150,2160,2170,2180,2190/;

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

*PARAMETER SETENDDIM(ENDDIM) /2070 71
*                             2080 81
*                             2090 91
*                             2100 101
*                             2110 111
*                             2120 121
*                             2130 131
*                             2140 141
*                             2150 151
*                             2160 161
*                             2170 171
*                             2180 181
*                             2190 191/;
  
PARAMETER SETENDDIM(ENDDIM) /2190 191
                             2180 181
                             2170 171
                             2160 161
                             2150 151
                             2140 141
                             2130 131
                             2120 121
/;
                            
PARAMETER SETENDT2MDIM(T2MDIM) /15os 1/;
*                                20os  1/;

Parameter SETT2MDIM(T2MDIM) /15os 1.5/;
*                             20os  2/;

PARAMETERS

SCNT2MVALUE
SETENDDIMVALUE
SETT2MDIMVALUE
SETENDT2MDIMVALUE

T2M1(T,T2MDIM,ENDDIM)
ANTCO21(T,T2MDIM,ENDDIM)
ANTCH41(T,T2MDIM,ENDDIM)
;

VARIABLES

T2M2(T,T2MDIM,ENDDIM)
ANTCO22(T,T2MDIM,ENDDIM)
ANTCH42(T,T2MDIM,ENDDIM)
;


*** INCLUDE section (pre-run) ------------------------------------------------------------------------------------------

$INCLUDE ACC2_set_common.inc
$INCLUDE ACC2_constants_common.inc
$INCLUDE ACC2_variables_common.inc
$INCLUDE ACC2_scenarios_future_2300.inc
$INCLUDE ACC2_settings_future.inc
$INCLUDE ACC2_bounds_future.inc
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

LOOP((T2MDIM,ENDDIM),

EXECUTE_LOADPOINT 'FUTURE_p.gdx';

SOLVE FUTURE USING NLP MINIMIZING CF;


ANTCO21(T,T2MDIM,ENDDIM) = ANT.L(T,'CO2');
ANTCH41(T,T2MDIM,ENDDIM) = ANT.L(T,'CH4');
T2M1(T,T2MDIM,ENDDIM) = T2M.L(T);


SETENDT2MDIMVALUE=SETENDT2MDIM(T2MDIM);
display SETENDT2MDIMVALUE;
SETENDDIMVALUE=SETENDDIM(ENDDIM);
display SETENDDIMVALUE;
SETT2MDIMVALUE=SETT2MDIM(T2MDIM);
display SETT2MDIMVALUE;

*T2M.UP(T)$(ORD(T) GE SETENDDIMVALUE*SETENDT2MDIMVALUE) = SETT2MDIMVALUE;
T2M.UP(T)$(ORD(T) GE SETENDDIMVALUE) = SETT2MDIMVALUE;


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

);

*** Data transfer from GDX to EXCEL (post-run) ------------------------------------------------------------------------
* See the GAMS Documentation section: Writing Parameter to Spreadsheet including Zero Values
* Moving data from parameters storing the loop results to variables
ANTCO22.L(T,T2MDIM,ENDDIM) = ANTCO21(T,T2MDIM,ENDDIM) ;
ANTCH42.L(T,T2MDIM,ENDDIM) = ANTCH41(T,T2MDIM,ENDDIM) ;
T2M2.L(T,T2MDIM,ENDDIM) = T2M1(T,T2MDIM,ENDDIM);

* In order to write every entry including zeros in the spreadsheet, one must allocate
* a non-zero value to one of the variable attributes .m, .lo or .up
ANTCO22.M(T,T2MDIM,ENDDIM) = 1;
ANTCH42.M(T,T2MDIM,ENDDIM) = 1;
T2M2.M(T,T2MDIM,ENDDIM) = 1;

* Storing the loop results kept in variables in gdx file
execute_unload 'out_future_loop.gdx' ANTCO22,ANTCH42,T2M2;

* Exporting the loop results in gdx file to excel
$onecho > ACC2_excel_future.txt

var = ANTCO22.L RDIM=1 CDIM=3 RNG=ANTCO2!
var = ANTCH42.L RDIM=1 CDIM=3 RNG=ANTCH4!
var = T2M.L RDIM=1 CDIM=3 RNG=T2M!

$offecho

execute 'gdxxrw.exe i=out_future_loop.gdx o=out_future_loop_CO2_15.xlsx squeeze=n @ACC2_excel_future.txt';

*execute 'o=out_future_loop_overshoot.xlsx ';
*o=out_future_loop_overshoot.xlsx squeeze=n @ACC2_excel_future.txt';

*** INCLUDE section (post-run) -----------------------------------------------------------------------------------------

$INCLUDE ACC2_outputs_future.inc
$INCLUDE ACC2_alert_future.inc
*$INCLUDE ACC2_transfer_metric_future.inc
