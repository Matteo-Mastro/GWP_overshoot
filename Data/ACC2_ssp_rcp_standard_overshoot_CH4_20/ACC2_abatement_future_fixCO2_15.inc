********** ACC2 version 4.2: ACC2_abatement_future.inc *****************************************************************
*
****** Stabilization target ============================================================================================

$if not setglobal TARGET $setglobal TARGET 101
$if not setglobal TEMP $setglobal TEMP 1.5
T2M.UP(T)$(Ord(T) GE %TARGET% ) = %TEMP%;

* The new CO2 MAC should only be applied to fossil fuel CO2 emissions and not deforestation CO2. However, the CH4 and
* N2O should be applied to all anthropogenic emissions (Daniel Johansson, personal communication, 12 June 2013).
* Choose one from below:
* A2rGB: GGI A2r baseline
* A2rG9: GGI A2r 970 ppm in CO2-eq stabilization
* A2rG6: GGI A2r 670 ppm in CO2-eq stabilization
* A2rG4: GGI A2r 480 ppm in CO2-eq stabilization

* For CO2
ENG.FX(T,FOB)$(Ord(FOB) EQ 1) = (ENGSCN(T,FOB,'ssp119'))$(Ord(T) LT 101)+(ENGSCN('2100',FOB,'ssp119'))$((Ord(T) GE 101) AND (Ord(T) LT 201)) + 0$(Ord(T) GE 201);

NEG.FX(T,FOB)$(Ord(FOB) EQ 1) = (NEGSCN(T,FOB,'ssp119'))$(Ord(T) LT 101)+(NEGSCN('2100',FOB,'ssp119'))$((Ord(T) GE 101) AND (Ord(T) LT 201)) + 0$(Ord(T) GE 201);

* For CH4
ENG.FX(T,FOB)$(Ord(FOB) EQ 2) = (ENGSCN(T,FOB,'A2rGB'))$(Ord(T) LT 101)+(ENGSCN('2100',FOB,'A2rGB'))$(Ord(T) GE 101);

NEG.FX(T,FOB)$(Ord(FOB) EQ 2) = (NEGSCN(T,FOB,'A2rGB'))$(Ord(T) LT 101)+(NEGSCN('2100',FOB,'A2rGB'))$(Ord(T) GE 101);

* For N2O
ENG.FX(T,FOB)$(Ord(FOB) EQ 3) = (ENGSCN(T,FOB,'ssp119'))$(Ord(T) LT 101)+(ENGSCN('2100',FOB,'ssp119'))$(Ord(T) GE 101);

NEG.FX(T,FOB)$(Ord(FOB) EQ 3) = (NEGSCN(T,FOB,'ssp119'))$(Ord(T) LT 101)+(NEGSCN('2100',FOB,'ssp119'))$(Ord(T) GE 101);

* For agents other than CO2, CH4, and N2O
ENG.FX(T,FOB)$(Ord(FOB) GT 3) = (ENGSCN(T,FOB,'ssp119'))$(Ord(T) LT 101)+(ENGSCN('2100',FOB,'ssp119'))$(Ord(T) GE 101);

NEG.FX(T,FOB)$(Ord(FOB) GT 3) = (NEGSCN(T,FOB,'ssp119'))$(Ord(T) LT 101)+(NEGSCN('2100',FOB,'ssp119'))$(Ord(T) GE 101);


*ENG.FX(T,FOB) = 0;

*NEG.FX(T,FOB) = 0;

* For all agents
*ENG.FX(T,FOB) = ENGSCN(T,FOB,'ssp460');

*NEG.FX(T,FOB) = 0;


****** Constant declarations and assignments ===========================================================================

SCALARS

DSC                      Discount rate                                          [%]  /4/
* 2% Stern (2006)
* 4% Johansson (2011, Climatic Change) reference case
* 6% Nordhaus (2008)

IFL                      Adjustment for inflation between 2000 and 2010              /0.7897/
* The data from the GET model used to derive the CO2 MAC function is based on the US dollar in 2010. Other MAC functions
* have a reference in 2000. Thus, the dollar unit for the other MAC function needs adjustment.
* The inflation estimate is obtained on 13 June 2013 from http://www.bls.gov/data/inflation_calculator.htm

ZEY                      Year to meet the net zero GHG emission target;


PARAMETERS

MACA(FOH)                Constant a in the marginal abatement cost function
*                        That for CO2 is derived from an output from GET by K. Tanaka (11 June 2013).
*                        Those for CH4 and N2O are taken from Johansson (2011, Climatic Change, Table 1).
                         /CO2 10.53
                          CH4 0.0
                          N2O 0.0/

MACB(FOH)                Constant b in the marginal abatement cost function
                         /CO2 0.004707
                          CH4 34.4
                          N2O 127/

MACC(FOH)                Constant c in the marginal abatement cost function
                         /CO2 0.1173
                          CH4 0.10
                          N2O 0.16/

MACALTREC(T,FOH)         Marginal abatement costs alternatively obtained from emission equations (for each agent)
*                        [US$ in year T])

GCPCH4REC(T)             GCP (price ratio) of CH4 each year

GCPN2OREC(T)             GCP (price ratio) of N2O each year

ABECO2EQREC(T)           CO2-equivalent abated emissions each year (based on GCP)

GWP100(FOA)              GWP100 taken from from AR5 (Table 8A)
                         /CO2            1
                          CH4            28
                          N2O            265
                          CF4            6630
                          C2F6           11100
                          C4F10          9200
                          HFC23          12400
                          HFC32          677
                          HFC4310mee     1650
                          HFC125         3170
                          HFC134a        1300
                          HFC143a        4800
                          HFC152a        138
                          HFC227ea       3350
                          HFC236fa       8060
                          HFC245ca       716
                          SF6            23500/

GWP020(FOA)              GWP20 taken from from AR5 (Table 8A)
                         /CO2            1
                          CH4            84
                          N2O            264
                          CF4            4880
                          C2F6           8210
                          C4F10          6870
                          HFC23          10800
                          HFC32          2430
                          HFC4310mee     4310
                          HFC125         6090
                          HFC134a        3710
                          HFC143a        6940
                          HFC152a        506
                          HFC227ea       5360
                          HFC236fa       6940
                          HFC245ca       2510
                          SF6            17500/

GTP100(FOA)              GTP100 taken from from AR5 (Table 8A)
                         /CO2            1
                          CH4            4
                          N2O            234
                          CF4            8040
                          C2F6           13500
                          C4F10          11000
                          HFC23          12700
                          HFC32          94
                          HFC4310mee     281
                          HFC125         967
                          HFC134a        201
                          HFC143a        2500
                          HFC152a        19
                          HFC227ea       1460
                          HFC236fa       8380
                          HFC245ca       100
                          SF6            28200/

GTP020(FOA)              GTP20 taken from from AR5 (Table 8A)
                         /CO2            1
                          CH4            67
                          N2O            277
                          CF4            5270
                          C2F6           8880
                          C4F10          7420
                          HFC23          11500
                          HFC32          1360
                          HFC4310mee     3720
                          HFC125         5800
                          HFC134a        3050
                          HFC143a        6960
                          HFC152a        174
                          HFC227ea       5280
                          HFC236fa       7400
                          HFC245ca       1570
                          SF6            18900/;

****** Variable declarations ===========================================================================================

VARIABLES

ABL(T,FOH)               Abatement level relative to the prescribed scenario
*                        CO2             [%]
*                        CH4             [%]
*                        N2O             [%]

ABLFIR(T,FOH)            First order derivative in the abatement level [%]

ABLSEC(T,FOH)            Second order derivative in the abatement level [%]

MAC(T,FOH)               Marginal abatement costs (for each agent)
*                        [US$ in year T]

ABC(T,FOH)               Abatement costs each year (for each agent)
*                        [Billion US$ in year T]

ABCINT(FOH)              Abatement costs integrated over the time horizon considered (for each agent)
*                        [Billion US$ in year 2000]

ABCTTL                   Total abatement costs integrated over the time horizon considered (for all the agents)
*                        [Billion US$ in year 2000]

DSCCUR(T)                The current effect of discount rate relative to present [%]

METRIC(FOA)              Emission metric used to compute the CO2eq emissions

ANTCO2EQ(T)              Abated CO2 equivalent emissions of Kyoto Gases [GtCO2eq]

BSLCO2EQ(T)              Baseline CO2 equivalent emissions of Kyoto Gases [GtCO2eq]

ANTCO2EQFIR(T)           First derivative of GHG emissions [GtCO2eq]

ANTCO2EQSEC(T)           Second derivative of GHG emissions [GtCO2eq]

ANTCO2FIR(T)             First derivative of CO2 emissions [GtC]

ANTCO2SEC(T)             Second derivative of CO2 emissions [GtC]

ANTCH4FIR(T)             First derivative of CH4 emissions [MtCH4]

ANTCH4SEC(T)             Second derivative of CH4 emissions [MtCH4]

ANTN2OFIR(T)             First derivative of N2O emissions [MtN]

ANTN2OSEC(T)             Second derivative of N2O emissions [MtN];


****** Equation declarations ===========================================================================================

EQUATIONS

ANTQ1(T)

ANTQ2(T)

ANTQ3(T)

ANTQ4(T,FOB)

ABLINIQ(T,FOH)

ABLFIRQ(T,FOH)

ABLSECQ(T,FOH)

MACQ1(T)

MACQ2(T)

MACQ3(T)

ABCQ1(T)

ABCQ2(T)

ABCQ3(T)

ABCINTQ(FOH)

ABCTTLQ

DSCCURQ1(T)

DSCCURQ2(T)

ANTCO2EQQ(T)

BSLCO2EQQ(T)

ANTCO2EQFIRQ(T)

ANTCO2EQSECQ(T)

ANTCO2FIRQ(T)

ANTCO2SECQ(T)

ANTCH4FIRQ(T)

ANTCH4SECQ(T)

ANTN2OFIRQ(T)

ANTN2OSECQ(T);


* Assigning the emission scenario with abatement for CO2, CH4, and N2O
* The new CO2 MAC should only be applied to fossil fuel CO2 emissions and not deforestation CO2. However, the CH4 and
* N2O should be applied to all anthropogenic emissions (Daniel Johansson, personal communication, 12 June 2013).
ANTQ1(T)                 .. ANT(T,'CO2') =E= ENG(T,'CO2')*(1-ABL(T,'CO2')/100)+NEG(T,'CO2');

ANTQ2(T)                 .. ANT(T,'CH4') =E= (ENG(T,'CH4')+NEG(T,'CH4'))*(1-ABL(T,'CH4')/100);

ANTQ3(T)                 .. ANT(T,'N2O') =E= (ENG(T,'N2O')+NEG(T,'N2O'))*(1-ABL(T,'N2O')/100);

ANTQ4(T,FOB)$(Ord(FOB) GE 4) .. ANT(T,FOB) =E= ENG(T,FOB)+NEG(T,FOB);

* Emission abatements are allowed only after year 2020.
ABLINIQ(T,FOH)$(Ord(T) LE 21) .. ABL(T,FOH) =E= 0;

* The maximum abatement levels for CO2, CH4, N2O, and BC are assumed 112%, 70%, 50%, and 23% relative to the baseline
* levels of each year, respectively.
* These assumptions seem fine with Daniel Johannson (31 March 2015, personal communication)

*ABL.UP(T,'CO2') = 112;
ABL.UP(T,'CH4') = 70;
*ABL.UP(T,'N2O') = 50;

ABL.FX(T,'CO2') = 0;
*ABL.FX(T,'CH4') = 0;
ABL.FX(T,'N2O') = 0;

* The change in abatement level is kept less than 4% per year.
* The following first order condition should not be used to compute metric values (e.g. price ratio)
* because long-term effects (technological change and learning by doing) should not influence the metric values.
* (Daniel Johansson, June 16, 2011, personal communication)
ABLFIRQ(T,FOH)$(Ord(T) GT 1) .. ABLFIR(T,FOH) =E= ABL(T,FOH)-ABL(T-1,FOH);

ABLFIR.UP(T,FOH) = 4;
ABLFIR.LO(T,FOH) = -4;

* It takes 10 years until the abatement level reaches 4% from the baseline,
* which is consistent with Daniel's approach.
* A smoothing constraint to avoid fluctuations in the solution
* The change in abatement level is kept less than 4% per year.
* The following second order condition should not be used to compute metric values (e.g. price ratio)
* because long-term effects (technological change and learning by doing) should not influence the metric values.
* (Daniel Johansson, June 16, 2011, personal communication)
ABLSECQ(T,FOH)$(Ord(T) GT 2) .. ABLSEC(T,FOH) =E= ABLFIR(T,FOH)-ABLFIR(T-1,FOH);

ABLSEC.UP(T,FOH) = 0.4;
ABLSEC.LO(T,FOH) = -0.4;

* Marginal abatement costs
MACQ1(T) .. MAC(T,'CO2') =E= MACA('CO2')*ABL(T,'CO2')+MACB('CO2')*exp(MACC('CO2')*ABL(T,'CO2'))-MACB('CO2');

MACQ2(T) .. MAC(T,'CH4') =E= (MACA('CH4')*ABL(T,'CH4')+MACB('CH4')*exp(MACC('CH4')*ABL(T,'CH4'))-MACB('CH4'))/IFL;

MACQ3(T) .. MAC(T,'N2O') =E= (MACA('N2O')*ABL(T,'N2O')+MACB('N2O')*exp(MACC('N2O')*ABL(T,'N2O'))-MACB('N2O'))/IFL;

* Abatement costs each year (for each agent)
* Integrating of the marginal abatement cost function with respect to the abatement amount
* (not with respect to the abatement percentage)
* Note: There are three different unit systems used to express CO2, CH4, and N2O emissions relevant to
* the price ratio calculations (confirmed with Daniel Johansson (personal correspondence, October 4, 2011)):
*        MAC functions (Johansson, 2010, Climatic Change)        C,   CH4, N2O
*        ACC2 and SRES                                           C,   CH4, N2
*        Price ratios                                            CO2, CH4, N2O
* The difference in the unit systems explains why, while there is a multiplication term of 44/28 to adjust for
* the effect of molecular weight in ABCQ3(T), there is no multiplication term of 44/12 in equation ABCQ1(T).
* This also explains why 44/12 is multiplied when the price ratio is computed in the separate Excel file.

ABCQ1(T)                 .. ABC(T,'CO2') =E=
                            1/2*MACA('CO2')*ENG(T,'CO2')/100*ABL(T,'CO2')*ABL(T,'CO2')
                            +MACB('CO2')/MACC('CO2')*ENG(T,'CO2')/100*exp(MACC('CO2')*ABL(T,'CO2'))
                            -MACB('CO2')*ENG(T,'CO2')*ABL(T,'CO2')/100
                            -MACB('CO2')/MACC('CO2')*ENG(T,'CO2')/100;

ABCQ2(T)                 .. ABC(T,'CH4') =E=
                            (1/2*MACA('CH4')*(ENG(T,'CH4')+NEG(T,'CH4'))/1000/100*ABL(T,'CH4')*ABL(T,'CH4')
                            +MACB('CH4')/MACC('CH4')*(ENG(T,'CH4')+NEG(T,'CH4'))/1000/100*exp(MACC('CH4')*ABL(T,'CH4'))
                            -MACB('CH4')*(ENG(T,'CH4')+NEG(T,'CH4'))*ABL(T,'CH4')/1000/100
                            -MACB('CH4')/MACC('CH4')*(ENG(T,'CH4')+NEG(T,'CH4'))/1000/100)/IFL;

ABCQ3(T)                 .. ABC(T,'N2O') =E=
                            (1/2*MACA('N2O')*(ENG(T,'N2O')+NEG(T,'N2O'))/1000*(44/28)/100*ABL(T,'N2O')*ABL(T,'N2O')
                            +MACB('N2O')/MACC('N2O')*(ENG(T,'N2O')+NEG(T,'N2O'))/1000*(44/28)/100*exp(MACC('N2O')*ABL(T,'N2O'))
                            -MACB('N2O')*(ENG(T,'N2O')+NEG(T,'N2O'))*ABL(T,'N2O')/1000*(44/28)/100
                            -MACB('N2O')/MACC('N2O')*(ENG(T,'N2O')+NEG(T,'N2O'))/1000*(44/28)/100)/IFL;

* Abatement costs integrated over the time horizon considered (for each agent)
ABCINTQ(FOH)             .. ABCINT(FOH) =E= SUM(T$(Ord(T) LE 201),1/(1+DSC/100)**(Ord(T)-1)*ABC(T,FOH));

* Total abatement costs integrated over the time horizon considered (for all the agents)
ABCTTLQ                  .. ABCTTL =E= SUM(FOH,ABCINT(FOH));

* The current effect of discount rate relative to present
DSCCURQ1(T)$(Ord(T) EQ 1)        .. DSCCUR(T) =E= 1;
DSCCURQ2(T)$(Ord(T) GE 2)        .. DSCCUR(T) =E= DSCCUR(T-1)*(1+DSC/100);

* CO2 equivalent emissions using a metric
ANTCO2EQQ(T)             .. ANTCO2EQ(T)       =E= ANT(T,'CO2')*44/12
                                                 +METRIC('CH4')*ANT(T,'CH4')/1000
                                                 +METRIC('N2O')*ANT(T,'N2O')/1000*44/28
                                                 +METRIC('CF4')*ANT(T,'CF4')/1000/1000
                                                 +METRIC('C2F6')*ANT(T,'C2F6')/1000/1000
                                                 +METRIC('C4F10')*ANT(T,'C4F10')/1000/1000
                                                 +METRIC('HFC23')*ANT(T,'HFC23')/1000/1000
                                                 +METRIC('HFC32')*ANT(T,'HFC32')/1000/1000
                                                 +METRIC('HFC4310mee')*ANT(T,'HFC4310mee')/1000/1000
                                                 +METRIC('HFC125')*ANT(T,'HFC125')/1000/1000
                                                 +METRIC('HFC134a')*ANT(T,'HFC134a')/1000/1000
                                                 +METRIC('HFC143a')*ANT(T,'HFC143a')/1000/1000
                                                 +METRIC('HFC152a')*ANT(T,'HFC152a')/1000/1000
                                                 +METRIC('HFC227ea')*ANT(T,'HFC227ea')/1000/1000
                                                 +METRIC('HFC236fa')*ANT(T,'HFC236fa')/1000/1000
                                                 +METRIC('HFC245ca')*ANT(T,'HFC245ca')/1000/1000
                                                 +METRIC('SF6')*ANT(T,'SF6')/1000/1000;

BSLCO2EQQ(T)             .. BSLCO2EQ(T)       =E= (ENG(T,'CO2')+NEG(T,'CO2'))*44/12
                                                 +METRIC('CH4')*(ENG(T,'CH4')+NEG(T,'CH4'))/1000
                                                 +METRIC('N2O')*(ENG(T,'N2O')+NEG(T,'N2O'))/1000*44/28
                                                 +METRIC('CF4')*((ENGSCN(T,'CF4','A2rGB'))$(Ord(T) LT 101)
                                                 +(ENGSCN('2100','CF4','A2rGB'))$(Ord(T) GE 101))/1000/1000
                                                 +METRIC('C2F6')*((ENGSCN(T,'C2F6','A2rGB'))$(Ord(T) LT 101)
                                                 +(ENGSCN('2100','C2F6','A2rGB'))$(Ord(T) GE 101))/1000/1000
                                                 +METRIC('C4F10')*((ENGSCN(T,'C4F10','A2rGB'))$(Ord(T) LT 101)
                                                 +(ENGSCN('2100','C4F10','A2rGB'))$(Ord(T) GE 101))/1000/1000
                                                 +METRIC('HFC23')*((ENGSCN(T,'HFC23','A2rGB'))$(Ord(T) LT 101)
                                                 +(ENGSCN('2100','HFC23','A2rGB'))$(Ord(T) GE 101))/1000/1000
                                                 +METRIC('HFC32')*((ENGSCN(T,'HFC32','A2rGB'))$(Ord(T) LT 101)
                                                 +(ENGSCN('2100','HFC32','A2rGB'))$(Ord(T) GE 101))/1000/1000
                                                 +METRIC('HFC4310mee')*((ENGSCN(T,'HFC4310mee','A2rGB'))$(Ord(T) LT 101)
                                                 +(ENGSCN('2100','HFC4310mee','A2rGB'))$(Ord(T) GE 101))/1000/1000
                                                 +METRIC('HFC125')*((ENGSCN(T,'HFC125','A2rGB'))$(Ord(T) LT 101)
                                                 +(ENGSCN('2100','HFC125','A2rGB'))$(Ord(T) GE 101))/1000/1000
                                                 +METRIC('HFC134a')*((ENGSCN(T,'HFC134a','A2rGB'))$(Ord(T) LT 101)
                                                 +(ENGSCN('2100','HFC134a','A2rGB'))$(Ord(T) GE 101))/1000/1000
                                                 +METRIC('HFC143a')*((ENGSCN(T,'HFC143a','A2rGB'))$(Ord(T) LT 101)
                                                 +(ENGSCN('2100','HFC143a','A2rGB'))$(Ord(T) GE 101))/1000/1000
                                                 +METRIC('HFC152a')*((ENGSCN(T,'HFC152a','A2rGB'))$(Ord(T) LT 101)
                                                 +(ENGSCN('2100','HFC152a','A2rGB'))$(Ord(T) GE 101))/1000/1000
                                                 +METRIC('HFC227ea')*((ENGSCN(T,'HFC227ea','A2rGB'))$(Ord(T) LT 101)
                                                 +(ENGSCN('2100','HFC227ea','A2rGB'))$(Ord(T) GE 101))/1000/1000
                                                 +METRIC('HFC236fa')*((ENGSCN(T,'HFC236fa','A2rGB'))$(Ord(T) LT 101)
                                                 +(ENGSCN('2100','HFC236fa','A2rGB'))$(Ord(T) GE 101))/1000/1000
                                                 +METRIC('HFC245ca')*((ENGSCN(T,'HFC245ca','A2rGB'))$(Ord(T) LT 101)
                                                 +(ENGSCN('2100','HFC245ca','A2rGB'))$(Ord(T) GE 101))/1000/1000
                                                 +METRIC('SF6')*((ENGSCN(T,'SF6','A2rGB'))$(Ord(T) LT 101)
                                                 +(ENGSCN('2100','SF6','A2rGB'))$(Ord(T) GE 101))/1000/1000;

*First derivative of GHG emissions (annual change)
ANTCO2EQFIRQ(T)$(Ord(T) GT 1) .. ANTCO2EQFIR(T) =E= ANTCO2EQ(T)-ANTCO2EQ(T-1);

*Second derivative of GHG emissions
ANTCO2EQSECQ(T)$(Ord(T) GT 2) .. ANTCO2EQSEC(T) =E= ANTCO2EQFIR(T)-ANTCO2EQFIR(T-1);

*First derivative of CO2 emissions (annual change)
ANTCO2FIRQ(T)$(Ord(T) GT 1) .. ANTCO2FIR(T) =E= ANT(T,'CO2')-ANT(T-1,'CO2');

*Second derivative of CO2 emissions
ANTCO2SECQ(T)$(Ord(T) GT 2) .. ANTCO2SEC(T) =E= ANTCO2FIR(T)-ANTCO2FIR(T-1);

*First derivative of CH4 emissions (annual change)
ANTCH4FIRQ(T)$(Ord(T) GT 1) .. ANTCH4FIR(T) =E= ANT(T,'CH4')-ANT(T-1,'CH4');

*Second derivative of CH4 emissions
ANTCH4SECQ(T)$(Ord(T) GT 2) .. ANTCH4SEC(T) =E= ANTCH4FIR(T)-ANTCH4FIR(T-1);

*First derivative of N2O emissions (annual change)
ANTN2OFIRQ(T)$(Ord(T) GT 1) .. ANTN2OFIR(T) =E= ANT(T,'N2O')-ANT(T-1,'N2O');

*Second derivative of N2O emissions
ANTN2OSECQ(T)$(Ord(T) GT 2) .. ANTN2OSEC(T) =E= ANTN2OFIR(T)-ANTN2OFIR(T-1);


****** Emission targets ================================================================================================

*2030 NDC constraint
*ANTCO2EQ.UP(T)$(Ord(T) EQ 31) = 55.5;

*Net zero GHG emission target yera
ZEY = 61;

*Net zero GHG emission constraint
*ANTCO2EQ.FX(T)$(Ord(T) GE ZEY) = 0;

*Peaking out the GHG emissions collectively and individually by 2030
*ANTCO2EQFIR.UP(T)$(Ord(T) GE 32) = 0;
*ANTCO2FIR.UP(T)$(Ord(T) GE 32) = 0;
*ANTCH4FIR.UP(T)$(Ord(T) GE 32) = 0;
*ANTN2OFIR.UP(T)$(Ord(T) GE 32) = 0;

*Constant GHG emissions reduction collectively and individually from 2040
*ANTCO2EQSEC.FX(T)$((Ord(T) GE 37) and (Ord(T) LE ZEY-5)) = 0;
*ANTCO2SEC.FX(T)$((Ord(T) GE 37) and (Ord(T) LE ZEY-5)) = 0;
*ANTCH4SEC.FX(T)$((Ord(T) GE 37) and (Ord(T) LE ZEY-5)) = 0;
*ANTN2OSEC.FX(T)$((Ord(T) GE 37) and (Ord(T) LE ZEY-5)) = 0;

*Specifying the emission metric
METRIC.FX(FOA) = GWP100(FOA);


****** Variable bounds =================================================================================================

*ABL.UP(T,FOH) = 100;
ABL.LO(T,FOH) = 0;

$ontext
****** Variable scalings ===============================================================================================

MAC.SCALE(T,FOH) = 1000;
ABC.SCALE(T,FOH) = 1000;
ABCINT.SCALE(FOH) = 10000;
ABCTTL.SCALE = 10000;
CF.SCALE = 10000;


****** Equation scalings ===============================================================================================

MACQ1.SCALE(T) = 1000;
MACQ2.SCALE(T) = 1000;
MACQ3.SCALE(T) = 1000;
ABCQ1.SCALE(T) = 10000;
ABCQ2.SCALE(T) = 10000;
ABCQ3.SCALE(T) = 10000;
ABCINTQ.SCALE(FOH) = 10000;
ABCTTLQ.SCALE = 10000;
$offtext
