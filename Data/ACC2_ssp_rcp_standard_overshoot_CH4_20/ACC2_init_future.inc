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
 BETA.FX('2000')              =    0.9000000000;

* Q10 for heterotrophic respiration
 Q10BAS.FX                    =    1.3000000000;

* Temperature sensitivity to Q10
* (not used in the standard inversion.)
 Q10TMP.FX                    =    0.0460000000;

* Temperature sensitivity to NPP
* (not used in the standard inversion.)
 TNPPSENS.FX                  =    0.0160000000;

* Preindustrial ocean CO2 uptake
 CO2OCNUPPRE.FX('2000')       =   -0.4351132724;

* Preindustrial land CO2 uptake
 CO2LNDUPPRE.FX('2000')       =   -0.0048644873;

* Natural methane emission
 NATCH4.FX                    =  319.8976970270;

* Natural nitrous oxide emission
 NATN2O.FX                    =   11.3426207745;

* CH4 lifetime with respect to OH depletion
 TAUCH4OH.FX                  =    8.5348039481;

* N2O lifetime
 TAUN2O.FX                    =  114.3500146401;

* Radiative forcing of a doubling of atmospheric CO2 concentration
* (fixed in the standard inversion)
 RAD2XCO2.FX                  =    3.7000000000;

* Forcing scaling factor
* (fixed in the standard inversion)
 SCLRAD.FX                    =    1.1046724311;

* Climate sensitivity
 SENS.FX('2000')              =    3.8000000000;

* Ocean vertical diffusivity
* (fixed in the standard inversion)
 KAPPA.FX                     =    0.5500000000;

* Propagator for forcing residual AR(1) model
* (not used in the standard inversion)
 RADMISAR1PRO.FX              =    0.5000000000;

* Propagator for temperature residual AR(1) model
* (not used in the standard inversion)
 T2MAR1PRO.FX                 =    0.4500000000;


*** Data and other variables:

* Preindustrial carbon content in the ocean mixed layer
 COC1MIXBAC.FX                =  505.2191379927;

* Preindustrial carbon content in the atmosphere
 COC1ATMBAC.FX                =  587.6357142857;

* Carbon content change in each layer of the ocean carbon cycle model in 2000
 COC.FX('2000','1')           =  211.2695465612;
 COC.FX('2000','2')           =   54.9926821263;
 COC.FX('2000','3')           =   42.2254755536;
 COC.FX('2000','4')           =   13.3578930198;

* Carbon content change in the entire ocean in 1999
 SUMCOC1999.FX                =  123.7287165541;

* Carbon content change in the ocean mixed layer in 2000
 COC1MIX.FX('2000')           =   15.1573547394;

* Carbon content change in each layer of the land carbon cycle model in 2000
 CTR.FX('2000','1')           = -250.4742925698;
 CTR.FX('2000','2')           =  419.2462323407;
 CTR.FX('2000','3')           =  377.6472534873;
 CTR.FX('2000','4')           = 1920.5321062144;

* Carbon content change in the entire land biosphere in 1999
 SUMCTR1999.FX                = 2463.2692450560;

* Carbon content in the land biosphere in 1750
 CTR1750.FX('1')              = -205.2067684571;
 CTR1750.FX('2')              =  343.9893984777;
 CTR1750.FX('3')              =  321.9360000000;
 CTR1750.FX('4')              = 1759.3800000000;

* Carbon content change in the atmosphere in 2000
 COC1ATM.FX('2000')           =  196.1121918219;

* Cumulative anthropogenic CO2 emissions in 1999
 CUMCO21999.FX                =  454.4050786254;

* Atmospheric CO2 concentration in 2000
* (Do not transfer CON('2000','CO2') because of the different ENSO treatment in the past and future modes)
 CONCO2INT.FX('2000')         =  369.4434574245;

* Change in the atmospheric CO2 concentration due to ENSO in 2000
 DCONCO2ENSO2000.FX           =    0.2509629338;

* Atmospheric concentrations of non-CO2 gases in 2000
 CON.FX('2000','CH4')         = 1752.1096778748;
 CON.FX('2000','N2O')         =  317.0441334096;
 CON.FX('2000','CF4')         =   77.2000000000;
 CON.FX('2000','C2F6')        =    3.1700000000;
 CON.FX('2000','C4F10')       =    0.0000000000;
 CON.FX('2000','HFC23')       =   15.5000000000;
 CON.FX('2000','HFC32')       =    0.0000000000;
 CON.FX('2000','HFC4310mee')  =    0.0000000000;
 CON.FX('2000','HFC125')      =    1.4000000000;
 CON.FX('2000','HFC134a')     =   13.8000000000;
 CON.FX('2000','HFC143a')     =    0.0000000000;
 CON.FX('2000','HFC152a')     =    1.7000000000;
 CON.FX('2000','HFC227ea')    =    0.0000000000;
 CON.FX('2000','HFC236fa')    =    0.0000000000;
 CON.FX('2000','HFC245ca')    =    0.0000000000;
 CON.FX('2000','CFC11')       =  262.6000000000;
 CON.FX('2000','CFC12')       =  540.5000000000;
 CON.FX('2000','CFC113')      =   81.7000000000;
 CON.FX('2000','CFC114')      =   16.4000000000;
 CON.FX('2000','CFC115')      =    9.0000000000;
 CON.FX('2000','CCl4')        =   97.9000000000;
 CON.FX('2000','CH3CCl3')     =   49.9000000000;
 CON.FX('2000','Halon1211')   =    3.9700000000;
 CON.FX('2000','Halon1301')   =    2.8200000000;
 CON.FX('2000','Halon2402')   =    0.4100000000;
 CON.FX('2000','HCFC22')      =  141.6000000000;
 CON.FX('2000','HCFC141b')    =   10.8000000000;
 CON.FX('2000','HCFC142b')    =   10.5000000000;
 CON.FX('2000','HCFC123')     =    0.0000000000;
 CON.FX('2000','CH3Cl')       =  536;
 CON.FX('2000','CH3Br')       =  8.1;
 CON.FX('2000','SF6')         =    4.5300000000;
 CON.FX('2000','O3TRP')       =  34.0;

* Concentrations prior to 2000 used for the calculation of EESCl
 CON1999.FX('CFC11')          =  264.7000000000;
 CON1999.FX('CFC12')          =  539.5000000000;
 CON1999.FX('CFC113')         =   82.3000000000;
 CON1999.FX('CFC114')         =   16.3000000000;
 CON1999.FX('CFC115')         =    8.8000000000;
 CON1999.FX('CCl4')           =   99.2000000000;
 CON1999.FX('CH3CCl3')        =   58.6000000000;
 CON1999.FX('Halon1211')      =    3.8300000000;
 CON1999.FX('Halon1301')      =    2.7300000000;
 CON1999.FX('Halon2402')      =    0.4200000000;
 CON1999.FX('HCFC22')         =  136.1000000000;
 CON1999.FX('HCFC141b')       =    9.0000000000;
 CON1999.FX('HCFC142b')       =    9.4000000000;
 CON1999.FX('HCFC123')        =    0.0000000000;
 CON1999.FX('CH3Cl')          =  537;
 CON1999.FX('CH3Br')          =  8.25;

 CON1998.FX('CFC11')          =  266.6000000000;
 CON1998.FX('CFC12')          =  537.4000000000;
 CON1998.FX('CFC113')         =   82.6000000000;
 CON1998.FX('CFC114')         =   16.3000000000;
 CON1998.FX('CFC115')         =    8.6000000000;
 CON1998.FX('CCl4')           =  100.5000000000;
 CON1998.FX('CH3CCl3')        =   69.3000000000;
 CON1998.FX('Halon1211')      =    3.6900000000;
 CON1998.FX('Halon1301')      =    2.6400000000;
 CON1998.FX('Halon2402')      =    0.4300000000;
 CON1998.FX('HCFC22')         =  130.5000000000;
 CON1998.FX('HCFC141b')       =    7.1000000000;
 CON1998.FX('HCFC142b')       =    8.5000000000;
 CON1998.FX('HCFC123')        =    0.0000000000;
 CON1998.FX('CH3Cl')          =  538;
 CON1998.FX('CH3Br')          =  8.4;

 CON1997.FX('CFC11')          =  268.2000000000;
 CON1997.FX('CFC12')          =  534.2000000000;
 CON1997.FX('CFC113')         =   82.7000000000;
 CON1997.FX('CFC114')         =   16.3000000000;
 CON1997.FX('CFC115')         =    8.3000000000;
 CON1997.FX('CCl4')           =  101.6000000000;
 CON1997.FX('CH3CCl3')        =   82.9000000000;
 CON1997.FX('Halon1211')      =    3.5500000000;
 CON1997.FX('Halon1301')      =    2.5500000000;
 CON1997.FX('Halon2402')      =    0.4300000000;
 CON1997.FX('HCFC22')         =  125.0000000000;
 CON1997.FX('HCFC141b')       =    5.5000000000;
 CON1997.FX('HCFC142b')       =    7.4000000000;
 CON1997.FX('HCFC123')        =    0.0000000000;
 CON1997.FX('CH3Cl')          =  539;
 CON1997.FX('CH3Br')          =  8.55;

* Solar forcing in 2000
 RADSOL2000.FX                =    0.1340000000;

* Solar forcing averaged between 1901 and 2000
 RADSOL20AVE.FX               =    0.0331400000;

* Average aerosol forcing between 1901 and 2000
 RADAERAVE.FX                 =   -1.1114834838;

* Average missing forcing between 1901 and 2000
 RADMISAVE.FX                 =    0.0000000000;

* Principal components of temperature change
 PCTRES.FX('2000','1')        =    0.0268998118;
 PCTRES.FX('2000','2')        =    0.5455788152;

* Global-mean temperature change in 1750 (needed for temperature feedback to land CO2 uptake)
 T2M1750.FX                   =    0.0000000000;

* Global-mean temperature change in 1990 (needed for sea level calculation)
 T2M1990.FX                   =    0.6274777317;

* Land surface temperature in 1750 defined in DOECLIM
 TAL1750.FX                   =    0.0000000000;

* Land surface temperature in 2000 defined in DOECLIM
 TAL.FX('2000')               =    1.0126198732;

* Sea surface temperature in 1750 defined in DOECLIM
 TSS1750.FX                   =    0.0000000000;

* Sea surface temperature in 2000 defined in DOECLIM
 TSS.FX('2000')               =    0.6573629177;

* Temperature data transfer for DOECLIM memory calculation
* Note that the index is used in an opposite way of TSS(T).
 TSSMM.FX('1')                =    0.6573629177;
 TSSMM.FX('2')                =    0.6264275166;
 TSSMM.FX('3')                =    0.5947607802;
 TSSMM.FX('4')                =    0.5638831968;
 TSSMM.FX('5')                =    0.5336953791;
 TSSMM.FX('6')                =    0.4987882246;
 TSSMM.FX('7')                =    0.4546955162;
 TSSMM.FX('8')                =    0.4164218449;
 TSSMM.FX('9')                =    0.4445839894;
 TSSMM.FX('10')               =    0.5145667405;
 TSSMM.FX('11')               =    0.5278267087;
 TSSMM.FX('12')               =    0.4952067224;
 TSSMM.FX('13')               =    0.4623000200;
 TSSMM.FX('14')               =    0.4331527899;
 TSSMM.FX('15')               =    0.4047377765;
 TSSMM.FX('16')               =    0.3734350169;
 TSSMM.FX('17')               =    0.3444601645;
 TSSMM.FX('18')               =    0.3533973337;
 TSSMM.FX('19')               =    0.3992846950;
 TSSMM.FX('20')               =    0.4121376896;
 TSSMM.FX('21')               =    0.3871028569;
 TSSMM.FX('22')               =    0.3598554175;
 TSSMM.FX('23')               =    0.3331854216;
 TSSMM.FX('24')               =    0.3080812119;
 TSSMM.FX('25')               =    0.2892926956;
 TSSMM.FX('26')               =    0.2896555333;
 TSSMM.FX('27')               =    0.2911146952;
 TSSMM.FX('28')               =    0.2780902406;
 TSSMM.FX('29')               =    0.2614047717;
 TSSMM.FX('30')               =    0.2410776615;
 TSSMM.FX('31')               =    0.2208926051;
 TSSMM.FX('32')               =    0.2163102988;
 TSSMM.FX('33')               =    0.2215668754;
 TSSMM.FX('34')               =    0.2135335703;
 TSSMM.FX('35')               =    0.1924910156;
 TSSMM.FX('36')               =    0.1731834665;
 TSSMM.FX('37')               =    0.1893917803;
 TSSMM.FX('38')               =    0.2573345865;
 TSSMM.FX('39')               =    0.3020465323;
 TSSMM.FX('40')               =    0.3002198789;
 TSSMM.FX('41')               =    0.2932992386;
 TSSMM.FX('42')               =    0.2853212620;
 TSSMM.FX('43')               =    0.2761326568;
 TSSMM.FX('44')               =    0.2660555936;
 TSSMM.FX('45')               =    0.2594585726;
 TSSMM.FX('46')               =    0.2549163416;
 TSSMM.FX('47')               =    0.2541959765;
 TSSMM.FX('48')               =    0.2588201023;
 TSSMM.FX('49')               =    0.2627554796;
 TSSMM.FX('50')               =    0.2658279491;
 TSSMM.FX('51')               =    0.2662145390;
 TSSMM.FX('52')               =    0.2622324378;
 TSSMM.FX('53')               =    0.2542604049;
 TSSMM.FX('54')               =    0.2459611994;
 TSSMM.FX('55')               =    0.2417894910;
 TSSMM.FX('56')               =    0.2378964863;
 TSSMM.FX('57')               =    0.2347515053;
 TSSMM.FX('58')               =    0.2345828879;
 TSSMM.FX('59')               =    0.2351981898;
 TSSMM.FX('60')               =    0.2339544528;
 TSSMM.FX('61')               =    0.2293541038;
 TSSMM.FX('62')               =    0.2215660737;
 TSSMM.FX('63')               =    0.2121234515;
 TSSMM.FX('64')               =    0.2021406325;
 TSSMM.FX('65')               =    0.1903064800;
 TSSMM.FX('66')               =    0.1781047916;
 TSSMM.FX('67')               =    0.1676149432;
 TSSMM.FX('68')               =    0.1574971031;
 TSSMM.FX('69')               =    0.1480706939;
 TSSMM.FX('70')               =    0.1434021278;
 TSSMM.FX('71')               =    0.1410212751;
 TSSMM.FX('72')               =    0.1493337487;
 TSSMM.FX('73')               =    0.1607218128;
 TSSMM.FX('74')               =    0.1567215966;
 TSSMM.FX('75')               =    0.1483271711;
 TSSMM.FX('76')               =    0.1422386177;
 TSSMM.FX('77')               =    0.1370785862;
 TSSMM.FX('78')               =    0.1323336051;
 TSSMM.FX('79')               =    0.1350880946;
 TSSMM.FX('80')               =    0.1386709256;
 TSSMM.FX('81')               =    0.1341020187;
 TSSMM.FX('82')               =    0.1265037108;
 TSSMM.FX('83')               =    0.1143498170;
 TSSMM.FX('84')               =    0.0993431680;
 TSSMM.FX('85')               =    0.0822997281;
 TSSMM.FX('86')               =    0.0638188904;
 TSSMM.FX('87')               =    0.0462268208;
 TSSMM.FX('88')               =    0.0398275072;
 TSSMM.FX('89')               =    0.0498356212;
 TSSMM.FX('90')               =    0.0519498448;
 TSSMM.FX('91')               =    0.0398066528;
 TSSMM.FX('92')               =    0.0248694011;
 TSSMM.FX('93')               =    0.0070329464;
 TSSMM.FX('94')               =   -0.0087819822;
 TSSMM.FX('95')               =   -0.0327038746;
 TSSMM.FX('96')               =   -0.0696104409;
 TSSMM.FX('97')               =   -0.0878864808;
 TSSMM.FX('98')               =   -0.0131258030;
 TSSMM.FX('99')               =    0.0745532446;
 TSSMM.FX('100')              =    0.0876453781;
 TSSMM.FX('101')              =    0.0844167332;
 TSSMM.FX('102')              =    0.0796336635;
 TSSMM.FX('103')              =    0.0734768582;
 TSSMM.FX('104')              =    0.0657639337;
 TSSMM.FX('105')              =    0.0542805647;
 TSSMM.FX('106')              =    0.0377986055;
 TSSMM.FX('107')              =    0.0164624085;
 TSSMM.FX('108')              =   -0.0092768826;
 TSSMM.FX('109')              =   -0.0347234399;
 TSSMM.FX('110')              =   -0.0435176259;
 TSSMM.FX('111')              =   -0.0206152600;
 TSSMM.FX('112')              =   -0.0122230919;
 TSSMM.FX('113')              =   -0.0385265064;
 TSSMM.FX('114')              =   -0.0744232735;
 TSSMM.FX('115')              =   -0.1200130350;
 TSSMM.FX('116')              =   -0.1466676142;
 TSSMM.FX('117')              =   -0.0620059744;
 TSSMM.FX('118')              =    0.0572394536;
 TSSMM.FX('119')              =    0.0890202490;
 TSSMM.FX('120')              =    0.0864276366;
 TSSMM.FX('121')              =    0.0848693383;
 TSSMM.FX('122')              =    0.0842621150;
 TSSMM.FX('123')              =    0.0835517738;
 TSSMM.FX('124')              =    0.0832755365;
 TSSMM.FX('125')              =    0.0829972483;
 TSSMM.FX('126')              =    0.0817677123;
 TSSMM.FX('127')              =    0.0789805306;
 TSSMM.FX('128')              =    0.0744834608;
 TSSMM.FX('129')              =    0.0689102784;
 TSSMM.FX('130')              =    0.0622321195;
 TSSMM.FX('131')              =    0.0556646022;
 TSSMM.FX('132')              =    0.0504315107;
 TSSMM.FX('133')              =    0.0461098753;
 TSSMM.FX('134')              =    0.0429515069;
 TSSMM.FX('135')              =    0.0403541247;
 TSSMM.FX('136')              =    0.0376121137;
 TSSMM.FX('137')              =    0.0344509523;
 TSSMM.FX('138')              =    0.0302194875;
 TSSMM.FX('139')              =    0.0240599702;
 TSSMM.FX('140')              =    0.0159227592;
 TSSMM.FX('141')              =    0.0067757694;
 TSSMM.FX('142')              =   -0.0027895427;
 TSSMM.FX('143')              =   -0.0127263381;
 TSSMM.FX('144')              =   -0.0166124985;
 TSSMM.FX('145')              =   -0.0176593780;
 TSSMM.FX('146')              =   -0.0223970543;
 TSSMM.FX('147')              =   -0.0178221999;
 TSSMM.FX('148')              =   -0.0137043254;
 TSSMM.FX('149')              =   -0.0203046824;
 TSSMM.FX('150')              =   -0.0284769775;
 TSSMM.FX('151')              =   -0.0394635060;
 TSSMM.FX('152')              =   -0.0532520853;
 TSSMM.FX('153')              =   -0.0681578181;
 TSSMM.FX('154')              =   -0.0715802331;
 TSSMM.FX('155')              =   -0.0714747913;
 TSSMM.FX('156')              =   -0.0813394667;
 TSSMM.FX('157')              =   -0.0931761315;
 TSSMM.FX('158')              =   -0.1057357445;
 TSSMM.FX('159')              =   -0.1205383958;
 TSSMM.FX('160')              =   -0.1403227932;
 TSSMM.FX('161')              =   -0.1673897997;
 TSSMM.FX('162')              =   -0.2037405944;
 TSSMM.FX('163')              =   -0.2533853787;
 TSSMM.FX('164')              =   -0.3098766787;
 TSSMM.FX('165')              =   -0.3304358793;
 TSSMM.FX('166')              =   -0.2642218599;
 TSSMM.FX('167')              =   -0.2263435997;
 TSSMM.FX('168')              =   -0.2717175043;
 TSSMM.FX('169')              =   -0.3053330624;
 TSSMM.FX('170')              =   -0.2628950504;
 TSSMM.FX('171')              =   -0.1875586344;
 TSSMM.FX('172')              =   -0.1758397966;
 TSSMM.FX('173')              =   -0.2034493318;
 TSSMM.FX('174')              =   -0.2370842545;
 TSSMM.FX('175')              =   -0.2784790520;
 TSSMM.FX('176')              =   -0.3284986847;
 TSSMM.FX('177')              =   -0.3678073594;
 TSSMM.FX('178')              =   -0.3454230467;
 TSSMM.FX('179')              =   -0.3287082656;
 TSSMM.FX('180')              =   -0.3753259386;
 TSSMM.FX('181')              =   -0.4427184749;
 TSSMM.FX('182')              =   -0.5321478144;
 TSSMM.FX('183')              =   -0.6545114373;
 TSSMM.FX('184')              =   -0.7791826276;
 TSSMM.FX('185')              =   -0.7153203063;
 TSSMM.FX('186')              =   -0.4319435679;
 TSSMM.FX('187')              =   -0.2933113167;
 TSSMM.FX('188')              =   -0.3370185463;
 TSSMM.FX('189')              =   -0.3775102851;
 TSSMM.FX('190')              =   -0.4378719001;
 TSSMM.FX('191')              =   -0.4207524695;
 TSSMM.FX('192')              =   -0.1763972239;
 TSSMM.FX('193')              =    0.0130321312;
 TSSMM.FX('194')              =    0.0139876617;
 TSSMM.FX('195')              =    0.0140649565;
 TSSMM.FX('196')              =    0.0131463541;
 TSSMM.FX('197')              =    0.0126384396;
 TSSMM.FX('198')              =    0.0120938964;
 TSSMM.FX('199')              =    0.0098311431;
 TSSMM.FX('200')              =    0.0071593704;
 TSSMM.FX('201')              =    0.0056708839;
 TSSMM.FX('202')              =    0.0043954525;
 TSSMM.FX('203')              =    0.0030296556;
 TSSMM.FX('204')              =    0.0011401091;
 TSSMM.FX('205')              =   -0.0013164001;
 TSSMM.FX('206')              =   -0.0049749310;
 TSSMM.FX('207')              =   -0.0100392923;
 TSSMM.FX('208')              =   -0.0171852915;
 TSSMM.FX('209')              =   -0.0271487244;
 TSSMM.FX('210')              =   -0.0420316655;
 TSSMM.FX('211')              =   -0.0585199729;
 TSSMM.FX('212')              =   -0.0629342999;
 TSSMM.FX('213')              =   -0.0405081243;
 TSSMM.FX('214')              =   -0.0268169145;
 TSSMM.FX('215')              =   -0.0399011683;
 TSSMM.FX('216')              =   -0.0538192452;
 TSSMM.FX('217')              =   -0.0481635630;
 TSSMM.FX('218')              =   -0.0152558623;
 TSSMM.FX('219')              =    0.0031445124;
 TSSMM.FX('220')              =    0.0007366920;
 TSSMM.FX('221')              =   -0.0025933644;
 TSSMM.FX('222')              =   -0.0091782599;
 TSSMM.FX('223')              =   -0.0189307129;
 TSSMM.FX('224')              =   -0.0288855390;
 TSSMM.FX('225')              =   -0.0212860915;
 TSSMM.FX('226')              =   -0.0085100201;
 TSSMM.FX('227')              =   -0.0051867678;
 TSSMM.FX('228')              =   -0.0026817233;
 TSSMM.FX('229')              =   -0.0027007792;
 TSSMM.FX('230')              =   -0.0054265282;
 TSSMM.FX('231')              =   -0.0110996625;
 TSSMM.FX('232')              =   -0.0182134399;
 TSSMM.FX('233')              =   -0.0253318948;
 TSSMM.FX('234')              =   -0.0294529919;
 TSSMM.FX('235')              =   -0.0242878171;
 TSSMM.FX('236')              =   -0.0179728601;
 TSSMM.FX('237')              =   -0.0162490181;
 TSSMM.FX('238')              =   -0.0150097930;
 TSSMM.FX('239')              =   -0.0153029917;
 TSSMM.FX('240')              =   -0.0152250570;
 TSSMM.FX('241')              =   -0.0160146285;
 TSSMM.FX('242')              =   -0.0179281987;
 TSSMM.FX('243')              =   -0.0187047137;
 TSSMM.FX('244')              =   -0.0187907529;
 TSSMM.FX('245')              =   -0.0174677740;
 TSSMM.FX('246')              =   -0.0146296750;
 TSSMM.FX('247')              =   -0.0114549928;
 TSSMM.FX('248')              =   -0.0086986396;
 TSSMM.FX('249')              =   -0.0054776087;
 TSSMM.FX('250')              =   -0.0019092852;
 TSSMM.FX('251')              =    0.0000000000;

* Principal components of sea level rise
 PCSRES.FX('2000','1')        =    0.0198425209;
 PCSRES.FX('2000','2')        =    0.0146507239;

* Sea level rise in 1990 due to thermal expansion
 PCS1990.FX                   =    0.0254665409;

* Sea level rise in 2000 due to loss of mass of glaciers and small ice caps without glacier area contraction
 SLRGLWT.FX('2000')           =    0.0330762255;

* Sea level rise in 2000 due to loss of mass of glaciers and small ice caps with glacier area contraction
 SLRGLWC.FX('2000')           =    0.0296186419;

* Sea level rise in 1990 due to loss of mass of glaciers and small ice caps with glacier area contraction
 SLRGLWC1990.FX               =    0.0243687150;

* Greenland ice sheet contribution to sea level rise in 2000
 SLRGR.FX('2000')             =    0.0042331257;

* Antarctic ice sheet contribution to sea level rise in 2000
 SLRAN.FX('2000')             =   -0.0118425538;

* Sea level rise in 1990 due to four major factors relative to 1750
 SLR1990.FX                   =    0.0458352559;

