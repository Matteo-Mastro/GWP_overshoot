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
 BETA.FX('2000')              =    0.6000000000;

* Q10 for heterotrophic respiration
 Q10BAS.FX                    =    1.2000000000;

* Temperature sensitivity to Q10
* (not used in the standard inversion.)
 Q10TMP.FX                    =    0.0460000000;

* Temperature sensitivity to NPP
* (not used in the standard inversion.)
 TNPPSENS.FX                  =    0.0160000000;

* Preindustrial ocean CO2 uptake
 CO2OCNUPPRE.FX('2000')       =   -0.3635692617;

* Preindustrial land CO2 uptake
 CO2LNDUPPRE.FX('2000')       =    0.1069230783;

* Natural methane emission
 NATCH4.FX                    =  319.9005825195;

* Natural nitrous oxide emission
 NATN2O.FX                    =   11.3519874508;

* CH4 lifetime with respect to OH depletion
 TAUCH4OH.FX                  =    8.5347020252;

* N2O lifetime
 TAUN2O.FX                    =  114.2556629170;

* Radiative forcing of a doubling of atmospheric CO2 concentration
* (fixed in the standard inversion)
 RAD2XCO2.FX                  =    3.7000000000;

* Forcing scaling factor
* (fixed in the standard inversion)
 SCLRAD.FX                    =    1.1308428367;

* Climate sensitivity
 SENS.FX('2000')              =    4.0000000000;

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
 COC.FX('2000','1')           =  211.7638396497;
 COC.FX('2000','2')           =   55.0553925326;
 COC.FX('2000','3')           =   42.2550916019;
 COC.FX('2000','4')           =   13.3581878918;

* Carbon content change in the entire ocean in 1999
 SUMCOC1999.FX                =  123.8366518135;

* Carbon content change in the ocean mixed layer in 2000
 COC1MIX.FX('2000')           =   15.1873155821;

* Carbon content change in each layer of the land carbon cycle model in 2000
 CTR.FX('2000','1')           = -239.4482085772;
 CTR.FX('2000','2')           =  400.7450848193;
 CTR.FX('2000','3')           =  362.9814451066;
 CTR.FX('2000','4')           = 1874.6583388733;

* Carbon content change in the entire land biosphere in 1999
 SUMCTR1999.FX                = 2395.9377456254;

* Carbon content in the land biosphere in 1750
 CTR1750.FX('1')              = -205.2067684571;
 CTR1750.FX('2')              =  343.9893984777;
 CTR1750.FX('3')              =  321.9360000000;
 CTR1750.FX('4')              = 1759.3800000000;

* Carbon content change in the atmosphere in 2000
 COC1ATM.FX('2000')           =  196.5765240676;

* Cumulative anthropogenic CO2 emissions in 1999
 CUMCO21999.FX                =  432.9380707274;

* Atmospheric CO2 concentration in 2000
* (Do not transfer CON('2000','CO2') because of the different ENSO treatment in the past and future modes)
 CONCO2INT.FX('2000')         =  369.6623345773;

* Change in the atmospheric CO2 concentration due to ENSO in 2000
 DCONCO2ENSO2000.FX           =    0.0663240942;

* Atmospheric concentrations of non-CO2 gases in 2000
 CON.FX('2000','CH4')         = 1752.0985264728;
 CON.FX('2000','N2O')         =  317.0390785562;
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
 RADAERAVE.FX                 =   -1.1378152476;

* Average missing forcing between 1901 and 2000
 RADMISAVE.FX                 =    0.0000000000;

* Principal components of temperature change
 PCTRES.FX('2000','1')        =    0.0271034402;
 PCTRES.FX('2000','2')        =    0.5377453261;

* Global-mean temperature change in 1750 (needed for temperature feedback to land CO2 uptake)
 T2M1750.FX                   =    0.0000000000;

* Global-mean temperature change in 1990 (needed for sea level calculation)
 T2M1990.FX                   =    0.6247183250;

* Land surface temperature in 1750 defined in DOECLIM
 TAL1750.FX                   =    0.0000000000;

* Land surface temperature in 2000 defined in DOECLIM
 TAL.FX('2000')               =    1.0134213810;

* Sea surface temperature in 1750 defined in DOECLIM
 TSS1750.FX                   =    0.0000000000;

* Sea surface temperature in 2000 defined in DOECLIM
 TSS.FX('2000')               =    0.6578122657;

* Temperature data transfer for DOECLIM memory calculation
* Note that the index is used in an opposite way of TSS(T).
 TSSMM.FX('1')                =    0.6578122657;
 TSSMM.FX('2')                =    0.6261395345;
 TSSMM.FX('3')                =    0.5938451453;
 TSSMM.FX('4')                =    0.5624223960;
 TSSMM.FX('5')                =    0.5317518522;
 TSSMM.FX('6')                =    0.4964278348;
 TSSMM.FX('7')                =    0.4520479648;
 TSSMM.FX('8')                =    0.4137875099;
 TSSMM.FX('9')                =    0.4423871477;
 TSSMM.FX('10')               =    0.5125980346;
 TSSMM.FX('11')               =    0.5255524921;
 TSSMM.FX('12')               =    0.4924654184;
 TSSMM.FX('13')               =    0.4591426070;
 TSSMM.FX('14')               =    0.4296578097;
 TSSMM.FX('15')               =    0.4009812228;
 TSSMM.FX('16')               =    0.3695121773;
 TSSMM.FX('17')               =    0.3405449070;
 TSSMM.FX('18')               =    0.3497451080;
 TSSMM.FX('19')               =    0.3958126110;
 TSSMM.FX('20')               =    0.4084627357;
 TSSMM.FX('21')               =    0.3831287880;
 TSSMM.FX('22')               =    0.3556264709;
 TSSMM.FX('23')               =    0.3287740183;
 TSSMM.FX('24')               =    0.3035658547;
 TSSMM.FX('25')               =    0.2847765466;
 TSSMM.FX('26')               =    0.2852284495;
 TSSMM.FX('27')               =    0.2867043065;
 TSSMM.FX('28')               =    0.2736627755;
 TSSMM.FX('29')               =    0.2569811676;
 TSSMM.FX('30')               =    0.2367108702;
 TSSMM.FX('31')               =    0.2166985855;
 TSSMM.FX('32')               =    0.2124238553;
 TSSMM.FX('33')               =    0.2179804084;
 TSSMM.FX('34')               =    0.2101821974;
 TSSMM.FX('35')               =    0.1894153662;
 TSSMM.FX('36')               =    0.1705513026;
 TSSMM.FX('37')               =    0.1874917973;
 TSSMM.FX('38')               =    0.2562122190;
 TSSMM.FX('39')               =    0.3012694655;
 TSSMM.FX('40')               =    0.2995830844;
 TSSMM.FX('41')               =    0.2927669170;
 TSSMM.FX('42')               =    0.2848494432;
 TSSMM.FX('43')               =    0.2756741014;
 TSSMM.FX('44')               =    0.2656454085;
 TSSMM.FX('45')               =    0.2591351133;
 TSSMM.FX('46')               =    0.2547018299;
 TSSMM.FX('47')               =    0.2541266337;
 TSSMM.FX('48')               =    0.2589004805;
 TSSMM.FX('49')               =    0.2629517560;
 TSSMM.FX('50')               =    0.2661029179;
 TSSMM.FX('51')               =    0.2665019626;
 TSSMM.FX('52')               =    0.2624955974;
 TSSMM.FX('53')               =    0.2544825430;
 TSSMM.FX('54')               =    0.2461653587;
 TSSMM.FX('55')               =    0.2419884838;
 TSSMM.FX('56')               =    0.2380864322;
 TSSMM.FX('57')               =    0.2349417234;
 TSSMM.FX('58')               =    0.2347754966;
 TSSMM.FX('59')               =    0.2353678379;
 TSSMM.FX('60')               =    0.2340557131;
 TSSMM.FX('61')               =    0.2293301970;
 TSSMM.FX('62')               =    0.2213727413;
 TSSMM.FX('63')               =    0.2117320258;
 TSSMM.FX('64')               =    0.2015273016;
 TSSMM.FX('65')               =    0.1894415312;
 TSSMM.FX('66')               =    0.1769804886;
 TSSMM.FX('67')               =    0.1662349201;
 TSSMM.FX('68')               =    0.1558851479;
 TSSMM.FX('69')               =    0.1462728570;
 TSSMM.FX('70')               =    0.1414866413;
 TSSMM.FX('71')               =    0.1390318496;
 TSSMM.FX('72')               =    0.1473530316;
 TSSMM.FX('73')               =    0.1586863809;
 TSSMM.FX('74')               =    0.1545644578;
 TSSMM.FX('75')               =    0.1460521900;
 TSSMM.FX('76')               =    0.1398605055;
 TSSMM.FX('77')               =    0.1346103711;
 TSSMM.FX('78')               =    0.1297912723;
 TSSMM.FX('79')               =    0.1325056138;
 TSSMM.FX('80')               =    0.1359900275;
 TSSMM.FX('81')               =    0.1312724512;
 TSSMM.FX('82')               =    0.1234895761;
 TSSMM.FX('83')               =    0.1111146640;
 TSSMM.FX('84')               =    0.0958777964;
 TSSMM.FX('85')               =    0.0785983832;
 TSSMM.FX('86')               =    0.0599118454;
 TSSMM.FX('87')               =    0.0421748139;
 TSSMM.FX('88')               =    0.0357563468;
 TSSMM.FX('89')               =    0.0457975351;
 TSSMM.FX('90')               =    0.0478559880;
 TSSMM.FX('91')               =    0.0356220813;
 TSSMM.FX('92')               =    0.0205793515;
 TSSMM.FX('93')               =    0.0026266346;
 TSSMM.FX('94')               =   -0.0132905305;
 TSSMM.FX('95')               =   -0.0373432555;
 TSSMM.FX('96')               =   -0.0742600966;
 TSSMM.FX('97')               =   -0.0921722452;
 TSSMM.FX('98')               =   -0.0165528579;
 TSSMM.FX('99')               =    0.0715773092;
 TSSMM.FX('100')              =    0.0848088186;
 TSSMM.FX('101')              =    0.0816373483;
 TSSMM.FX('102')              =    0.0768712428;
 TSSMM.FX('103')              =    0.0706971556;
 TSSMM.FX('104')              =    0.0629231443;
 TSSMM.FX('105')              =    0.0513216908;
 TSSMM.FX('106')              =    0.0346684210;
 TSSMM.FX('107')              =    0.0131212703;
 TSSMM.FX('108')              =   -0.0128346272;
 TSSMM.FX('109')              =   -0.0384365836;
 TSSMM.FX('110')              =   -0.0472495129;
 TSSMM.FX('111')              =   -0.0243107493;
 TSSMM.FX('112')              =   -0.0161636559;
 TSSMM.FX('113')              =   -0.0428378735;
 TSSMM.FX('114')              =   -0.0791160201;
 TSSMM.FX('115')              =   -0.1249431355;
 TSSMM.FX('116')              =   -0.1513588185;
 TSSMM.FX('117')              =   -0.0657690722;
 TSSMM.FX('118')              =    0.0540733422;
 TSSMM.FX('119')              =    0.0860125557;
 TSSMM.FX('120')              =    0.0834619316;
 TSSMM.FX('121')              =    0.0819268414;
 TSSMM.FX('122')              =    0.0813153275;
 TSSMM.FX('123')              =    0.0805732448;
 TSSMM.FX('124')              =    0.0802515086;
 TSSMM.FX('125')              =    0.0799174049;
 TSSMM.FX('126')              =    0.0786200570;
 TSSMM.FX('127')              =    0.0757495181;
 TSSMM.FX('128')              =    0.0711678182;
 TSSMM.FX('129')              =    0.0655264241;
 TSSMM.FX('130')              =    0.0587844971;
 TSSMM.FX('131')              =    0.0521566804;
 TSSMM.FX('132')              =    0.0468626465;
 TSSMM.FX('133')              =    0.0424793652;
 TSSMM.FX('134')              =    0.0392573578;
 TSSMM.FX('135')              =    0.0365877696;
 TSSMM.FX('136')              =    0.0337608592;
 TSSMM.FX('137')              =    0.0304968020;
 TSSMM.FX('138')              =    0.0261321740;
 TSSMM.FX('139')              =    0.0198010611;
 TSSMM.FX('140')              =    0.0114664410;
 TSSMM.FX('141')              =    0.0021199586;
 TSSMM.FX('142')              =   -0.0076376221;
 TSSMM.FX('143')              =   -0.0177460138;
 TSSMM.FX('144')              =   -0.0217484803;
 TSSMM.FX('145')              =   -0.0229072725;
 TSSMM.FX('146')              =   -0.0277398031;
 TSSMM.FX('147')              =   -0.0232029767;
 TSSMM.FX('148')              =   -0.0191874739;
 TSSMM.FX('149')              =   -0.0259166041;
 TSSMM.FX('150')              =   -0.0342222255;
 TSSMM.FX('151')              =   -0.0453555334;
 TSSMM.FX('152')              =   -0.0592899542;
 TSSMM.FX('153')              =   -0.0743116081;
 TSSMM.FX('154')              =   -0.0777668514;
 TSSMM.FX('155')              =   -0.0777362294;
 TSSMM.FX('156')              =   -0.0877102011;
 TSSMM.FX('157')              =   -0.0996698984;
 TSSMM.FX('158')              =   -0.1123803857;
 TSSMM.FX('159')              =   -0.1273839077;
 TSSMM.FX('160')              =   -0.1474454123;
 TSSMM.FX('161')              =   -0.1748753884;
 TSSMM.FX('162')              =   -0.2116674193;
 TSSMM.FX('163')              =   -0.2617752480;
 TSSMM.FX('164')              =   -0.3185697654;
 TSSMM.FX('165')              =   -0.3390172978;
 TSSMM.FX('166')              =   -0.2723840620;
 TSSMM.FX('167')              =   -0.2346054640;
 TSSMM.FX('168')              =   -0.2801558293;
 TSSMM.FX('169')              =   -0.3136241601;
 TSSMM.FX('170')              =   -0.2705726600;
 TSSMM.FX('171')              =   -0.1948277626;
 TSSMM.FX('172')              =   -0.1831129931;
 TSSMM.FX('173')              =   -0.2108772338;
 TSSMM.FX('174')              =   -0.2447359908;
 TSSMM.FX('175')              =   -0.2863986994;
 TSSMM.FX('176')              =   -0.3366583913;
 TSSMM.FX('177')              =   -0.3760132986;
 TSSMM.FX('178')              =   -0.3534579944;
 TSSMM.FX('179')              =   -0.3370315790;
 TSSMM.FX('180')              =   -0.3842850634;
 TSSMM.FX('181')              =   -0.4525276453;
 TSSMM.FX('182')              =   -0.5429515816;
 TSSMM.FX('183')              =   -0.6662295086;
 TSSMM.FX('184')              =   -0.7910769848;
 TSSMM.FX('185')              =   -0.7258100835;
 TSSMM.FX('186')              =   -0.4407489435;
 TSSMM.FX('187')              =   -0.3019665286;
 TSSMM.FX('188')              =   -0.3458926180;
 TSSMM.FX('189')              =   -0.3865047987;
 TSSMM.FX('190')              =   -0.4466312115;
 TSSMM.FX('191')              =   -0.4281566457;
 TSSMM.FX('192')              =   -0.1815697424;
 TSSMM.FX('193')              =    0.0087578746;
 TSSMM.FX('194')              =    0.0101089645;
 TSSMM.FX('195')              =    0.0105105832;
 TSSMM.FX('196')              =    0.0098624433;
 TSSMM.FX('197')              =    0.0095825251;
 TSSMM.FX('198')              =    0.0092234243;
 TSSMM.FX('199')              =    0.0071082012;
 TSSMM.FX('200')              =    0.0045637495;
 TSSMM.FX('201')              =    0.0031848182;
 TSSMM.FX('202')              =    0.0019986259;
 TSSMM.FX('203')              =    0.0007022225;
 TSSMM.FX('204')              =   -0.0011410758;
 TSSMM.FX('205')              =   -0.0035747168;
 TSSMM.FX('206')              =   -0.0072422498;
 TSSMM.FX('207')              =   -0.0123456371;
 TSSMM.FX('208')              =   -0.0195689138;
 TSSMM.FX('209')              =   -0.0296390066;
 TSSMM.FX('210')              =   -0.0446464102;
 TSSMM.FX('211')              =   -0.0612020697;
 TSSMM.FX('212')              =   -0.0655590446;
 TSSMM.FX('213')              =   -0.0429765520;
 TSSMM.FX('214')              =   -0.0293048077;
 TSSMM.FX('215')              =   -0.0424514439;
 TSSMM.FX('216')              =   -0.0563574180;
 TSSMM.FX('217')              =   -0.0505080006;
 TSSMM.FX('218')              =   -0.0173558523;
 TSSMM.FX('219')              =    0.0011108703;
 TSSMM.FX('220')              =   -0.0012885343;
 TSSMM.FX('221')              =   -0.0046417701;
 TSSMM.FX('222')              =   -0.0112676541;
 TSSMM.FX('223')              =   -0.0210608152;
 TSSMM.FX('224')              =   -0.0309965082;
 TSSMM.FX('225')              =   -0.0232639956;
 TSSMM.FX('226')              =   -0.0104016693;
 TSSMM.FX('227')              =   -0.0070183205;
 TSSMM.FX('228')              =   -0.0044723254;
 TSSMM.FX('229')              =   -0.0044750853;
 TSSMM.FX('230')              =   -0.0072080629;
 TSSMM.FX('231')              =   -0.0129054829;
 TSSMM.FX('232')              =   -0.0200397999;
 TSSMM.FX('233')              =   -0.0271581551;
 TSSMM.FX('234')              =   -0.0312288037;
 TSSMM.FX('235')              =   -0.0259596064;
 TSSMM.FX('236')              =   -0.0195628914;
 TSSMM.FX('237')              =   -0.0177648119;
 TSSMM.FX('238')              =   -0.0164560901;
 TSSMM.FX('239')              =   -0.0166828452;
 TSSMM.FX('240')              =   -0.0165327066;
 TSSMM.FX('241')              =   -0.0172566742;
 TSSMM.FX('242')              =   -0.0190940214;
 TSSMM.FX('243')              =   -0.0197797489;
 TSSMM.FX('244')              =   -0.0197582371;
 TSSMM.FX('245')              =   -0.0183071510;
 TSSMM.FX('246')              =   -0.0153259636;
 TSSMM.FX('247')              =   -0.0120007476;
 TSSMM.FX('248')              =   -0.0090897086;
 TSSMM.FX('249')              =   -0.0057111020;
 TSSMM.FX('250')              =   -0.0020003137;
 TSSMM.FX('251')              =    0.0000000000;

* Principal components of sea level rise
 PCSRES.FX('2000','1')        =    0.0196854921;
 PCSRES.FX('2000','2')        =    0.0148876630;

* Sea level rise in 1990 due to thermal expansion
 PCS1990.FX                   =    0.0253264155;

* Sea level rise in 2000 due to loss of mass of glaciers and small ice caps without glacier area contraction
 SLRGLWT.FX('2000')           =    0.0330797306;

* Sea level rise in 2000 due to loss of mass of glaciers and small ice caps with glacier area contraction
 SLRGLWC.FX('2000')           =    0.0296216454;

* Sea level rise in 1990 due to loss of mass of glaciers and small ice caps with glacier area contraction
 SLRGLWC1990.FX               =    0.0243687150;

* Greenland ice sheet contribution to sea level rise in 2000
 SLRGR.FX('2000')             =    0.0042332852;

* Antarctic ice sheet contribution to sea level rise in 2000
 SLRAN.FX('2000')             =   -0.0118451840;

* Sea level rise in 1990 due to four major factors relative to 1750
 SLR1990.FX                   =    0.0456951305;

