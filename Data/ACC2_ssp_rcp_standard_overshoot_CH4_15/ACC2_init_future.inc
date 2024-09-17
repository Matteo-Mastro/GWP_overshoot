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
 Q10BAS.FX                    =    2.0000000000;

* Temperature sensitivity to Q10
* (not used in the standard inversion.)
 Q10TMP.FX                    =    0.0460000000;

* Temperature sensitivity to NPP
* (not used in the standard inversion.)
 TNPPSENS.FX                  =    0.0160000000;

* Preindustrial ocean CO2 uptake
 CO2OCNUPPRE.FX('2000')       =   -0.4213980682;

* Preindustrial land CO2 uptake
 CO2LNDUPPRE.FX('2000')       =    0.0165655196;

* Natural methane emission
 NATCH4.FX                    =  319.8716210678;

* Natural nitrous oxide emission
 NATN2O.FX                    =   11.3267812707;

* CH4 lifetime with respect to OH depletion
 TAUCH4OH.FX                  =    8.5357566592;

* N2O lifetime
 TAUN2O.FX                    =  114.5099230420;

* Radiative forcing of a doubling of atmospheric CO2 concentration
* (fixed in the standard inversion)
 RAD2XCO2.FX                  =    3.7000000000;

* Forcing scaling factor
* (fixed in the standard inversion)
 SCLRAD.FX                    =    0.8750493350;

* Climate sensitivity
 SENS.FX('2000')              =    2.4000000000;

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
 COC.FX('2000','1')           =  211.9881324135;
 COC.FX('2000','2')           =   55.0382129265;
 COC.FX('2000','3')           =   42.2159211171;
 COC.FX('2000','4')           =   13.2898432096;

* Carbon content change in the entire ocean in 1999
 SUMCOC1999.FX                =  123.7090891399;

* Carbon content change in the ocean mixed layer in 2000
 COC1MIX.FX('2000')           =   15.2009036538;

* Carbon content change in each layer of the land carbon cycle model in 2000
 CTR.FX('2000','1')           = -241.8033953695;
 CTR.FX('2000','2')           =  405.0869801428;
 CTR.FX('2000','3')           =  368.6536329425;
 CTR.FX('2000','4')           = 1900.0809246847;

* Carbon content change in the entire land biosphere in 1999
 SUMCTR1999.FX                = 2429.3638254947;

* Carbon content in the land biosphere in 1750
 CTR1750.FX('1')              = -205.2067684571;
 CTR1750.FX('2')              =  343.9893984777;
 CTR1750.FX('3')              =  321.9360000000;
 CTR1750.FX('4')              = 1759.3800000000;

* Carbon content change in the atmosphere in 2000
 COC1ATM.FX('2000')           =  196.7872287597;

* Cumulative anthropogenic CO2 emissions in 1999
 CUMCO21999.FX                =  429.1537966830;

* Atmospheric CO2 concentration in 2000
* (Do not transfer CON('2000','CO2') because of the different ENSO treatment in the past and future modes)
 CONCO2INT.FX('2000')         =  369.7616566544;

* Change in the atmospheric CO2 concentration due to ENSO in 2000
 DCONCO2ENSO2000.FX           =    0.1628577668;

* Atmospheric concentrations of non-CO2 gases in 2000
 CON.FX('2000','CH4')         = 1752.2306938441;
 CON.FX('2000','N2O')         =  317.0832658512;
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
 RADAERAVE.FX                 =   -0.8804446059;

* Average missing forcing between 1901 and 2000
 RADMISAVE.FX                 =    0.0000000000;

* Principal components of temperature change
 PCTRES.FX('2000','1')        =    0.0210128329;
 PCTRES.FX('2000','2')        =    0.5719577976;

* Global-mean temperature change in 1750 (needed for temperature feedback to land CO2 uptake)
 T2M1750.FX                   =    0.0000000000;

* Global-mean temperature change in 1990 (needed for sea level calculation)
 T2M1990.FX                   =    0.6024632130;

* Land surface temperature in 1750 defined in DOECLIM
 TAL1750.FX                   =    0.0000000000;

* Land surface temperature in 2000 defined in DOECLIM
 TAL.FX('2000')               =    0.9393982581;

* Sea surface temperature in 1750 defined in DOECLIM
 TSS1750.FX                   =    0.0000000000;

* Sea surface temperature in 2000 defined in DOECLIM
 TSS.FX('2000')               =    0.6080870939;

* Temperature data transfer for DOECLIM memory calculation
* Note that the index is used in an opposite way of TSS(T).
 TSSMM.FX('1')                =    0.6080870939;
 TSSMM.FX('2')                =    0.5833537967;
 TSSMM.FX('3')                =    0.5576526356;
 TSSMM.FX('4')                =    0.5322340712;
 TSSMM.FX('5')                =    0.5070547218;
 TSSMM.FX('6')                =    0.4765922143;
 TSSMM.FX('7')                =    0.4355518830;
 TSSMM.FX('8')                =    0.3969651635;
 TSSMM.FX('9')                =    0.4202967076;
 TSSMM.FX('10')               =    0.4881823883;
 TSSMM.FX('11')               =    0.5053438852;
 TSSMM.FX('12')               =    0.4781822476;
 TSSMM.FX('13')               =    0.4500788670;
 TSSMM.FX('14')               =    0.4249436975;
 TSSMM.FX('15')               =    0.3998096413;
 TSSMM.FX('16')               =    0.3707660141;
 TSSMM.FX('17')               =    0.3420188978;
 TSSMM.FX('18')               =    0.3481558024;
 TSSMM.FX('19')               =    0.3921480510;
 TSSMM.FX('20')               =    0.4075727476;
 TSSMM.FX('21')               =    0.3863457606;
 TSSMM.FX('22')               =    0.3625321157;
 TSSMM.FX('23')               =    0.3385960545;
 TSSMM.FX('24')               =    0.3153618591;
 TSSMM.FX('25')               =    0.2971894112;
 TSSMM.FX('26')               =    0.2971345137;
 TSSMM.FX('27')               =    0.2992329958;
 TSSMM.FX('28')               =    0.2874338646;
 TSSMM.FX('29')               =    0.2718026547;
 TSSMM.FX('30')               =    0.2521077473;
 TSSMM.FX('31')               =    0.2313983344;
 TSSMM.FX('32')               =    0.2247795106;
 TSSMM.FX('33')               =    0.2280979100;
 TSSMM.FX('34')               =    0.2188312382;
 TSSMM.FX('35')               =    0.1958541400;
 TSSMM.FX('36')               =    0.1724352572;
 TSSMM.FX('37')               =    0.1812295489;
 TSSMM.FX('38')               =    0.2415534637;
 TSSMM.FX('39')               =    0.2839088322;
 TSSMM.FX('40')               =    0.2820118451;
 TSSMM.FX('41')               =    0.2754125903;
 TSSMM.FX('42')               =    0.2677612994;
 TSSMM.FX('43')               =    0.2588584519;
 TSSMM.FX('44')               =    0.2487468412;
 TSSMM.FX('45')               =    0.2416653318;
 TSSMM.FX('46')               =    0.2363817804;
 TSSMM.FX('47')               =    0.2345394089;
 TSSMM.FX('48')               =    0.2380606431;
 TSSMM.FX('49')               =    0.2413334404;
 TSSMM.FX('50')               =    0.2441392774;
 TSSMM.FX('51')               =    0.2449497995;
 TSSMM.FX('52')               =    0.2417981435;
 TSSMM.FX('53')               =    0.2348374306;
 TSSMM.FX('54')               =    0.2272338348;
 TSSMM.FX('55')               =    0.2235385803;
 TSSMM.FX('56')               =    0.2201258167;
 TSSMM.FX('57')               =    0.2173699521;
 TSSMM.FX('58')               =    0.2176219021;
 TSSMM.FX('59')               =    0.2189786900;
 TSSMM.FX('60')               =    0.2189379632;
 TSSMM.FX('61')               =    0.2160270491;
 TSSMM.FX('62')               =    0.2102070929;
 TSSMM.FX('63')               =    0.2027997033;
 TSSMM.FX('64')               =    0.1948345264;
 TSSMM.FX('65')               =    0.1849998819;
 TSSMM.FX('66')               =    0.1744929388;
 TSSMM.FX('67')               =    0.1654046763;
 TSSMM.FX('68')               =    0.1563990588;
 TSSMM.FX('69')               =    0.1476752077;
 TSSMM.FX('70')               =    0.1431042893;
 TSSMM.FX('71')               =    0.1405640211;
 TSSMM.FX('72')               =    0.1481079503;
 TSSMM.FX('73')               =    0.1598017141;
 TSSMM.FX('74')               =    0.1570694319;
 TSSMM.FX('75')               =    0.1498809322;
 TSSMM.FX('76')               =    0.1448110902;
 TSSMM.FX('77')               =    0.1405537134;
 TSSMM.FX('78')               =    0.1365775188;
 TSSMM.FX('79')               =    0.1398087591;
 TSSMM.FX('80')               =    0.1446655289;
 TSSMM.FX('81')               =    0.1420164606;
 TSSMM.FX('82')               =    0.1367087014;
 TSSMM.FX('83')               =    0.1271681802;
 TSSMM.FX('84')               =    0.1147315844;
 TSSMM.FX('85')               =    0.1001465821;
 TSSMM.FX('86')               =    0.0836295259;
 TSSMM.FX('87')               =    0.0672726447;
 TSSMM.FX('88')               =    0.0608717443;
 TSSMM.FX('89')               =    0.0707139257;
 TSSMM.FX('90')               =    0.0740973312;
 TSSMM.FX('91')               =    0.0637154813;
 TSSMM.FX('92')               =    0.0504630029;
 TSSMM.FX('93')               =    0.0340456354;
 TSSMM.FX('94')               =    0.0190162597;
 TSSMM.FX('95')               =   -0.0043627275;
 TSSMM.FX('96')               =   -0.0426825930;
 TSSMM.FX('97')               =   -0.0667669313;
 TSSMM.FX('98')               =   -0.0025305653;
 TSSMM.FX('99')               =    0.0808505648;
 TSSMM.FX('100')              =    0.0941992667;
 TSSMM.FX('101')              =    0.0921642411;
 TSSMM.FX('102')              =    0.0889163929;
 TSSMM.FX('103')              =    0.0845852508;
 TSSMM.FX('104')              =    0.0790435486;
 TSSMM.FX('105')              =    0.0701257557;
 TSSMM.FX('106')              =    0.0564220179;
 TSSMM.FX('107')              =    0.0377816578;
 TSSMM.FX('108')              =    0.0141829724;
 TSSMM.FX('109')              =   -0.0104147192;
 TSSMM.FX('110')              =   -0.0203108952;
 TSSMM.FX('111')              =    0.0007605408;
 TSSMM.FX('112')              =    0.0103567437;
 TSSMM.FX('113')              =   -0.0139903269;
 TSSMM.FX('114')              =   -0.0489140824;
 TSSMM.FX('115')              =   -0.0960797219;
 TSSMM.FX('116')              =   -0.1298230961;
 TSSMM.FX('117')              =   -0.0587493683;
 TSSMM.FX('118')              =    0.0529727876;
 TSSMM.FX('119')              =    0.0837374333;
 TSSMM.FX('120')              =    0.0814930954;
 TSSMM.FX('121')              =    0.0802069188;
 TSSMM.FX('122')              =    0.0799759726;
 TSSMM.FX('123')              =    0.0798203258;
 TSSMM.FX('124')              =    0.0802215866;
 TSSMM.FX('125')              =    0.0807740107;
 TSSMM.FX('126')              =    0.0805771654;
 TSSMM.FX('127')              =    0.0790611659;
 TSSMM.FX('128')              =    0.0759577416;
 TSSMM.FX('129')              =    0.0717525920;
 TSSMM.FX('130')              =    0.0664915146;
 TSSMM.FX('131')              =    0.0612654945;
 TSSMM.FX('132')              =    0.0572833169;
 TSSMM.FX('133')              =    0.0541463155;
 TSSMM.FX('134')              =    0.0521551101;
 TSSMM.FX('135')              =    0.0508037639;
 TSSMM.FX('136')              =    0.0494396663;
 TSSMM.FX('137')              =    0.0478221366;
 TSSMM.FX('138')              =    0.0453707058;
 TSSMM.FX('139')              =    0.0412276241;
 TSSMM.FX('140')              =    0.0352018361;
 TSSMM.FX('141')              =    0.0280981178;
 TSSMM.FX('142')              =    0.0204503887;
 TSSMM.FX('143')              =    0.0121819807;
 TSSMM.FX('144')              =    0.0094185128;
 TSSMM.FX('145')              =    0.0096313074;
 TSSMM.FX('146')              =    0.0061874666;
 TSSMM.FX('147')              =    0.0116878290;
 TSSMM.FX('148')              =    0.0178154860;
 TSSMM.FX('149')              =    0.0138527861;
 TSSMM.FX('150')              =    0.0085564528;
 TSSMM.FX('151')              =    0.0006562254;
 TSSMM.FX('152')              =   -0.0100569169;
 TSSMM.FX('153')              =   -0.0221593222;
 TSSMM.FX('154')              =   -0.0235017899;
 TSSMM.FX('155')              =   -0.0204893929;
 TSSMM.FX('156')              =   -0.0267541078;
 TSSMM.FX('157')              =   -0.0347017056;
 TSSMM.FX('158')              =   -0.0430128148;
 TSSMM.FX('159')              =   -0.0530807201;
 TSSMM.FX('160')              =   -0.0675537435;
 TSSMM.FX('161')              =   -0.0888889088;
 TSSMM.FX('162')              =   -0.1194463459;
 TSSMM.FX('163')              =   -0.1640509522;
 TSSMM.FX('164')              =   -0.2182241088;
 TSSMM.FX('165')              =   -0.2415231265;
 TSSMM.FX('166')              =   -0.1808581676;
 TSSMM.FX('167')              =   -0.1414498173;
 TSSMM.FX('168')              =   -0.1840137246;
 TSSMM.FX('169')              =   -0.2185898968;
 TSSMM.FX('170')              =   -0.1815499613;
 TSSMM.FX('171')              =   -0.1076433269;
 TSSMM.FX('172')              =   -0.0914336041;
 TSSMM.FX('173')              =   -0.1128163118;
 TSSMM.FX('174')              =   -0.1400185877;
 TSSMM.FX('175')              =   -0.1751813533;
 TSSMM.FX('176')              =   -0.2200250925;
 TSSMM.FX('177')              =   -0.2569371289;
 TSSMM.FX('178')              =   -0.2347897334;
 TSSMM.FX('179')              =   -0.2131827840;
 TSSMM.FX('180')              =   -0.2520418064;
 TSSMM.FX('181')              =   -0.3116805774;
 TSSMM.FX('182')              =   -0.3948085960;
 TSSMM.FX('183')              =   -0.5151103556;
 TSSMM.FX('184')              =   -0.6485879635;
 TSSMM.FX('185')              =   -0.6112730681;
 TSSMM.FX('186')              =   -0.3541423597;
 TSSMM.FX('187')              =   -0.2207393558;
 TSSMM.FX('188')              =   -0.2652719083;
 TSSMM.FX('189')              =   -0.3094863178;
 TSSMM.FX('190')              =   -0.3784942835;
 TSSMM.FX('191')              =   -0.3818711291;
 TSSMM.FX('192')              =   -0.1645540808;
 TSSMM.FX('193')              =    0.0172182777;
 TSSMM.FX('194')              =    0.0179770685;
 TSSMM.FX('195')              =    0.0180497888;
 TSSMM.FX('196')              =    0.0172339971;
 TSSMM.FX('197')              =    0.0168906427;
 TSSMM.FX('198')              =    0.0166477820;
 TSSMM.FX('199')              =    0.0148213899;
 TSSMM.FX('200')              =    0.0125567208;
 TSSMM.FX('201')              =    0.0114703326;
 TSSMM.FX('202')              =    0.0106727138;
 TSSMM.FX('203')              =    0.0098808776;
 TSSMM.FX('204')              =    0.0087073706;
 TSSMM.FX('205')              =    0.0070986270;
 TSSMM.FX('206')              =    0.0044843572;
 TSSMM.FX('207')              =    0.0005903462;
 TSSMM.FX('208')              =   -0.0052320261;
 TSSMM.FX('209')              =   -0.0138874625;
 TSSMM.FX('210')              =   -0.0276525517;
 TSSMM.FX('211')              =   -0.0440173458;
 TSSMM.FX('212')              =   -0.0498532315;
 TSSMM.FX('213')              =   -0.0297679173;
 TSSMM.FX('214')              =   -0.0161487172;
 TSSMM.FX('215')              =   -0.0288609871;
 TSSMM.FX('216')              =   -0.0434793675;
 TSSMM.FX('217')              =   -0.0404953236;
 TSSMM.FX('218')              =   -0.0103624585;
 TSSMM.FX('219')              =    0.0078054325;
 TSSMM.FX('220')              =    0.0059511879;
 TSSMM.FX('221')              =    0.0033637151;
 TSSMM.FX('222')              =   -0.0024711395;
 TSSMM.FX('223')              =   -0.0116797829;
 TSSMM.FX('224')              =   -0.0218898713;
 TSSMM.FX('225')              =   -0.0157570302;
 TSSMM.FX('226')              =   -0.0036407370;
 TSSMM.FX('227')              =   -0.0004454429;
 TSSMM.FX('228')              =    0.0022228037;
 TSSMM.FX('229')              =    0.0026642087;
 TSSMM.FX('230')              =    0.0006095390;
 TSSMM.FX('231')              =   -0.0043177848;
 TSSMM.FX('232')              =   -0.0108741953;
 TSSMM.FX('233')              =   -0.0177875595;
 TSSMM.FX('234')              =   -0.0223066776;
 TSSMM.FX('235')              =   -0.0180521347;
 TSSMM.FX('236')              =   -0.0122155876;
 TSSMM.FX('237')              =   -0.0107467317;
 TSSMM.FX('238')              =   -0.0096493484;
 TSSMM.FX('239')              =   -0.0100122539;
 TSSMM.FX('240')              =   -0.0100517185;
 TSSMM.FX('241')              =   -0.0108745188;
 TSSMM.FX('242')              =   -0.0129393883;
 TSSMM.FX('243')              =   -0.0140220321;
 TSSMM.FX('244')              =   -0.0145600176;
 TSSMM.FX('245')              =   -0.0138472322;
 TSSMM.FX('246')              =   -0.0116874145;
 TSSMM.FX('247')              =   -0.0091677157;
 TSSMM.FX('248')              =   -0.0070332029;
 TSSMM.FX('249')              =   -0.0044473818;
 TSSMM.FX('250')              =   -0.0014371904;
 TSSMM.FX('251')              =    0.0000000000;

* Principal components of sea level rise
 PCSRES.FX('2000','1')        =    0.0188166691;
 PCSRES.FX('2000','2')        =    0.0120963853;

* Sea level rise in 1990 due to thermal expansion
 PCS1990.FX                   =    0.0238554517;

* Sea level rise in 2000 due to loss of mass of glaciers and small ice caps without glacier area contraction
 SLRGLWT.FX('2000')           =    0.0330393115;

* Sea level rise in 2000 due to loss of mass of glaciers and small ice caps with glacier area contraction
 SLRGLWC.FX('2000')           =    0.0295870074;

* Sea level rise in 1990 due to loss of mass of glaciers and small ice caps with glacier area contraction
 SLRGLWC1990.FX               =    0.0243687150;

* Greenland ice sheet contribution to sea level rise in 2000
 SLRGR.FX('2000')             =    0.0042314451;

* Antarctic ice sheet contribution to sea level rise in 2000
 SLRAN.FX('2000')             =   -0.0118148538;

* Sea level rise in 1990 due to four major factors relative to 1750
 SLR1990.FX                   =    0.0442241667;

