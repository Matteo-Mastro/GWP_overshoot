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
 BETA.FX('2000')              =    0.3000000000;

* Q10 for heterotrophic respiration
 Q10BAS.FX                    =    1.2000000000;

* Temperature sensitivity to Q10
* (not used in the standard inversion.)
 Q10TMP.FX                    =    0.0460000000;

* Temperature sensitivity to NPP
* (not used in the standard inversion.)
 TNPPSENS.FX                  =    0.0160000000;

* Preindustrial ocean CO2 uptake
 CO2OCNUPPRE.FX('2000')       =   -0.2912647760;

* Preindustrial land CO2 uptake
 CO2LNDUPPRE.FX('2000')       =    0.2198988101;

* Natural methane emission
 NATCH4.FX                    =  319.8989953789;

* Natural nitrous oxide emission
 NATN2O.FX                    =   11.3590617717;

* CH4 lifetime with respect to OH depletion
 TAUCH4OH.FX                  =    8.5347464087;

* N2O lifetime
 TAUN2O.FX                    =  114.1845055245;

* Radiative forcing of a doubling of atmospheric CO2 concentration
* (fixed in the standard inversion)
 RAD2XCO2.FX                  =    3.7000000000;

* Forcing scaling factor
* (fixed in the standard inversion)
 SCLRAD.FX                    =    1.1377978986;

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
 COC.FX('2000','1')           =  212.9978063292;
 COC.FX('2000','2')           =   55.1406282175;
 COC.FX('2000','3')           =   42.2779105857;
 COC.FX('2000','4')           =   13.3129422731;

* Carbon content change in the entire ocean in 1999
 SUMCOC1999.FX                =  123.9181997300;

* Carbon content change in the ocean mixed layer in 2000
 COC1MIX.FX('2000')           =   15.2620167484;

* Carbon content change in each layer of the land carbon cycle model in 2000
 CTR.FX('2000','1')           = -218.7470024157;
 CTR.FX('2000','2')           =  366.5518939307;
 CTR.FX('2000','3')           =  338.8984648669;
 CTR.FX('2000','4')           = 1809.6606358728;

* Carbon content change in the entire land biosphere in 1999
 SUMCTR1999.FX                = 2295.3538376545;

* Carbon content in the land biosphere in 1750
 CTR1750.FX('1')              = -205.2067684571;
 CTR1750.FX('2')              =  343.9893984777;
 CTR1750.FX('3')              =  321.9360000000;
 CTR1750.FX('4')              = 1759.3800000000;

* Carbon content change in the atmosphere in 2000
 COC1ATM.FX('2000')           =  197.7357895808;

* Cumulative anthropogenic CO2 emissions in 1999
 CUMCO21999.FX                =  378.3635436993;

* Atmospheric CO2 concentration in 2000
* (Do not transfer CON('2000','CO2') because of the different ENSO treatment in the past and future modes)
 CONCO2INT.FX('2000')         =  370.2087897014;

* Change in the atmospheric CO2 concentration due to ENSO in 2000
 DCONCO2ENSO2000.FX           =   -0.1731642676;

* Atmospheric concentrations of non-CO2 gases in 2000
 CON.FX('2000','CH4')         = 1752.0971683220;
 CON.FX('2000','N2O')         =  317.0380080874;
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
 RADAERAVE.FX                 =   -1.1448131922;

* Average missing forcing between 1901 and 2000
 RADMISAVE.FX                 =    0.0000000000;

* Principal components of temperature change
 PCTRES.FX('2000','1')        =    0.0266312329;
 PCTRES.FX('2000','2')        =    0.5332337859;

* Global-mean temperature change in 1750 (needed for temperature feedback to land CO2 uptake)
 T2M1750.FX                   =    0.0000000000;

* Global-mean temperature change in 1990 (needed for sea level calculation)
 T2M1990.FX                   =    0.6194503727;

* Land surface temperature in 1750 defined in DOECLIM
 TAL1750.FX                   =    0.0000000000;

* Land surface temperature in 2000 defined in DOECLIM
 TAL.FX('2000')               =    1.0066937590;

* Sea surface temperature in 1750 defined in DOECLIM
 TSS1750.FX                   =    0.0000000000;

* Sea surface temperature in 2000 defined in DOECLIM
 TSS.FX('2000')               =    0.6530870074;

* Temperature data transfer for DOECLIM memory calculation
* Note that the index is used in an opposite way of TSS(T).
 TSSMM.FX('1')                =    0.6530870074;
 TSSMM.FX('2')                =    0.6211402219;
 TSSMM.FX('3')                =    0.5888046277;
 TSSMM.FX('4')                =    0.5573998423;
 TSSMM.FX('5')                =    0.5267952415;
 TSSMM.FX('6')                =    0.4915596047;
 TSSMM.FX('7')                =    0.4472917341;
 TSSMM.FX('8')                =    0.4091710530;
 TSSMM.FX('9')                =    0.4378870700;
 TSSMM.FX('10')               =    0.5081246063;
 TSSMM.FX('11')               =    0.5210463873;
 TSSMM.FX('12')               =    0.4879850525;
 TSSMM.FX('13')               =    0.4547516327;
 TSSMM.FX('14')               =    0.4253564527;
 TSSMM.FX('15')               =    0.3967824372;
 TSSMM.FX('16')               =    0.3654476930;
 TSSMM.FX('17')               =    0.3366240999;
 TSSMM.FX('18')               =    0.3458853616;
 TSSMM.FX('19')               =    0.3918925496;
 TSSMM.FX('20')               =    0.4044425978;
 TSSMM.FX('21')               =    0.3790647447;
 TSSMM.FX('22')               =    0.3515998138;
 TSSMM.FX('23')               =    0.3248271863;
 TSSMM.FX('24')               =    0.2996966477;
 TSSMM.FX('25')               =    0.2809691288;
 TSSMM.FX('26')               =    0.2814621419;
 TSSMM.FX('27')               =    0.2829824468;
 TSSMM.FX('28')               =    0.2699890638;
 TSSMM.FX('29')               =    0.2533662216;
 TSSMM.FX('30')               =    0.2332179250;
 TSSMM.FX('31')               =    0.2133612229;
 TSSMM.FX('32')               =    0.2092222959;
 TSSMM.FX('33')               =    0.2149025294;
 TSSMM.FX('34')               =    0.2072409926;
 TSSMM.FX('35')               =    0.1866183168;
 TSSMM.FX('36')               =    0.1679168623;
 TSSMM.FX('37')               =    0.1850030399;
 TSSMM.FX('38')               =    0.2537882097;
 TSSMM.FX('39')               =    0.2988630077;
 TSSMM.FX('40')               =    0.2972184448;
 TSSMM.FX('41')               =    0.2905209587;
 TSSMM.FX('42')               =    0.2827311359;
 TSSMM.FX('43')               =    0.2736668036;
 TSSMM.FX('44')               =    0.2637770924;
 TSSMM.FX('45')               =    0.2574109519;
 TSSMM.FX('46')               =    0.2531196197;
 TSSMM.FX('47')               =    0.2526791085;
 TSSMM.FX('48')               =    0.2575739989;
 TSSMM.FX('49')               =    0.2617191225;
 TSSMM.FX('50')               =    0.2649127499;
 TSSMM.FX('51')               =    0.2653090624;
 TSSMM.FX('52')               =    0.2613035992;
 TSSMM.FX('53')               =    0.2533179685;
 TSSMM.FX('54')               =    0.2450388651;
 TSSMM.FX('55')               =    0.2408878878;
 TSSMM.FX('56')               =    0.2369945246;
 TSSMM.FX('57')               =    0.2338564772;
 TSSMM.FX('58')               =    0.2336974434;
 TSSMM.FX('59')               =    0.2342866517;
 TSSMM.FX('60')               =    0.2329571026;
 TSSMM.FX('61')               =    0.2282072593;
 TSSMM.FX('62')               =    0.2202260937;
 TSSMM.FX('63')               =    0.2105627107;
 TSSMM.FX('64')               =    0.2003367589;
 TSSMM.FX('65')               =    0.1882276411;
 TSSMM.FX('66')               =    0.1757355658;
 TSSMM.FX('67')               =    0.1649612934;
 TSSMM.FX('68')               =    0.1546013641;
 TSSMM.FX('69')               =    0.1449999530;
 TSSMM.FX('70')               =    0.1402334398;
 TSSMM.FX('71')               =    0.1377935284;
 TSSMM.FX('72')               =    0.1461140039;
 TSSMM.FX('73')               =    0.1574253310;
 TSSMM.FX('74')               =    0.1532681980;
 TSSMM.FX('75')               =    0.1447190579;
 TSSMM.FX('76')               =    0.1384941583;
 TSSMM.FX('77')               =    0.1332105524;
 TSSMM.FX('78')               =    0.1283461557;
 TSSMM.FX('79')               =    0.1309911099;
 TSSMM.FX('80')               =    0.1343745999;
 TSSMM.FX('81')               =    0.1295348879;
 TSSMM.FX('82')               =    0.1216322246;
 TSSMM.FX('83')               =    0.1091567719;
 TSSMM.FX('84')               =    0.0938473656;
 TSSMM.FX('85')               =    0.0765289590;
 TSSMM.FX('86')               =    0.0578397062;
 TSSMM.FX('87')               =    0.0401344176;
 TSSMM.FX('88')               =    0.0337707059;
 TSSMM.FX('89')               =    0.0438744107;
 TSSMM.FX('90')               =    0.0460060808;
 TSSMM.FX('91')               =    0.0338705915;
 TSSMM.FX('92')               =    0.0189532478;
 TSSMM.FX('93')               =    0.0011545423;
 TSSMM.FX('94')               =   -0.0145739516;
 TSSMM.FX('95')               =   -0.0383993184;
 TSSMM.FX('96')               =   -0.0750541170;
 TSSMM.FX('97')               =   -0.0927182219;
 TSSMM.FX('98')               =   -0.0169738297;
 TSSMM.FX('99')               =    0.0711133316;
 TSSMM.FX('100')              =    0.0842214928;
 TSSMM.FX('101')              =    0.0809234547;
 TSSMM.FX('102')              =    0.0760433346;
 TSSMM.FX('103')              =    0.0697704090;
 TSSMM.FX('104')              =    0.0619116973;
 TSSMM.FX('105')              =    0.0502437763;
 TSSMM.FX('106')              =    0.0335501897;
 TSSMM.FX('107')              =    0.0119952495;
 TSSMM.FX('108')              =   -0.0139325331;
 TSSMM.FX('109')              =   -0.0394751658;
 TSSMM.FX('110')              =   -0.0482195651;
 TSSMM.FX('111')              =   -0.0252312664;
 TSSMM.FX('112')              =   -0.0170396507;
 TSSMM.FX('113')              =   -0.0436140018;
 TSSMM.FX('114')              =   -0.0797052860;
 TSSMM.FX('115')              =   -0.1252742348;
 TSSMM.FX('116')              =   -0.1514275762;
 TSSMM.FX('117')              =   -0.0657149128;
 TSSMM.FX('118')              =    0.0540373756;
 TSSMM.FX('119')              =    0.0857830474;
 TSSMM.FX('120')              =    0.0830532104;
 TSSMM.FX('121')              =    0.0813700428;
 TSSMM.FX('122')              =    0.0806292913;
 TSSMM.FX('123')              =    0.0797709103;
 TSSMM.FX('124')              =    0.0793446005;
 TSSMM.FX('125')              =    0.0789163240;
 TSSMM.FX('126')              =    0.0775340043;
 TSSMM.FX('127')              =    0.0745854558;
 TSSMM.FX('128')              =    0.0699323778;
 TSSMM.FX('129')              =    0.0642216631;
 TSSMM.FX('130')              =    0.0574075012;
 TSSMM.FX('131')              =    0.0507110287;
 TSSMM.FX('132')              =    0.0453564748;
 TSSMM.FX('133')              =    0.0409173007;
 TSSMM.FX('134')              =    0.0376400930;
 TSSMM.FX('135')              =    0.0349137479;
 TSSMM.FX('136')              =    0.0320285593;
 TSSMM.FX('137')              =    0.0287065876;
 TSSMM.FX('138')              =    0.0242873668;
 TSSMM.FX('139')              =    0.0179108962;
 TSSMM.FX('140')              =    0.0095437052;
 TSSMM.FX('141')              =    0.0001735605;
 TSSMM.FX('142')              =   -0.0096033053;
 TSSMM.FX('143')              =   -0.0197297897;
 TSSMM.FX('144')              =   -0.0237572335;
 TSSMM.FX('145')              =   -0.0249558195;
 TSSMM.FX('146')              =   -0.0298446512;
 TSSMM.FX('147')              =   -0.0253864448;
 TSSMM.FX('148')              =   -0.0214755636;
 TSSMM.FX('149')              =   -0.0283229184;
 TSSMM.FX('150')              =   -0.0367481907;
 TSSMM.FX('151')              =   -0.0479961189;
 TSSMM.FX('152')              =   -0.0620376980;
 TSSMM.FX('153')              =   -0.0771639667;
 TSSMM.FX('154')              =   -0.0807346797;
 TSSMM.FX('155')              =   -0.0808398004;
 TSSMM.FX('156')              =   -0.0909584083;
 TSSMM.FX('157')              =   -0.1030578648;
 TSSMM.FX('158')              =   -0.1158958978;
 TSSMM.FX('159')              =   -0.1310031543;
 TSSMM.FX('160')              =   -0.1511271406;
 TSSMM.FX('161')              =   -0.1785566408;
 TSSMM.FX('162')              =   -0.2152639438;
 TSSMM.FX('163')              =   -0.2651894728;
 TSSMM.FX('164')              =   -0.3217159630;
 TSSMM.FX('165')              =   -0.3418752819;
 TSSMM.FX('166')              =   -0.2750456805;
 TSSMM.FX('167')              =   -0.2371825439;
 TSSMM.FX('168')              =   -0.2826623021;
 TSSMM.FX('169')              =   -0.3160560762;
 TSSMM.FX('170')              =   -0.2730282542;
 TSSMM.FX('171')              =   -0.1974683929;
 TSSMM.FX('172')              =   -0.1860293933;
 TSSMM.FX('173')              =   -0.2140519074;
 TSSMM.FX('174')              =   -0.2481063709;
 TSSMM.FX('175')              =   -0.2898905267;
 TSSMM.FX('176')              =   -0.3401927112;
 TSSMM.FX('177')              =   -0.3795305468;
 TSSMM.FX('178')              =   -0.3569509528;
 TSSMM.FX('179')              =   -0.3404886300;
 TSSMM.FX('180')              =   -0.3875842435;
 TSSMM.FX('181')              =   -0.4554445074;
 TSSMM.FX('182')              =   -0.5452133752;
 TSSMM.FX('183')              =   -0.6675516421;
 TSSMM.FX('184')              =   -0.7912668782;
 TSSMM.FX('185')              =   -0.7249681508;
 TSSMM.FX('186')              =   -0.4392873360;
 TSSMM.FX('187')              =   -0.3002046113;
 TSSMM.FX('188')              =   -0.3438049289;
 TSSMM.FX('189')              =   -0.3839490060;
 TSSMM.FX('190')              =   -0.4435449602;
 TSSMM.FX('191')              =   -0.4246562349;
 TSSMM.FX('192')              =   -0.1780385290;
 TSSMM.FX('193')              =    0.0119289483;
 TSSMM.FX('194')              =    0.0128241413;
 TSSMM.FX('195')              =    0.0128327183;
 TSSMM.FX('196')              =    0.0118430332;
 TSSMM.FX('197')              =    0.0112637763;
 TSSMM.FX('198')              =    0.0106398705;
 TSSMM.FX('199')              =    0.0082895856;
 TSSMM.FX('200')              =    0.0055349905;
 TSSMM.FX('201')              =    0.0039649668;
 TSSMM.FX('202')              =    0.0026021586;
 TSSMM.FX('203')              =    0.0011413684;
 TSSMM.FX('204')              =   -0.0008543741;
 TSSMM.FX('205')              =   -0.0034273812;
 TSSMM.FX('206')              =   -0.0072183029;
 TSSMM.FX('207')              =   -0.0124243369;
 TSSMM.FX('208')              =   -0.0197237346;
 TSSMM.FX('209')              =   -0.0298374324;
 TSSMM.FX('210')              =   -0.0448518430;
 TSSMM.FX('211')              =   -0.0613828468;
 TSSMM.FX('212')              =   -0.0657066827;
 TSSMM.FX('213')              =   -0.0431168429;
 TSSMM.FX('214')              =   -0.0294630567;
 TSSMM.FX('215')              =   -0.0426160252;
 TSSMM.FX('216')              =   -0.0565075901;
 TSSMM.FX('217')              =   -0.0506560562;
 TSSMM.FX('218')              =   -0.0175516626;
 TSSMM.FX('219')              =    0.0008245013;
 TSSMM.FX('220')              =   -0.0016670168;
 TSSMM.FX('221')              =   -0.0050932772;
 TSSMM.FX('222')              =   -0.0117700223;
 TSSMM.FX('223')              =   -0.0215910159;
 TSSMM.FX('224')              =   -0.0315392941;
 TSSMM.FX('225')              =   -0.0238266088;
 TSSMM.FX('226')              =   -0.0110070549;
 TSSMM.FX('227')              =   -0.0076791690;
 TSSMM.FX('228')              =   -0.0051886959;
 TSSMM.FX('229')              =   -0.0052409530;
 TSSMM.FX('230')              =   -0.0080103686;
 TSSMM.FX('231')              =   -0.0137256724;
 TSSMM.FX('232')              =   -0.0208572004;
 TSSMM.FX('233')              =   -0.0279548292;
 TSSMM.FX('234')              =   -0.0319953161;
 TSSMM.FX('235')              =   -0.0266996941;
 TSSMM.FX('236')              =   -0.0202861537;
 TSSMM.FX('237')              =   -0.0184743244;
 TSSMM.FX('238')              =   -0.0171487631;
 TSSMM.FX('239')              =   -0.0173529223;
 TSSMM.FX('240')              =   -0.0171726799;
 TSSMM.FX('241')              =   -0.0178578285;
 TSSMM.FX('242')              =   -0.0196463104;
 TSSMM.FX('243')              =   -0.0202742174;
 TSSMM.FX('244')              =   -0.0201886038;
 TSSMM.FX('245')              =   -0.0186699999;
 TSSMM.FX('246')              =   -0.0156207854;
 TSSMM.FX('247')              =   -0.0122285184;
 TSSMM.FX('248')              =   -0.0092516620;
 TSSMM.FX('249')              =   -0.0058091886;
 TSSMM.FX('250')              =   -0.0020370928;
 TSSMM.FX('251')              =    0.0000000000;

* Principal components of sea level rise
 PCSRES.FX('2000','1')        =    0.0193153436;
 PCSRES.FX('2000','2')        =    0.0147577623;

* Sea level rise in 1990 due to thermal expansion
 PCS1990.FX                   =    0.0248887934;

* Sea level rise in 2000 due to loss of mass of glaciers and small ice caps without glacier area contraction
 SLRGLWT.FX('2000')           =    0.0330775504;

* Sea level rise in 2000 due to loss of mass of glaciers and small ice caps with glacier area contraction
 SLRGLWC.FX('2000')           =    0.0296197772;

* Sea level rise in 1990 due to loss of mass of glaciers and small ice caps with glacier area contraction
 SLRGLWC1990.FX               =    0.0243687150;

* Greenland ice sheet contribution to sea level rise in 2000
 SLRGR.FX('2000')             =    0.0042331860;

* Antarctic ice sheet contribution to sea level rise in 2000
 SLRAN.FX('2000')             =   -0.0118435480;

* Sea level rise in 1990 due to four major factors relative to 1750
 SLR1990.FX                   =    0.0452575084;

