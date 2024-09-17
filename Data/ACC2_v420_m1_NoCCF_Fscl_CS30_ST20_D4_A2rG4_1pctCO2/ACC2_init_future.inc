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
* BETA.FX('2000')              =    0.6;
 BETA.FX('2000')              =    0.5239643167;
 
* Q10 for heterotrophic respiration
* Q10BAS.FX                    =    1.2000000000;
 Q10BAS.FX                    =    2.0000000000;

* Temperature sensitivity to Q10
* (not used in the standard inversion.)
 Q10TMP.FX                    =    0.0460000000;

* Temperature sensitivity to NPP
* (not used in the standard inversion.)
 TNPPSENS.FX                  =    0.0160000000;

* Preindustrial ocean CO2 uptake
 CO2OCNUPPRE.FX('2000')       =   -0.2400000000;

* Preindustrial land CO2 uptake
 CO2LNDUPPRE.FX('2000')       =    0.3000000000;

* Natural methane emission
 NATCH4.FX                    =  226.3089521219;

* Natural nitrous oxide emission
 NATN2O.FX                    =   10.8837305766;

* CH4 lifetime with respect to OH depletion
 TAUCH4OH.FX                  =    9.7019949431;

* N2O lifetime
 TAUN2O.FX                    =  117.9826946717;

* Radiative forcing of a doubling of atmospheric CO2 concentration
* (fixed in the standard inversion)
 RAD2XCO2.FX                  =    3.7000000000;

* Forcing scaling factor
* (fixed in the standard inversion)
 SCLRAD.FX                    =    0.9990283568;

* Climate sensitivity
* SENS.FX('2000')              =    4.0000000000;
 SENS.FX('2000')              =    3.0000000000;

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
 COC.FX('2000','1')           =    0.0000000000;
 COC.FX('2000','2')           =    0.0000000000;
 COC.FX('2000','3')           =    0.0000000000;
 COC.FX('2000','4')           =    0.0000000000;

* Carbon content change in the entire ocean in 1999
 SUMCOC1999.FX                =    0.0000000000;

* Carbon content change in the ocean mixed layer in 2000
 COC1MIX.FX('2000')           =    0.0000000000;

* Carbon content change in each layer of the land carbon cycle model in 2000
 CTR.FX('2000','1')           = -205.2067684571;
 CTR.FX('2000','2')           =  343.9893984777;
 CTR.FX('2000','3')           =  321.9360000000;
 CTR.FX('2000','4')           = 1759.3800000000;

* Carbon content change in the entire land biosphere in 1999
 SUMCTR1999.FX                = 2220.0986300205;

* Carbon content in the land biosphere in 1750
 CTR1750.FX('1')              = -205.2067684571;
 CTR1750.FX('2')              =  343.9893984777;
 CTR1750.FX('3')              =  321.9360000000;
 CTR1750.FX('4')              = 1759.3800000000;

* Carbon content change in the atmosphere in 2000
 COC1ATM.FX('2000')           =    0.0000000000;

* Cumulative anthropogenic CO2 emissions in 1999
 CUMCO21999.FX                =    0.0000000000;

* Atmospheric CO2 concentration in 2000
* (Do not transfer CON('2000','CO2') because of the different ENSO treatment in the past and future modes)
 CONCO2INT.FX('2000')         =  277.0000000000;

* Change in the atmospheric CO2 concentration due to ENSO in 2000
 DCONCO2ENSO2000.FX           =    0.0000000000;

* Atmospheric concentrations of non-CO2 gases in 2000
 CON.FX('2000','CH4')         =  699.9999999935;
 CON.FX('2000','N2O')         =  270.0000000001;
 CON.FX('2000','CF4')         =   39.0000000000;
 CON.FX('2000','C2F6')        =    0.0000000000;
 CON.FX('2000','C4F10')       =    0.0000000000;
 CON.FX('2000','HFC23')       =    0.0000000000;
 CON.FX('2000','HFC32')       =    0.0000000000;
 CON.FX('2000','HFC4310mee')  =    0.0000000000;
 CON.FX('2000','HFC125')      =    0.0000000000;
 CON.FX('2000','HFC134a')     =    0.0000000000;
 CON.FX('2000','HFC143a')     =    0.0000000000;
 CON.FX('2000','HFC152a')     =    0.0000000000;
 CON.FX('2000','HFC227ea')    =    0.0000000000;
 CON.FX('2000','HFC236fa')    =    0.0000000000;
 CON.FX('2000','HFC245ca')    =    0.0000000000;
 CON.FX('2000','CFC11')       =    0.0000000000;
 CON.FX('2000','CFC12')       =    0.0000000000;
 CON.FX('2000','CFC113')      =    0.0000000000;
 CON.FX('2000','CFC114')      =    0.0000000000;
 CON.FX('2000','CFC115')      =    0.0000000000;
 CON.FX('2000','CCl4')        =    0.0000000000;
 CON.FX('2000','CH3CCl3')     =    0.0000000000;
 CON.FX('2000','Halon1211')   =    0.0000000000;
 CON.FX('2000','Halon1301')   =    0.0000000000;
 CON.FX('2000','Halon2402')   =    0.0000000000;
 CON.FX('2000','HCFC22')      =    0.0000000000;
 CON.FX('2000','HCFC141b')    =    0.0000000000;
 CON.FX('2000','HCFC142b')    =    0.0000000000;
 CON.FX('2000','HCFC123')     =    0.0000000000;
 CON.FX('2000','CH3Cl')       =  550.59;
 CON.FX('2000','CH3Br')       =  8.2454;
 CON.FX('2000','SF6')         =    0.0100000000;
 CON.FX('2000','O3TRP')       =  25.0;

* Concentrations prior to 2000 used for the calculation of EESCl
 CON1999.FX('CFC11')          =    0.0000000000;
 CON1999.FX('CFC12')          =    0.0000000000;
 CON1999.FX('CFC113')         =    0.0000000000;
 CON1999.FX('CFC114')         =    0.0000000000;
 CON1999.FX('CFC115')         =    0.0000000000;
 CON1999.FX('CCl4')           =    0.0000000000;
 CON1999.FX('CH3CCl3')        =    0.0000000000;
 CON1999.FX('Halon1211')      =    0.0000000000;
 CON1999.FX('Halon1301')      =    0.0000000000;
 CON1999.FX('Halon2402')      =    0.0000000000;
 CON1999.FX('HCFC22')         =    0.0000000000;
 CON1999.FX('HCFC141b')       =    0.0000000000;
 CON1999.FX('HCFC142b')       =    0.0000000000;
 CON1999.FX('HCFC123')        =    0.0000000000;
 CON1999.FX('CH3Cl')          =  550.59;
 CON1999.FX('CH3Br')          =  8.2454;

 CON1998.FX('CFC11')          =    0.0000000000;
 CON1998.FX('CFC12')          =    0.0000000000;
 CON1998.FX('CFC113')         =    0.0000000000;
 CON1998.FX('CFC114')         =    0.0000000000;
 CON1998.FX('CFC115')         =    0.0000000000;
 CON1998.FX('CCl4')           =    0.0000000000;
 CON1998.FX('CH3CCl3')        =    0.0000000000;
 CON1998.FX('Halon1211')      =    0.0000000000;
 CON1998.FX('Halon1301')      =    0.0000000000;
 CON1998.FX('Halon2402')      =    0.0000000000;
 CON1998.FX('HCFC22')         =    0.0000000000;
 CON1998.FX('HCFC141b')       =    0.0000000000;
 CON1998.FX('HCFC142b')       =    0.0000000000;
 CON1998.FX('HCFC123')        =    0.0000000000;
 CON1998.FX('CH3Cl')          =  550.59;
 CON1998.FX('CH3Br')          =  8.2454;

 CON1997.FX('CFC11')          =    0.0000000000;
 CON1997.FX('CFC12')          =    0.0000000000;
 CON1997.FX('CFC113')         =    0.0000000000;
 CON1997.FX('CFC114')         =    0.0000000000;
 CON1997.FX('CFC115')         =    0.0000000000;
 CON1997.FX('CCl4')           =    0.0000000000;
 CON1997.FX('CH3CCl3')        =    0.0000000000;
 CON1997.FX('Halon1211')      =    0.0000000000;
 CON1997.FX('Halon1301')      =    0.0000000000;
 CON1997.FX('Halon2402')      =    0.0000000000;
 CON1997.FX('HCFC22')         =    0.0000000000;
 CON1997.FX('HCFC141b')       =    0.0000000000;
 CON1997.FX('HCFC142b')       =    0.0000000000;
 CON1997.FX('HCFC123')        =    0.0000000000;
 CON1997.FX('CH3Cl')          =  550.59;
 CON1997.FX('CH3Br')          =  8.2454;

* Solar forcing in 2000
 RADSOL2000.FX                =    0.0000000000;

* Solar forcing averaged between 1901 and 2000
 RADSOL20AVE.FX               =    0.0000000000;

* Average aerosol forcing between 1901 and 2000
 RADAERAVE.FX                 =    0.0000000000;

* Average missing forcing between 1901 and 2000
 RADMISAVE.FX                 =    0.0000000000;

* Principal components of temperature change
 PCTRES.FX('2000','1')        =    0.0000000000;
 PCTRES.FX('2000','2')        =    0.0000000000;

* Global-mean temperature change in 1750 (needed for temperature feedback to land CO2 uptake)
 T2M1750.FX                   =    0.0000000000;

* Global-mean temperature change in 1990 (needed for sea level calculation)
 T2M1990.FX                   =    0.0000000000;

* Land surface temperature in 1750 defined in DOECLIM
 TAL1750.FX                   =    0.0000000000;

* Land surface temperature in 2000 defined in DOECLIM
 TAL.FX('2000')               =    0.0000000000;

* Sea surface temperature in 1750 defined in DOECLIM
 TSS1750.FX                   =    0.0000000000;

* Sea surface temperature in 2000 defined in DOECLIM
 TSS.FX('2000')               =    0.0000000000;

* Temperature data transfer for DOECLIM memory calculation
* Note that the index is used in an opposite way of TSS(T).
 TSSMM.FX('1')                =    0.0000000000;
 TSSMM.FX('2')                =    0.0000000000;
 TSSMM.FX('3')                =    0.0000000000;
 TSSMM.FX('4')                =    0.0000000000;
 TSSMM.FX('5')                =    0.0000000000;
 TSSMM.FX('6')                =    0.0000000000;
 TSSMM.FX('7')                =    0.0000000000;
 TSSMM.FX('8')                =    0.0000000000;
 TSSMM.FX('9')                =    0.0000000000;
 TSSMM.FX('10')               =    0.0000000000;
 TSSMM.FX('11')               =    0.0000000000;
 TSSMM.FX('12')               =    0.0000000000;
 TSSMM.FX('13')               =    0.0000000000;
 TSSMM.FX('14')               =    0.0000000000;
 TSSMM.FX('15')               =    0.0000000000;
 TSSMM.FX('16')               =    0.0000000000;
 TSSMM.FX('17')               =    0.0000000000;
 TSSMM.FX('18')               =    0.0000000000;
 TSSMM.FX('19')               =    0.0000000000;
 TSSMM.FX('20')               =    0.0000000000;
 TSSMM.FX('21')               =    0.0000000000;
 TSSMM.FX('22')               =    0.0000000000;
 TSSMM.FX('23')               =    0.0000000000;
 TSSMM.FX('24')               =    0.0000000000;
 TSSMM.FX('25')               =    0.0000000000;
 TSSMM.FX('26')               =    0.0000000000;
 TSSMM.FX('27')               =    0.0000000000;
 TSSMM.FX('28')               =    0.0000000000;
 TSSMM.FX('29')               =    0.0000000000;
 TSSMM.FX('30')               =    0.0000000000;
 TSSMM.FX('31')               =    0.0000000000;
 TSSMM.FX('32')               =    0.0000000000;
 TSSMM.FX('33')               =    0.0000000000;
 TSSMM.FX('34')               =    0.0000000000;
 TSSMM.FX('35')               =    0.0000000000;
 TSSMM.FX('36')               =    0.0000000000;
 TSSMM.FX('37')               =    0.0000000000;
 TSSMM.FX('38')               =    0.0000000000;
 TSSMM.FX('39')               =    0.0000000000;
 TSSMM.FX('40')               =    0.0000000000;
 TSSMM.FX('41')               =    0.0000000000;
 TSSMM.FX('42')               =    0.0000000000;
 TSSMM.FX('43')               =    0.0000000000;
 TSSMM.FX('44')               =    0.0000000000;
 TSSMM.FX('45')               =    0.0000000000;
 TSSMM.FX('46')               =    0.0000000000;
 TSSMM.FX('47')               =    0.0000000000;
 TSSMM.FX('48')               =    0.0000000000;
 TSSMM.FX('49')               =    0.0000000000;
 TSSMM.FX('50')               =    0.0000000000;
 TSSMM.FX('51')               =    0.0000000000;
 TSSMM.FX('52')               =    0.0000000000;
 TSSMM.FX('53')               =    0.0000000000;
 TSSMM.FX('54')               =    0.0000000000;
 TSSMM.FX('55')               =    0.0000000000;
 TSSMM.FX('56')               =    0.0000000000;
 TSSMM.FX('57')               =    0.0000000000;
 TSSMM.FX('58')               =    0.0000000000;
 TSSMM.FX('59')               =    0.0000000000;
 TSSMM.FX('60')               =    0.0000000000;
 TSSMM.FX('61')               =    0.0000000000;
 TSSMM.FX('62')               =    0.0000000000;
 TSSMM.FX('63')               =    0.0000000000;
 TSSMM.FX('64')               =    0.0000000000;
 TSSMM.FX('65')               =    0.0000000000;
 TSSMM.FX('66')               =    0.0000000000;
 TSSMM.FX('67')               =    0.0000000000;
 TSSMM.FX('68')               =    0.0000000000;
 TSSMM.FX('69')               =    0.0000000000;
 TSSMM.FX('70')               =    0.0000000000;
 TSSMM.FX('71')               =    0.0000000000;
 TSSMM.FX('72')               =    0.0000000000;
 TSSMM.FX('73')               =    0.0000000000;
 TSSMM.FX('74')               =    0.0000000000;
 TSSMM.FX('75')               =    0.0000000000;
 TSSMM.FX('76')               =    0.0000000000;
 TSSMM.FX('77')               =    0.0000000000;
 TSSMM.FX('78')               =    0.0000000000;
 TSSMM.FX('79')               =    0.0000000000;
 TSSMM.FX('80')               =    0.0000000000;
 TSSMM.FX('81')               =    0.0000000000;
 TSSMM.FX('82')               =    0.0000000000;
 TSSMM.FX('83')               =    0.0000000000;
 TSSMM.FX('84')               =    0.0000000000;
 TSSMM.FX('85')               =    0.0000000000;
 TSSMM.FX('86')               =    0.0000000000;
 TSSMM.FX('87')               =    0.0000000000;
 TSSMM.FX('88')               =    0.0000000000;
 TSSMM.FX('89')               =    0.0000000000;
 TSSMM.FX('90')               =    0.0000000000;
 TSSMM.FX('91')               =    0.0000000000;
 TSSMM.FX('92')               =    0.0000000000;
 TSSMM.FX('93')               =    0.0000000000;
 TSSMM.FX('94')               =    0.0000000000;
 TSSMM.FX('95')               =    0.0000000000;
 TSSMM.FX('96')               =    0.0000000000;
 TSSMM.FX('97')               =    0.0000000000;
 TSSMM.FX('98')               =    0.0000000000;
 TSSMM.FX('99')               =    0.0000000000;
 TSSMM.FX('100')              =    0.0000000000;
 TSSMM.FX('101')              =    0.0000000000;
 TSSMM.FX('102')              =    0.0000000000;
 TSSMM.FX('103')              =    0.0000000000;
 TSSMM.FX('104')              =    0.0000000000;
 TSSMM.FX('105')              =    0.0000000000;
 TSSMM.FX('106')              =    0.0000000000;
 TSSMM.FX('107')              =    0.0000000000;
 TSSMM.FX('108')              =    0.0000000000;
 TSSMM.FX('109')              =    0.0000000000;
 TSSMM.FX('110')              =    0.0000000000;
 TSSMM.FX('111')              =    0.0000000000;
 TSSMM.FX('112')              =    0.0000000000;
 TSSMM.FX('113')              =    0.0000000000;
 TSSMM.FX('114')              =    0.0000000000;
 TSSMM.FX('115')              =    0.0000000000;
 TSSMM.FX('116')              =    0.0000000000;
 TSSMM.FX('117')              =    0.0000000000;
 TSSMM.FX('118')              =    0.0000000000;
 TSSMM.FX('119')              =    0.0000000000;
 TSSMM.FX('120')              =    0.0000000000;
 TSSMM.FX('121')              =    0.0000000000;
 TSSMM.FX('122')              =    0.0000000000;
 TSSMM.FX('123')              =    0.0000000000;
 TSSMM.FX('124')              =    0.0000000000;
 TSSMM.FX('125')              =    0.0000000000;
 TSSMM.FX('126')              =    0.0000000000;
 TSSMM.FX('127')              =    0.0000000000;
 TSSMM.FX('128')              =    0.0000000000;
 TSSMM.FX('129')              =    0.0000000000;
 TSSMM.FX('130')              =    0.0000000000;
 TSSMM.FX('131')              =    0.0000000000;
 TSSMM.FX('132')              =    0.0000000000;
 TSSMM.FX('133')              =    0.0000000000;
 TSSMM.FX('134')              =    0.0000000000;
 TSSMM.FX('135')              =    0.0000000000;
 TSSMM.FX('136')              =    0.0000000000;
 TSSMM.FX('137')              =    0.0000000000;
 TSSMM.FX('138')              =    0.0000000000;
 TSSMM.FX('139')              =    0.0000000000;
 TSSMM.FX('140')              =    0.0000000000;
 TSSMM.FX('141')              =    0.0000000000;
 TSSMM.FX('142')              =    0.0000000000;
 TSSMM.FX('143')              =    0.0000000000;
 TSSMM.FX('144')              =    0.0000000000;
 TSSMM.FX('145')              =    0.0000000000;
 TSSMM.FX('146')              =    0.0000000000;
 TSSMM.FX('147')              =    0.0000000000;
 TSSMM.FX('148')              =    0.0000000000;
 TSSMM.FX('149')              =    0.0000000000;
 TSSMM.FX('150')              =    0.0000000000;
 TSSMM.FX('151')              =    0.0000000000;
 TSSMM.FX('152')              =    0.0000000000;
 TSSMM.FX('153')              =    0.0000000000;
 TSSMM.FX('154')              =    0.0000000000;
 TSSMM.FX('155')              =    0.0000000000;
 TSSMM.FX('156')              =    0.0000000000;
 TSSMM.FX('157')              =    0.0000000000;
 TSSMM.FX('158')              =    0.0000000000;
 TSSMM.FX('159')              =    0.0000000000;
 TSSMM.FX('160')              =    0.0000000000;
 TSSMM.FX('161')              =    0.0000000000;
 TSSMM.FX('162')              =    0.0000000000;
 TSSMM.FX('163')              =    0.0000000000;
 TSSMM.FX('164')              =    0.0000000000;
 TSSMM.FX('165')              =    0.0000000000;
 TSSMM.FX('166')              =    0.0000000000;
 TSSMM.FX('167')              =    0.0000000000;
 TSSMM.FX('168')              =    0.0000000000;
 TSSMM.FX('169')              =    0.0000000000;
 TSSMM.FX('170')              =    0.0000000000;
 TSSMM.FX('171')              =    0.0000000000;
 TSSMM.FX('172')              =    0.0000000000;
 TSSMM.FX('173')              =    0.0000000000;
 TSSMM.FX('174')              =    0.0000000000;
 TSSMM.FX('175')              =    0.0000000000;
 TSSMM.FX('176')              =    0.0000000000;
 TSSMM.FX('177')              =    0.0000000000;
 TSSMM.FX('178')              =    0.0000000000;
 TSSMM.FX('179')              =    0.0000000000;
 TSSMM.FX('180')              =    0.0000000000;
 TSSMM.FX('181')              =    0.0000000000;
 TSSMM.FX('182')              =    0.0000000000;
 TSSMM.FX('183')              =    0.0000000000;
 TSSMM.FX('184')              =    0.0000000000;
 TSSMM.FX('185')              =    0.0000000000;
 TSSMM.FX('186')              =    0.0000000000;
 TSSMM.FX('187')              =    0.0000000000;
 TSSMM.FX('188')              =    0.0000000000;
 TSSMM.FX('189')              =    0.0000000000;
 TSSMM.FX('190')              =    0.0000000000;
 TSSMM.FX('191')              =    0.0000000000;
 TSSMM.FX('192')              =    0.0000000000;
 TSSMM.FX('193')              =    0.0000000000;
 TSSMM.FX('194')              =    0.0000000000;
 TSSMM.FX('195')              =    0.0000000000;
 TSSMM.FX('196')              =    0.0000000000;
 TSSMM.FX('197')              =    0.0000000000;
 TSSMM.FX('198')              =    0.0000000000;
 TSSMM.FX('199')              =    0.0000000000;
 TSSMM.FX('200')              =    0.0000000000;
 TSSMM.FX('201')              =    0.0000000000;
 TSSMM.FX('202')              =    0.0000000000;
 TSSMM.FX('203')              =    0.0000000000;
 TSSMM.FX('204')              =    0.0000000000;
 TSSMM.FX('205')              =    0.0000000000;
 TSSMM.FX('206')              =    0.0000000000;
 TSSMM.FX('207')              =    0.0000000000;
 TSSMM.FX('208')              =    0.0000000000;
 TSSMM.FX('209')              =    0.0000000000;
 TSSMM.FX('210')              =    0.0000000000;
 TSSMM.FX('211')              =    0.0000000000;
 TSSMM.FX('212')              =    0.0000000000;
 TSSMM.FX('213')              =    0.0000000000;
 TSSMM.FX('214')              =    0.0000000000;
 TSSMM.FX('215')              =    0.0000000000;
 TSSMM.FX('216')              =    0.0000000000;
 TSSMM.FX('217')              =    0.0000000000;
 TSSMM.FX('218')              =    0.0000000000;
 TSSMM.FX('219')              =    0.0000000000;
 TSSMM.FX('220')              =    0.0000000000;
 TSSMM.FX('221')              =    0.0000000000;
 TSSMM.FX('222')              =    0.0000000000;
 TSSMM.FX('223')              =    0.0000000000;
 TSSMM.FX('224')              =    0.0000000000;
 TSSMM.FX('225')              =    0.0000000000;
 TSSMM.FX('226')              =    0.0000000000;
 TSSMM.FX('227')              =    0.0000000000;
 TSSMM.FX('228')              =    0.0000000000;
 TSSMM.FX('229')              =    0.0000000000;
 TSSMM.FX('230')              =    0.0000000000;
 TSSMM.FX('231')              =    0.0000000000;
 TSSMM.FX('232')              =    0.0000000000;
 TSSMM.FX('233')              =    0.0000000000;
 TSSMM.FX('234')              =    0.0000000000;
 TSSMM.FX('235')              =    0.0000000000;
 TSSMM.FX('236')              =    0.0000000000;
 TSSMM.FX('237')              =    0.0000000000;
 TSSMM.FX('238')              =    0.0000000000;
 TSSMM.FX('239')              =    0.0000000000;
 TSSMM.FX('240')              =    0.0000000000;
 TSSMM.FX('241')              =    0.0000000000;
 TSSMM.FX('242')              =    0.0000000000;
 TSSMM.FX('243')              =    0.0000000000;
 TSSMM.FX('244')              =    0.0000000000;
 TSSMM.FX('245')              =    0.0000000000;
 TSSMM.FX('246')              =    0.0000000000;
 TSSMM.FX('247')              =    0.0000000000;
 TSSMM.FX('248')              =    0.0000000000;
 TSSMM.FX('249')              =    0.0000000000;
 TSSMM.FX('250')              =    0.0000000000;
 TSSMM.FX('251')              =    0.0000000000;

* Principal components of sea level rise
 PCSRES.FX('2000','1')        =    0.0000000000;
 PCSRES.FX('2000','2')        =    0.0000000000;

* Sea level rise in 1990 due to thermal expansion
 PCS1990.FX                   =    0.0000000000;

* Sea level rise in 2000 due to loss of mass of glaciers and small ice caps without glacier area contraction
 SLRGLWT.FX('2000')           =    0.0329241000;

* Sea level rise in 2000 due to loss of mass of glaciers and small ice caps with glacier area contraction
 SLRGLWC.FX('2000')           =    0.0294882536;

* Sea level rise in 1990 due to loss of mass of glaciers and small ice caps with glacier area contraction
 SLRGLWC1990.FX               =    0.0243687150;

* Greenland ice sheet contribution to sea level rise in 2000
 SLRGR.FX('2000')             =    0.0042262000;

* Antarctic ice sheet contribution to sea level rise in 2000
 SLRAN.FX('2000')             =   -0.0117284000;

* Sea level rise in 1990 due to four major factors relative to 1750
 SLR1990.FX                   =    0.0203687150;

