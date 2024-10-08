********** ACC2 version 4.2: ACC2_bounds_future.inc ********************************************************************
********** Levels and bounds for valuables (future mode) ===============================================================

* Note that, when upper or lower bounds are added or removed below,
* they should be reflected in "OUT_FUTURE_ALERT" section in "ACC2_main_future.gms".

H.LO(T)         = 3;
H.L(T)          = 8.6;
H.UP(T)         = 100;

COC.LO(T,NOC)   = -100;
COC.L(T,NOC)    = 100;
COC.UP(T,NOC)   = 10000;

COC1ATM.LO(T)   = -100;
COC1ATM.L(T)    = 200;
COC1ATM.UP(T)   = 10000;

COC1MIX.LO(T)   = -100;
COC1MIX.L(T)    = 15;
COC1MIX.UP(T)   = 10000;

CTR.LO(T,NTR)   = -10000;
CTR.L(T,NTR)    = 10;
CTR.UP(T,NTR)   = 10000;

DIC.LO(T)       = 1.000E-6;
DIC.L(T)        = 2.000E-3;
DIC.UP(T)       = 1.000E-1;

REVELLE.LO(T)   = 1;
REVELLE.L(T)    = 10;
REVELLE.UP(T)   = 100;

AK1S.LO(T)      = 5.000E2;
AK1S.L(T)       = 1.240E3;
AK1S.UP(T)      = 3.000E3;

AK1P.LO(T)      = 1.000E2;
AK1P.L(T)       = 4.000E2;
AK1P.UP(T)      = 1.500E3;

AK2S.LO(T)      = 1.000E-1;
AK2S.L(T)       = 8.270E-1;
AK2S.UP(T)      = 1.500E0;

AK2P.LO(T)      = 5.000E-3;
AK2P.L(T)       = 3.990E-2;
AK2P.UP(T)      = 5.000E-1;

AKB.LO(T)       = 5.000E-1;
AKB.L(T)        = 1.881E0;
AKB.UP(T)       = 5.000E0;

AKW.LO(T)       = 1.000E3;
AKW.L(T)        = 2.506E4;
AKW.UP(T)       = 3.000E5;

ALPHAS.LO(T)    = 1.000E-2;
ALPHAS.L(T)     = 2.988E-2;
ALPHAS.UP(T)    = 5.000E-2;

COC1ALL.LO(T)   = 100;
COC1ALL.L(T)    = 1300;
COC1ALL.UP(T)   = 10000;

*TMIX.LO(T)      = 0;
TMIX.L(T)       = 19.3;
*TMIX.UP(T)      = 40;
* Adding bounds for TMIX lead to infeasibility for GAMSCHK.

TNPPFAC.LO(T)   = 0.8;
TNPPFAC.L(T)    = 1.0;
TNPPFAC.UP(T)   = 1.2;

TRESFAC.LO(T)   = 0.05;
TRESFAC.L(T)    = 1;
TRESFAC.UP(T)   = 6;

CONCO2INT.LO(T) = 150;
CONCO2INT.L(T)  = 380;
CONCO2INT.UP(T) = 5000;

*ANT.LO(T,'CO2') = -5;
ANT.LO(T,'CO2') = -50;
ANT.L(T,'CO2')  = 10;
ANT.UP(T,'CO2') = 50;

ANT.LO(T,'CH4') = -10;
ANT.L(T,'CH4')  = 350;
ANT.UP(T,'CH4') = 3000;

ANT.LO(T,'N2O') = -5;
ANT.L(T,'N2O')  = 7;
ANT.UP(T,'N2O') = 50;

ANT.L(T,'CO')   = 0.001;
* To prevent a zero division error

ANT.LO(T,'SO2') = 0;
ANT.L(T,'SO2')  = 70;
ANT.UP(T,'SO2') = 300;

CON.LO(T,'CO2') = 150;
CON.L(T,'CO2')  = 380;
CON.UP(T,'CO2') = 5000;

*CON.LO(T,'CH4') = 700;
CON.L(T,'CH4')  = 1760;
CON.UP(T,'CH4') = 10000;

CON.LO(T,'N2O') = 150;
CON.L(T,'N2O')  = 316;
CON.UP(T,'N2O') = 1000;

CON.LO(T,'rOH') = 0.1;
CON.L(T,'rOH')  = 1;
CON.UP(T,'rOH') = 10;

CON.LO(T,'CFC11')  = 0;
CON.L(T,'CFC11')   = 262;
CON.UP(T,'CFC11')  = 1000;

CON.LO(T,'CFC12')  = 0;
CON.L(T,'CFC12')   = 541;
CON.UP(T,'CFC12')  = 2000;

CON.LO(T,'CFC113') = 0;
CON.L(T,'CFC113')  = 81;
CON.UP(T,'CFC113') = 500;

CON.LO(T,'CFC114') = 0;
CON.L(T,'CFC114')  = 16.4;
CON.UP(T,'CFC114') = 50;

CON.LO(T,'CFC115') = 0;
CON.L(T,'CFC115')  = 9;
CON.UP(T,'CFC115') = 50;

CON.LO(T,'CCl4')   = 0;
CON.L(T,'CCl4')    = 97.9;
CON.UP(T,'CCl4')   = 500;

CON.LO(T,'HCFC22') = 0;
CON.L(T,'HCFC22')  = 141.6;
CON.UP(T,'HCFC22') = 1000;

CON.LO(T,'HCFC141b') = 0;
CON.L(T,'HCFC141b')  = 10.8;
CON.UP(T,'HCFC141b') = 100;

CON.LO(T,'HCFC142b') = 0;
CON.L(T,'HCFC142b')  = 10.5;
CON.UP(T,'HCFC142b') = 100;

CON.LO(T,'HCFC123')  = 0;
CON.L(T,'HCFC123')   = 0;
CON.UP(T,'HCFC123')  = 50;

RAD.LO(T,FOA)   = -20;
RAD.L(T,FOA)    = 0.1;
RAD.UP(T,FOA)   = 30;

RADTTL.LO(T)    = -20;
RADTTL.L(T)     = 1;
RADTTL.UP(T)    = 30;

RADHAL.LO(T)    = -5;
RADHAL.L(T)     = 1;
RADHAL.UP(T)    = 5;

RADODS.LO(T)    = -5;
RADODS.L(T)     = 1;
RADODS.UP(T)    = 5;

RADSOL.LO(T)    = -3;
RADSOL.L(T)     = 0;
RADSOL.UP(T)    = 3;

NORRAD.LO(T)    = -10;
NORRAD.L(T)     = 2;
NORRAD.UP(T)    = 10;

SCALEMEANT.L    = 4.67384/40.4741;
SCALEMEANS.L    = 1.52056/10.1464;
* The values of SCALEMEANT and SCALEMEANS are fixed in the equation file.
* The two lines above are used to prevent a zero division.

T2M.LO(T)       = -5;
T2M.L(T)        = 0.5;
T2M.UP(T)       = 50;

PCT.LO(T)       = -10;
PCT.L(T)        = 0;
PCT.UP(T)       = 50;

PCTRES.LO(T,NPCR) = -10;
PCTRES.L(T,NPCR)  = 0;
PCTRES.UP(T,NPCR) = 50;

PCS.LO(T)       = -100;
PCS.L(T)        = 0;
PCS.UP(T)       = 500;

PCSRES.LO(T,NPCR) = -10;
PCSRES.L(T,NPCR)  = 0;
PCSRES.UP(T,NPCR) = 500;

SLRGLWT.LO(T)   = -10;
SLRGLWT.L(T)    = 0.1;
SLRGLWT.UP(T)   = 10;

SLRGLWC.LO(T)   = -10;
SLRGLWC.L(T)    = 0.1;
SLRGLWC.UP(T)   = 10;

SLRGR.LO(T)     = -2;
SLRGR.L(T)      = 0.1;
SLRGR.UP(T)     = 2;

SLRAN.LO(T)     = -2;
SLRAN.L(T)      = 0;
SLRAN.UP(T)     = 2;

SLRRT1750.LO(T) = -5;
SLRRT1750.L(T)  = 0.1;
SLRRT1750.UP(T) = 10;


* Following bounds were chosen according to the bound of SENS: [0.5, 15]

SQROOT.LO(T)    = 6;
SQROOT.L(T)     = 8;
SQROOT.UP(T)    = 20;

RPCTT.LO(T,'1') = 300;
RPCTT.L(T,'1')  = 450;
RPCTT.UP(T,'1') = 1200;

RPCTT.LO(T,'2') = 3;
RPCTT.L(T,'2')  = 15;
RPCTT.UP(T,'2') = 40;

RPCTA.LO(T,'1') = 0.028;
RPCTA.L(T,'1')  = 2.639;
RPCTA.UP(T,'1') = 18.973;

RPCTA.LO(T,'2') = 0.329;
RPCTA.L(T,'2')  = 6.642;
RPCTA.UP(T,'2') = 33.528;

SENSSCALE.LO(T) = 0.211;
SENSSCALE.L(T)  = 2.316;
SENSSCALE.UP(T) = 6.370;


* DOECLIM

KDIF.L          = 3155.8;
* The level for KDIF was calculated from CKD*KAPPA.L.
* The statement above is to avoid an error caused by zero division.

RADTTLL.LO(T)   = -10;
RADTTLL.L(T)    = 0.001;
RADTTLL.UP(T)   = 30;

RADTTLO.LO(T)   = -10;
RADTTLO.L(T)    = 0.001;
RADTTLO.UP(T)   = 30;

TSS.LO(T)       = -10;
TSS.L(T)        = 0.001;
TSS.UP(T)       = 50;

TAL.LO(T)       = -10;
TAL.L(T)        = 0.001;
TAL.UP(T)       = 50;

GMT.LO(T)       = -10;
GMT.L(T)        = 0.001;
GMT.UP(T)       = 50;

TSSMM.LO(I)     = -1;
TSSMM.L(I)      = 0.001;
TSSMM.UP(I)     = 3;

TSUM.LO(T)      = -50;
TSUM.L(T)       = 0.001;
TSUM.UP(T)      = 100;
