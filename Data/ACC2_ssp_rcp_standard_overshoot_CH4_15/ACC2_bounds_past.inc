********** ACC2 version 4.2: ACC2_bounds_past.inc **********************************************************************
********** Levels and bounds for valuables (past mode) =================================================================

* Note that, when upper or lower bounds are added or removed below,
* they should be reflected in "OUT_PAST_ALERT" section in "ACC2_main_past.gms".

H.LO(T)         = 2;
H.L(T)          = 5;
H.UP(T)         = 20;

COC.LO(T,'1')   = -50;
COC.L(T,'1')    = 100;
COC.UP(T,'1')   = 500;

COC.LO(T,'2')   = -20;
COC.L(T,'2')    = 20;
COC.UP(T,'2')   = 100;

COC.LO(T,'3')   = -10;
COC.L(T,'3')    = 20;
COC.UP(T,'3')   = 100;

COC.LO(T,'4')   = -10;
COC.L(T,'4')    = 5;
COC.UP(T,'4')   = 100;

COC1ATM.LO(T)   = -50;
COC1ATM.L(T)    = 0.001;
COC1ATM.UP(T)   = 500;

COC1MIX.LO(T)   = -10;
COC1MIX.L(T)    = 0.001;
COC1MIX.UP(T)   = 100;

COC1ATMBAC.LO   = 500;
COC1ATMBAC.L    = 594;
COC1ATMBAC.UP   = 700;

COC1MIXBAC.LO   = 400;
COC1MIXBAC.L    = 504.025863;
COC1MIXBAC.UP   = 600;

COC1ALL.LO(T)   = 500;
COC1ALL.L(T)    = 1108;
COC1ALL.UP(T)   = 2000;

DIC.LO(T)       = 1.000E-4;
DIC.L(T)        = 2.000E-3;
DIC.UP(T)       = 1.000E-1;

REVELLE.LO(T)   = 1;
REVELLE.L(T)    = 10;
REVELLE.UP(T)   = 100;

PH.LO(T)        = 6.0;
PH.L(T)         = 8.15;
PH.UP(T)        = 10;

CO2AQ.LO(T)     = 1.0E-7;
CO2AQ.L(T)      = 1.0E-5;
CO2AQ.UP(T)     = 1.0E-3;

HCO3.LO(T)      = 1.0E-4;
HCO3.L(T)       = 0.002;
HCO3.UP(T)      = 1.0E-2;

CO32.LO(T)      = 1.0E-5;
CO32.L(T)       = 2.3E-4;
CO32.UP(T)      = 1.0E-3;

FRCO2AQ.LO(T)   = 1.0E-4;
FRCO2AQ.L(T)    = 0.005;
FRCO2AQ.UP(T)   = 0.05;

FRHCO3.LO(T)    = 0.01;
FRHCO3.L(T)     = 0.88;
FRHCO3.UP(T)    = 1.0;

FRCO32.LO(T)    = 0.01;
FRCO32.L(T)     = 0.1;
FRCO32.UP(T)    = 0.5;

CO2FRT.LO(T,'1') = -50;
CO2FRT.L(T,'1')  = -5.0;
CO2FRT.UP(T,'1') = 10;

CO2FRT.LO(T,'2') = -10;
CO2FRT.L(T,'2')  = 10;
CO2FRT.UP(T,'2') = 100;

CO2FRT.LO(T,'3') = -2.0;
CO2FRT.L(T,'3')  = 1.0;
CO2FRT.UP(T,'3') = 10;

CO2FRT.LO(T,'4') = -2.0;
CO2FRT.L(T,'4')  = -1.0;
CO2FRT.UP(T,'4') = 10;

HTRRES.LO(T,'1') = -500;
HTRRES.L(T,'1')  = -100;
HTRRES.UP(T,'1') = 0;

HTRRES.LO(T,'2') = 10;
HTRRES.L(T,'2')  = 130;
HTRRES.UP(T,'2') = 500;

HTRRES.LO(T,'3') = 5;
HTRRES.L(T,'3')  = 17;
HTRRES.UP(T,'3') = 50;

HTRRES.LO(T,'4') = 5;
HTRRES.L(T,'4')  = 18;
HTRRES.UP(T,'4') = 50;

FLXCOC1.LO(T)   =  -7;
FLXCOC1.L(T)    =  1;
FLXCOC1.UP(T)   =  10;

FLXCOC2.LO(T)   =  -3;
FLXCOC2.L(T)    =  0.5;
FLXCOC2.UP(T)   =  5;

FLXCOC3.LO(T)   =  -1;
FLXCOC3.L(T)    =  0.3;
FLXCOC3.UP(T)   =  5;

FLXCOC4.LO(T)   =  -1;
FLXCOC4.L(T)    =  0.05;
FLXCOC4.UP(T)   =  1;

FLXCTR.LO(T,'1') = -5;
FLXCTR.L(T,'1')  = -0.3;
FLXCTR.UP(T,'1') = 3;

FLXCTR.LO(T,'2') = -5;
FLXCTR.L(T,'2')  = 0.5;
FLXCTR.UP(T,'2') = 10;

FLXCTR.LO(T,'3') = -1;
FLXCTR.L(T,'3')  = 0.3;
FLXCTR.UP(T,'3') = 3;

FLXCTR.LO(T,'4') = -1;
FLXCTR.L(T,'4')  = 0.8;
FLXCTR.UP(T,'4') = 5;

CO2OCNUPANT.LO(T) = -5;
CO2OCNUPANT.L(T)  = 1.0;
CO2OCNUPANT.UP(T) = 5;

CO2LNDUPANT.LO(T) = -5;
CO2LNDUPANT.L(T)  = 1.0;
CO2LNDUPANT.UP(T) = 5;

CO2OCNUPTTL.LO(T) = -5;
CO2OCNUPTTL.L(T)  = 1.0;
CO2OCNUPTTL.UP(T) = 5;

CO2LNDUPTTL.LO(T) = -5;
CO2LNDUPTTL.L(T)  = 1.0;
CO2LNDUPTTL.UP(T) = 5;

SUMCOC.LO(T)      = -20;
SUMCOC.L(T)       = 10;
SUMCOC.UP(T)      = 200;

SUMCTR.LO(T)      = 1000;
SUMCTR.L(T)       = 2300;
SUMCTR.UP(T)      = 3000;

CTR.LO(T,'1')   = -300;
CTR.L(T,'1')    = -220;
CTR.UP(T,'1')   = -100;

CTR.LO(T,'2')   = 200;
CTR.L(T,'2')    = 360;
CTR.UP(T,'2')   = 1000;

CTR.LO(T,'3')   = 200;
CTR.L(T,'3')    = 350;
CTR.UP(T,'3')   = 1000;

CTR.LO(T,'4')   = 1000;
CTR.L(T,'4')    = 1800;
CTR.UP(T,'4')   = 5000;

CTR1750.LO('1') = -300;
CTR1750.L('1')  = -205;
CTR1750.UP('1') = -100;

CTR1750.LO('2') = 200;
CTR1750.L('2')  = 343;
CTR1750.UP('2') = 1000;

CTR1750.LO('3') = 200;
CTR1750.L('3')  = 321;
CTR1750.UP('3') = 1000;

CTR1750.LO('4') = 1000;
CTR1750.L('4')  = 1759;
CTR1750.UP('4') = 5000;

NPPBAC.LO(T,'1')= -200;
NPPBAC.L(T,'1') = -94.053;
NPPBAC.UP(T,'1')= 0;

NPPBAC.LO(T,'2')= 0;
NPPBAC.L(T,'2') = 120.379;
NPPBAC.UP(T,'2')= 500;

NPPBAC.LO(T,'3')= 0;
NPPBAC.L(T,'3') = 16.097;
NPPBAC.UP(T,'3')= 100;

NPPBAC.LO(T,'4')= 0;
NPPBAC.L(T,'4') = 17.594;
NPPBAC.UP(T,'4')= 100;

CTRREL.LO(T,'1')= -100;
CTRREL.L(T,'1') = -15;
CTRREL.UP(T,'1')= 50;

CTRREL.LO(T,'2')= -10;
CTRREL.L(T,'2') = 25;
CTRREL.UP(T,'2')= 100;

CTRREL.LO(T,'3')= -10;
CTRREL.L(T,'3') = 20;
CTRREL.UP(T,'3')= 100;

CTRREL.LO(T,'4')= -10;
CTRREL.L(T,'4') = 50;
CTRREL.UP(T,'4')= 500;

AK1S.LO(T)      = 5.000E2;
AK1S.L(T)       = 1.291E3;
AK1S.UP(T)      = 2.000E3;

AK1P.LO(T)      = 3.500E2;
AK1P.L(T)       = 4.160E2;
AK1P.UP(T)      = 5.000E2;

AK2S.LO(T)      = 3.000E-1;
AK2S.L(T)       = 8.800E-1;
AK2S.UP(T)      = 1.500E0;

AK2P.LO(T)      = 1.000E-2;
AK2P.L(T)       = 4.140E-2;
AK2P.UP(T)      = 8.000E-2;

AKB.LO(T)       = 1.000E0;
AKB.L(T)        = 2.256E0;
AKB.UP(T)       = 5.000E0;

AKW.LO(T)       = 1.000E4;
AKW.L(T)        = 3.600E4;
AKW.UP(T)       = 9.000E4;

ALPHAS.LO(T)    = 1.000E-2;
ALPHAS.L(T)     = 2.900E-2;
ALPHAS.UP(T)    = 5.000E-2;

*TMIX.LO(T)      = 15;
TMIX.L(T)       = 19.59;
*TMIX.UP(T)      = 25;
* Adding bounds for TMIX lead to infeasibility for GAMSCHK.

TNPPFAC.LO(T)   = 0.001;
TNPPFAC.L(T)    = 1;
TNPPFAC.UP(T)   = 3;

TRESFAC.LO(T)   = 0.001;
TRESFAC.L(T)    = 1;
TRESFAC.UP(T)   = 3;

CO2FF.LO(T)     = -5;
CO2FF.L(T)      = 0.001;
CO2FF.UP(T)     = 10;

CO2LU.LO(T)     = -5;
CO2LU.L(T)      = 0.001;
CO2LU.UP(T)     = 10;

CONCO2INT.LO(T) = 150;
CONCO2INT.L(T)  = 300;
CONCO2INT.UP(T) = 1000;

CON.LO(T,'CO2') = 250;
CON.L(T,'CO2')  = 300;
CON.UP(T,'CO2') = 500;

CON.LO(T,'CH4') = 300;
CON.L(T,'CH4')  = 750;
CON.UP(T,'CH4') = 2000;

CON.LO(T,'N2O') = 200;
CON.L(T,'N2O')  = 250;
CON.UP(T,'N2O') = 400;

CON.LO(T,'rOH') = 0;
CON.L(T,'rOH')  = 1;
CON.UP(T,'rOH') = 10;

ANT.LO(T,'CO2')  = -15;
ANT.L(T,'CO2')   = 0.1;
ANT.UP(T,'CO2')  = 20;

ANT.LO(T,'CH4')    = -100;
ANT.L(T,'CH4')     = 50;
ANT.UP(T,'CH4')    = 1000;

ANTPRI.LO(T,'CH4')    = -100;
ANTPRI.L(T,'CH4')     = 50;
ANTPRI.UP(T,'CH4')    = 1000;

ANT.LO(T,'N2O')    = -10;
ANT.L(T,'N2O')     = 1;
ANT.UP(T,'N2O')    = 30;

ANTPRI.LO(T,'N2O')    = -10;
ANTPRI.L(T,'N2O')     = 1;
ANTPRI.UP(T,'N2O')    = 30;

ANT.LO(T,'NOx')    = -10;
ANT.L(T,'NOx')     = 5;
ANT.UP(T,'NOx')    = 100;

ANT.LO(T,'CO')    = -100;
ANT.L(T,'CO')     = 50;
ANT.UP(T,'CO')    = 1000;

ANT.LO(T,'VOC')    = -200;
ANT.L(T,'VOC')     = 20;
ANT.UP(T,'VOC')    = 1000;

ANT.LO(T,'SO2')    = -40;
ANT.L(T,'SO2')     = 0.5;
ANT.UP(T,'SO2')    = 100;

DCONCO2ANT.LO(T) = -10;
DCONCO2ANT.L(T)  = 0;
DCONCO2ANT.UP(T) = 10;

DCONCO2NAT.LO(T) = -10;
DCONCO2NAT.L(T)  = 0;
DCONCO2NAT.UP(T) = 10;

DCONCO2TTL.LO(T) = -10;
DCONCO2TTL.L(T)  = 0;
DCONCO2TTL.UP(T) = 10;

DCONCO2ENSO.LO(T)= -3;
DCONCO2ENSO.L(T) = 0;
DCONCO2ENSO.UP(T)= 3;

AVENINO3.LO      = -3.0;
AVENINO3.L       = 0.01;
AVENINO3.UP      = 3.0;

AVEDCONCO2NAT.LO = -5.0;
AVEDCONCO2NAT.L  = 0.1;
AVEDCONCO2NAT.UP = 5.0;

SX2.LO           = 0.0;
SX2.L            = 1.0;
SX2.UP           = 100;

SXY.LO           = 0.0;
SXY.L            = 1.0;
SXY.UP           = 100;

YINTRCO2ENSO.LO = -5.0;
YINTRCO2ENSO.L  = 0.0;
YINTRCO2ENSO.UP = 5.0;

SLOPECO2ENSO.LO = -2.0;
SLOPECO2ENSO.L  = 1.0;
SLOPECO2ENSO.UP = 4.0;

RAD.LO(T,FOA)   = -15;
RAD.L(T,FOA)    = 0.001;
RAD.UP(T,FOA)   = 15;

RADHAL.LO(T)    = -1.0;
RADHAL.L(T)     = 0.03;
RADHAL.UP(T)    = 1.0;

RADODS.LO(T)    = -1.0;
RADODS.L(T)     = 0;
RADODS.UP(T)    = 1.0;

RADMIS.LO(T)    = -15;
RADMIS.L(T)     = 0.001;
RADMIS.UP(T)    = 15;

RADTTLPRI.LO(T) = -15;
RADTTLPRI.L(T)  = 0.001;
RADTTLPRI.UP(T) = 15;

RADTTL.LO(T)    = -15;
RADTTL.L(T)     = 0.001;
RADTTL.UP(T)    = 15;

SCALEMEANT.L    = 4.67384/40.4741;
* To prevent zero division errors upon execution.

T2M.LO(T)       = -3;
T2M.L(T)        = 0.001;
T2M.UP(T)       = 7;

T2M1990.LO      = 0;
T2M1990.L       = 0.7;
T2M1990.UP      = 1.5;

T2MRT1990.LO(T) = -7;
T2MRT1990.L(T)  = -0.3;
T2MRT1990.UP(T) = 7;

NORRAD.LO(T)    = -10;
NORRAD.L(T)     = 0;
NORRAD.UP(T)    = 10;

PCT.LO(T)       = -5;
PCT.L(T)        = 0.001;
PCT.UP(T)       = 5;

PCTRES.LO(T,NPCR) = -100;
PCTRES.L(T,NPCR)  = 0;
PCTRES.UP(T,NPCR) = 100;

PCS.LO(T)       = -1;
PCS.L(T)        = 0;
PCS.UP(T)       = 1;

PCSRES.LO(T,NPCR) = -1;
PCSRES.L(T,NPCR)  = 0;
PCSRES.UP(T,NPCR) = 1;

DT2MENSO.LO(T)  = -1;
DT2MENSO.L(T)   = 0;
DT2MENSO.UP(T)  = 1;


* DOECLIM

KSL.LO          = 0.001;
* To prevent zero division errors upon execution.

RADANTL.LO(T)   = -15;
RADANTL.L(T)    = 0.2;
RADANTL.UP(T)   = 15;

RADANTO.LO(T)   = -15;
RADANTO.L(T)    = 0.2;
RADANTO.UP(T)   = 15;

RADSOLL.LO(T)   = -2;
RADSOLL.L(T)    = 0;
RADSOLL.UP(T)   = 2;

RADSOLO.LO(T)   = -2;
RADSOLO.L(T)    = 0;
RADSOLO.UP(T)   = 2;

RADNUML.LO(T)   = -50;
RADNUML.L(T)    = 0;
RADNUML.UP(T)   = 50;

RADNUMO.LO(T)   = -50;
RADNUMO.L(T)    = 0;
RADNUMO.UP(T)   = 50;

RADTTLL.LO(T)   = -15;
RADTTLL.L(T)    = 0;
RADTTLL.UP(T)   = 15;

RADTTLO.LO(T)   = -15;
RADTTLO.L(T)    = 0;
RADTTLO.UP(T)   = 15;

TSS.LO(T)       = -3;
TSS.L(T)        = 0.001;
TSS.UP(T)       = 5;

TAL.LO(T)       = -3;
TAL.L(T)        = 0.001;
TAL.UP(T)       = 7;

GMT.LO(T)       = -3;
GMT.L(T)        = 0.001;
GMT.UP(T)       = 7;

KER.LO(I)       = -5;
KER.L(I)        = 0.1;
KER.UP(I)       = 5;

KER1.LO(I)      = -5;
KER1.L(I)       = 0.1;
KER1.UP(I)      = 5;

TSSMM.LO(I)     = -7;
TSSMM.L(I)      = 0.001;
TSSMM.UP(I)     = 7;

TSUM.LO(T)      = -3;
TSUM.L(T)       = 0.001;
TSUM.UP(T)      = 3;


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


* Uncertain parameters

BETA.LO(T)      = 0.0;
BETA.L(T)       = 0.4;
BETA.UP(T)      = 1.5;

TMIX1750.LO     = 10;
TMIX1750.L      = 19.6;
TMIX1750.UP     = 40;

DMPTMIX.LO      = -1;
DMPTMIX.L       = 0.5;
DMPTMIX.UP      = 3;

Q10EFF.LO(T)    = 0.001;
Q10EFF.L(T)     = 2.25;
Q10EFF.UP(T)    = 5.0;

Q10BAS.LO       = 0.001;
Q10BAS.L        = 2.25;
Q10BAS.UP       = 5.0;

Q10TMP.LO       = 0.001;
Q10TMP.L        = 0.046;
Q10TMP.UP       = 1.0;

CO2OCNUPPRE.LO(T) = -5;
CO2OCNUPPRE.L(T)  = 0.001;
CO2OCNUPPRE.UP(T) = 5;

CO2LNDUPPRE.LO(T) = -5;
CO2LNDUPPRE.L(T)  = 0.001;
CO2LNDUPPRE.UP(T) = 5;

NATCH4.LO       = 100;
NATCH4.L        = 210;
NATCH4.UP       = 500;

NATN2O.LO       = 5;
NATN2O.L        = 10.2;
NATN2O.UP       = 30;

TAUCH4OH.LO     = 5;
TAUCH4OH.L      = 9.6;
TAUCH4OH.UP     = 15;

TAUN2O.LO       = 60;
TAUN2O.L        = 110;
TAUN2O.UP       = 160;

RAD2XCO2.LO     = 2.5;
RAD2XCO2.L      = 3.7;
RAD2XCO2.UP     = 5.0;

SCLRAD.LO       = 0.01;
SCLRAD.L        = 1.0;
SCLRAD.UP       = 5.0;

SENS.LO(T)      = 0.1;
SENS.L(T)       = 3.5;
SENS.UP(T)      = 13;

* SENS.FX(T)      = 3.26;
* DOECLIM original setting

KAPPA.LO        = 0;
KAPPA.L         = 0.55;
KAPPA.UP        = 5.0;
* To constrain KDIF instead of KAPPA, look at below.

KDIF.LO         = 0;
KDIF.L          = 1735.69;
KDIF.UP         = 15179;
* The level for KDIF was calculated from CKD*KAPPA.L.
* To prevent zero division errors upon execution.

RADMISAR1PRO.LO    = -5.0;
RADMISAR1PRO.L     = 0.5;
RADMISAR1PRO.UP    = 5.0;

T2MAR1PRO.LO    = -2.0;
T2MAR1PRO.L     = 0.45;
T2MAR1PRO.UP    = 2.0;


* To check the regression line

RSQ1.LO         = -1.0;
RSQ1.L          =  0.0;
RSQ1.UP         =  1.0;

RSQ2.LO         = -1.0;
RSQ2.L          =  0.0;
RSQ2.UP         =  1.0;


* Associated with cost function

CO2FFRES.LO(T)  = -10;
CO2FFRES.L(T)   = 0;
CO2FFRES.UP(T)  = 10;

CO2FFSGMMOD.L(T)= 0.1;

CO2FFRESSUM.LO  = 0;
CO2FFRESSUM.L   = 0;
CO2FFRESSUM.UP  = 5000000;

CO2LURES.LO(T)  = -10;
CO2LURES.L(T)   = 0;
CO2LURES.UP(T)  = 10;

CO2LUSGMMOD.L(T)= 0.1;

CO2LURESSUM.LO  = 0;
CO2LURESSUM.L   = 0;
CO2LURESSUM.UP  = 5000000;

ANTRES.LO(T,'CH4') = -100;
ANTRES.L(T,'CH4')  = 0;
ANTRES.UP(T,'CH4') = 100;

ANTRES.LO(T,'N2O') = -6;
ANTRES.L(T,'N2O')  = 0;
ANTRES.UP(T,'N2O') = 6;

ANTSGMMOD.L(T,FOG) = 0.001;

ANTRESSUMCH4.LO = 0;
ANTRESSUMCH4.L  = 0;
ANTRESSUMCH4.UP = 100;

ANTRESSUMN2O.LO = 0;
ANTRESSUMN2O.L  = 0;
ANTRESSUMN2O.UP = 100;

RADMISRES.LO(T)    = -10;
RADMISRES.L(T)     = 0;
RADMISRES.UP(T)    = 10;

RADMISSGMMOD.L(T)  = 0.1;

RADMISRESSUM.LO    = 0;
RADMISRESSUM.L     = 0;
RADMISRESSUM.UP    = 100000;


TMIX1750RES.LO     = -5;
TMIX1750RES.L      = 0;
TMIX1750RES.UP     = 5;

TMIX1750RESSUM.LO  = 0;
TMIX1750RESSUM.L   = 0;
TMIX1750RESSUM.UP  = 1000;

DMPTMIXRES.LO      = -2;
DMPTMIXRES.L       = 0;
DMPTMIXRES.UP      = 2;

DMPTMIXRESSUM.LO   = 0;
DMPTMIXRESSUM.L    = 0;
DMPTMIXRESSUM.UP   = 64;

BETARES.LO(T)      = -0.9;
BETARES.L(T)       = 0;
BETARES.UP(T)      = 2.0;

BETARESSUM.LO      = 0;
BETARESSUM.L       = 0;
BETARESSUM.UP      = 100;

Q10BASRES.LO       = -5;
Q10BASRES.L        = 0;
Q10BASRES.UP       = 5;

Q10TMPRES.LO       = -1;
Q10TMPRES.L        = 0;
Q10TMPRES.UP       = 1;

Q10BASRESSUM.LO    = 0;
Q10BASRESSUM.L     = 0;
Q10BASRESSUM.UP    = 128;

Q10TMPRESSUM.LO    = 0;
Q10TMPRESSUM.L     = 0;
Q10TMPRESSUM.UP    = 128;

TNPPSENSRES.LO     = -1;
TNPPSENSRES.L      = 0;
TNPPSENSRES.UP     = 1;

TNPPSENSRESSUM.LO  = 0;
TNPPSENSRESSUM.L   = 0;
TNPPSENSRESSUM.UP  = 16;

CO2OCNUPPRERES.LO(T) = -3;
CO2OCNUPPRERES.L(T)  = 0;
CO2OCNUPPRERES.UP(T) = 3;

CO2OCNUPPRERESSUM.LO = 0;
CO2OCNUPPRERESSUM.L  = 0;
CO2OCNUPPRERESSUM.UP = 10000;

CO2LNDUPPRERES.LO(T) = -3;
CO2LNDUPPRERES.L(T)  = 0;
CO2LNDUPPRERES.UP(T) = 3;

CO2LNDUPPRERESSUM.LO = 0;
CO2LNDUPPRERESSUM.L  = 0;
CO2LNDUPPRERESSUM.UP = 10000;

NATCH4RES.LO    = -210;
NATCH4RES.L     = 0;
NATCH4RES.UP    = 210;

NATCH4RESSUM.LO = 0;
NATCH4RESSUM.L  = 0;
NATCH4RESSUM.UP = 10000;

NATN2ORES.LO    = -10;
NATN2ORES.L     = 0;
NATN2ORES.UP    = 10;

NATN2ORESSUM.LO = 0;
NATN2ORESSUM.L  = 0;
NATN2ORESSUM.UP = 10000;

TAUCH4OHRES.LO    = -5;
TAUCH4OHRES.L     = 0;
TAUCH4OHRES.UP    = 5;

TAUCH4OHRESSUM.LO = 0;
TAUCH4OHRESSUM.L  = 0;
TAUCH4OHRESSUM.UP = 1000;

TAUN2ORES.LO    = -50;
TAUN2ORES.L     = 0;
TAUN2ORES.UP    = 50;

TAUN2ORESSUM.LO = 0;
TAUN2ORESSUM.L  = 0;
TAUN2ORESSUM.UP = 1000;

RAD2XCO2RES.LO     = -2;
RAD2XCO2RES.L      = 0;
RAD2XCO2RES.UP     = 2;

RAD2XCO2RESSUM.LO  = 0;
RAD2XCO2RESSUM.L   = 0;
RAD2XCO2RESSUM.UP  = 100;

SCLRADRES.LO       = -1;
SCLRADRES.L        = 0;
SCLRADRES.UP       = 4;

SCLRADRESSUM.LO    = 0;
SCLRADRESSUM.L     = 0;
SCLRADRESSUM.UP    = 100;

SENSRES.LO(T)      = -15;
SENSRES.L(T)       = 0;
SENSRES.UP(T)      = 15;

SENSRESSUM.LO      = 0;
SENSRESSUM.L       = 0;
SENSRESSUM.UP      = 10000;

KAPPARES.LO        = -2.5;
KAPPARES.L         = 0;
KAPPARES.UP        = 6;

KAPPARESSUM.LO     = 0;
KAPPARESSUM.L      = 0;
KAPPARESSUM.UP     = 100;

RADMISAR1PRORES.LO = -5;
RADMISAR1PRORES.L  = 0;
RADMISAR1PRORES.UP = 5;

RADMISAR1PRORESSUM.LO  = 0;
RADMISAR1PRORESSUM.L   = 0;
RADMISAR1PRORESSUM.UP  = 100;

T2MAR1PRORES.LO = -2;
T2MAR1PRORES.L  = 0;
T2MAR1PRORES.UP = 2;

T2MAR1PRORESSUM.LO  = 0;
T2MAR1PRORESSUM.L   = 0;
T2MAR1PRORESSUM.UP  = 100;


CO2OCNUPANTRES.LO(T) = -5.0;
CO2OCNUPANTRES.L(T)  = 0;
CO2OCNUPANTRES.UP(T) = 5.0;

CO2OCNUPANTRESSUM.LO = 0;
CO2OCNUPANTRESSUM.L  = 0;
CO2OCNUPANTRESSUM.UP = 100000;

CO2LNDUPANTRES.LO(T) = -5.0;
CO2LNDUPANTRES.L(T)  = 0;
CO2LNDUPANTRES.UP(T) = 5.0;

CO2LNDUPANTRESSUM.LO = 0;
CO2LNDUPANTRESSUM.L  = 0;
CO2LNDUPANTRESSUM.UP = 100000;

CONCO2RES.LO(T) = -20;
CONCO2RES.L(T)  = 0;
CONCO2RES.UP(T) = 20;

CONCO2SGMMOD.L(T) = 0.1;

CONCO2RESSUM.LO = 0;
CONCO2RESSUM.L  = 0;
CONCO2RESSUM.UP = 100000;

CONCH4RES.LO(T) = -300;
CONCH4RES.L(T)  = 0;
CONCH4RES.UP(T) = 300;

CONCH4SGMMOD.L(T) = 0.1;

CONCH4RESSUM.LO = 0;
CONCH4RESSUM.L  = 0;
CONCH4RESSUM.UP = 100000;

CONN2ORES.LO(T) = -50;
CONN2ORES.L(T)  = 0;
CONN2ORES.UP(T) = 50;

CONN2OSGMMOD.L(T) = 0.1;

CONN2ORESSUM.LO = 0;
CONN2ORESSUM.L  = 0;
CONN2ORESSUM.UP = 100000;

T2MRES.LO(T)    = -7.0;
T2MRES.L(T)     = 0;
T2MRES.UP(T)    = 7.0;

T2MSGMMOD.L(T)  = 0.1;

T2MRESSUM.LO    = 0;
T2MRESSUM.L     = 0;
T2MRESSUM.UP    = 100000;
