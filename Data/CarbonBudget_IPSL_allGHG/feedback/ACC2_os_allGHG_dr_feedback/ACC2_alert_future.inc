********** ACC2 version 4.2: ACC2_alert_future.inc *********************************************************************
********** File output to alert variables stuck at their bounds (future mode) ==========================================

* Note that, when upper or lower bounds are added or removed
* in "ACC2_bounds_future.inc", they should be reflected below.
FILE OUT_FUTURE_ALERT /out_future_alert.dat/;
OUT_FUTURE_ALERT.AP=0;
PUT OUT_FUTURE_ALERT;
PUT "Check the following variables are stuck at their bounds: "/;
PUT "(If no variables are listed below, there is no variable stuck at its bounds.) "/;
LOOP(T$(ord(T) GT 1),
PUT$((COC1ATM.L(T) GE COC1ATM.UP(T)) or (COC1ATM.L(T) LE COC1ATM.LO(T))) T.TL:4, @15, 'COC1ATM'/;
PUT$((COC1MIX.L(T) GE COC1MIX.UP(T)) or (COC1MIX.L(T) LE COC1MIX.LO(T))) T.TL:4, @15, 'COC1MIX'/;
PUT$((CONCO2INT.L(T) GE CONCO2INT.UP(T)) or (CONCO2INT.L(T) LE CONCO2INT.LO(T))) T.TL:4, @15, 'CONCO2INT'/;
PUT$((PCT.L(T) GE PCT.UP(T)) or (PCT.L(T) LE PCT.LO(T))) T.TL:4, @15, 'PCT'/;
PUT$((TSS.L(T) GE TSS.UP(T)) or (TSS.L(T) LE TSS.LO(T))) T.TL:4, @15, 'TSS'/;
PUT$((TAL.L(T) GE TAL.UP(T)) or (TAL.L(T) LE TAL.LO(T))) T.TL:4, @15, 'TAL'/;
PUT$((GMT.L(T) GE GMT.UP(T)) or (GMT.L(T) LE GMT.LO(T))) T.TL:4, @15, 'GMT'/;
PUT$((SLRGLWT.L(T) GE SLRGLWT.UP(T)) or (SLRGLWT.L(T) LE SLRGLWT.LO(T))) T.TL:4, @15, 'SLRGLWT'/;
PUT$((SLRGLWC.L(T) GE SLRGLWC.UP(T)) or (SLRGLWC.L(T) LE SLRGLWC.LO(T))) T.TL:4, @15, 'SLRGLWC'/;
PUT$((SLRGR.L(T) GE SLRGR.UP(T)) or (SLRGR.L(T) LE SLRGR.LO(T))) T.TL:4, @15, 'SLRGR'/;
PUT$((SLRAN.L(T) GE SLRAN.UP(T)) or (SLRAN.L(T) LE SLRAN.LO(T))) T.TL:4, @15, 'SLRAN'/;);

LOOP(T,
PUT$((H.L(T) GE H.UP(T)) or (H.L(T) LE H.LO(T))) T.TL:4, @15, 'H'/;
PUT$((COC1ALL.L(T) GE COC1ALL.UP(T)) or (COC1ALL.L(T) LE COC1ALL.LO(T))) T.TL:4, @15, 'COC1ALL'/;
PUT$((DIC.L(T) GE DIC.UP(T)) or (DIC.L(T) LE DIC.LO(T))) T.TL:4, @15, 'DIC'/;
PUT$((REVELLE.L(T) GE REVELLE.UP(T)) or (REVELLE.L(T) LE REVELLE.LO(T))) T.TL:4, @15, 'REVELLE'/;
PUT$((AK1S.L(T) GE AK1S.UP(T)) or (AK1S.L(T) LE AK1S.LO(T))) T.TL:4, @15, 'AK1S'/;
PUT$((AK1P.L(T) GE AK1P.UP(T)) or (AK1P.L(T) LE AK1P.LO(T))) T.TL:4, @15, 'AK1P'/;
PUT$((AK2S.L(T) GE AK2S.UP(T)) or (AK2S.L(T) LE AK2S.LO(T))) T.TL:4, @15, 'AK2S'/;
PUT$((AK2P.L(T) GE AK2P.UP(T)) or (AK2P.L(T) LE AK2P.LO(T))) T.TL:4, @15, 'AK2P'/;
PUT$((AKB.L(T) GE AKB.UP(T)) or (AKB.L(T) LE AKB.LO(T))) T.TL:4, @15, 'AKB'/;
PUT$((AKW.L(T) GE AKW.UP(T)) or (AKW.L(T) LE AKW.LO(T))) T.TL:4, @15, 'AKW'/;
PUT$((ALPHAS.L(T) GE ALPHAS.UP(T)) or (ALPHAS.L(T) LE ALPHAS.LO(T))) T.TL:4, @15, 'ALPHAS'/;
PUT$((TNPPFAC.L(T) GE TNPPFAC.UP(T)) or (TNPPFAC.L(T) LE TNPPFAC.LO(T))) T.TL:4, @15, 'TNPPFAC'/;
PUT$((TRESFAC.L(T) GE TRESFAC.UP(T)) or (TRESFAC.L(T) LE TRESFAC.LO(T))) T.TL:4, @15, 'TRESFAC'/;
PUT$((RADTTL.L(T) GE RADTTL.UP(T)) or (RADTTL.L(T) LE RADTTL.LO(T))) T.TL:4, @15, 'RADTTL'/;
PUT$((T2M.L(T) GE T2M.UP(T)) or (T2M.L(T) LE T2M.LO(T))) T.TL:4, @15, 'T2M'/;
PUT$((NORRAD.L(T) GE NORRAD.UP(T)) or (NORRAD.L(T) LE NORRAD.LO(T))) T.TL:4, @15, 'NORRAD'/;
PUT$((PCS.L(T) GE PCS.UP(T)) or (PCS.L(T) LE PCS.LO(T))) T.TL:4, @15, 'PCS'/;
PUT$((RADTTLL.L(T) GE RADTTLL.UP(T)) or (RADTTLL.L(T) LE RADTTLL.LO(T))) T.TL:4, @15, 'RADTTLL'/;
PUT$((RADTTLO.L(T) GE RADTTLO.UP(T)) or (RADTTLO.L(T) LE RADTTLO.LO(T))) T.TL:4, @15, 'RADTTLO'/;
PUT$((RADHAL.L(T) GE RADHAL.UP(T)) or (RADHAL.L(T) LE RADHAL.LO(T))) T.TL:4, @15, 'RADHAL'/;
PUT$((RADODS.L(T) GE RADODS.UP(T)) or (RADODS.L(T) LE RADODS.LO(T))) T.TL:4, @15, 'RADODS'/;
PUT$((RADSOL.L(T) GE RADSOL.UP(T)) or (RADSOL.L(T) LE RADSOL.LO(T))) T.TL:4, @15, 'RADSOL'/;
PUT$((TSUM.L(T) GE TSUM.UP(T)) or (TSUM.L(T) LE TSUM.LO(T))) T.TL:4, @15, 'TSUM'/;
PUT$((SQROOT.L(T) GE SQROOT.UP(T)) or (SQROOT.L(T) LE SQROOT.LO(T))) T.TL:4, @15, 'SQROOT'/;
PUT$((RPCTT.L(T,'1') GE RPCTT.UP(T,'1')) or (RPCTT.L(T,'1') LE RPCTT.LO(T,'1'))) T.TL:4, @15, 'RPCTT'/;
PUT$((RPCTT.L(T,'2') GE RPCTT.UP(T,'2')) or (RPCTT.L(T,'2') LE RPCTT.LO(T,'2'))) T.TL:4, @15, 'RPCTT'/;
PUT$((RPCTA.L(T,'1') GE RPCTA.UP(T,'1')) or (RPCTA.L(T,'1') LE RPCTA.LO(T,'1'))) T.TL:4, @15, 'RPCTA'/;
PUT$((RPCTA.L(T,'2') GE RPCTA.UP(T,'2')) or (RPCTA.L(T,'2') LE RPCTA.LO(T,'2'))) T.TL:4, @15, 'RPCTA'/;
PUT$((SENSSCALE.L(T) GE SENSSCALE.UP(T)) or (SENSSCALE.L(T) LE SENSSCALE.LO(T))) T.TL:4, @15, 'SENSSCALE'/;
PUT$((SLRRT1750.L(T) GE SLRRT1750.UP(T)) or (SLRRT1750.L(T) LE SLRRT1750.LO(T))) T.TL:4, @15, 'SLRRT1750'/;);

LOOP((T,NOC)$(ord(T) GT 1),
PUT$((COC.L(T,NOC) GE COC.UP(T,NOC)) or (COC.L(T,NOC) LE COC.LO(T,NOC))) T.TL:4, @15, NOC.TL:10, @30, 'COC'/;);

LOOP((T,FOB)$(ord(T) GT 1),
PUT$((ANT.L(T,FOB) GE ANT.UP(T,FOB)) or (ANT.L(T,FOB) LE ANT.LO(T,FOB))) T.TL:4, @15, FOB.TL:10, @30, 'ANT'/;);

LOOP((T,FOB)$((ord(T) GT 1) and (ord(FOB) NE 7) and (ord(FOB) NE 14) and (ord(FOB) NE 16) and (ord(FOB) NE 17)
and (ord(FOB) NE 24) and (ord(FOB) NE 36) and (ord(FOB) NE 37) and (ord(FOB) NE 38) and (ord(FOB) NE 40)),
PUT$((CON.L(T,FOB) GE CON.UP(T,FOB)) or (CON.L(T,FOB) LE CON.LO(T,FOB))) T.TL:4, @15, FOB.TL:10, @30, 'CON'/;);

LOOP((T,FOA),
PUT$((RAD.L(T,FOA) GE RAD.UP(T,FOA)) or (RAD.L(T,FOA) LE RAD.LO(T,FOA))) T.TL:4, @15, FOA.TL:10, @30, 'RAD'/;);

LOOP((T,NPCR)$(ord(T) GT 1),
PUT$((PCTRES.L(T,NPCR) GE PCTRES.UP(T,NPCR))
or (PCTRES.L(T,NPCR) LE PCTRES.LO(T,NPCR))) T.TL:4, @15, NPCR.TL:10, @30, 'PCTRES'/;
PUT$((PCSRES.L(T,NPCR) GE PCSRES.UP(T,NPCR))
or (PCSRES.L(T,NPCR) LE PCSRES.LO(T,NPCR))) T.TL:4, @15, NPCR.TL:10, @30, 'PCSRES'/;);

LOOP(I$(ord(I) GT 251),
PUT$((TSSMM.L(I) GE TSSMM.UP(I)) or (TSSMM.L(I) LE TSSMM.LO(I))) I.TL:4, @15, 'TSSMM'/;);
