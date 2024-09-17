********** ACC2 version 4.2: ACC2_equations_common.inc *****************************************************************
********** Equation declarations and definitions (past and future modes) ===============================================

****** Equation declarations ===========================================================================================

EQUATIONS

COC1ALLQ(T)

COC1BALQ(T)

COC1ATMQ(T)

COC1MIXQ(T)

ALKLQ(T)


BUFFERQ(T)

REVELLEQ(T)


PHQ(T)

DICQ(T)

CO2AQQ(T)

HCO3Q(T)

CO32Q(T)

FRCO2AQQ(T)

FRHCO3Q(T)

FRCO32Q(T)


NPPBACQ(T,NTR)

CO2FRTQ(T,NTR)

HTRRESQ(T,NTR)


TMIXQ1(T)

TMIXQ2(T)

TMIXQ3(T)

AK1SQ(T)

AK1PQ(T)

AK1AQ(T)

AK1BQ(T)

AK1CQ(T)

AK2SQ(T)

AK2PQ(T)

AK2AQ(T)

AK2BQ(T)

AK2CQ(T)

AKBQ(T)

AKWQ(T)

ALPHASQ(T)


CO2OCNUPANTQ1(T)

CO2OCNUPANTQ2(T)

CO2OCNUPANTQ3(T)

CO2OCNUPANTQ4(T)

CO2OCNUPTTLQ1(T)

CO2OCNUPTTLQ2(T)

CO2OCNUPTTLQ3(T)

CO2OCNUPTTLQ4(T)

CO2OCNUPTTLQ5(T)

CO2OCNUPTTLQ6(T)


CO2LNDUPANTQ1(T)

CO2LNDUPANTQ2(T)

CO2LNDUPANTQ3(T)

CO2LNDUPANTQ4(T)

CO2LNDUPTTLQ1(T)

CO2LNDUPTTLQ2(T)

CO2LNDUPTTLQ3(T)

CO2LNDUPTTLQ4(T)

CO2LNDUPTTLQ5(T)

CO2LNDUPTTLQ6(T)


CO2OCNUPPRECNSTQ(T)

CO2LNDUPPRECNSTQ(T)


COCEQ1(T)

COCEQ2(T)

COCEQ3(T)

COCEQ4(T)

CTREQ(T,NTR)

FLXCOC1Q(T)

FLXCOC2Q(T)

FLXCOC3Q(T)

FLXCOC4Q(T)

FLXCTRQ(T,NTR)

SUMCOCQ(T)

SUMCTRQ(T)

CTRRELQ(T,NTR)


CUMCO2Q1(T)

CUMCO2Q2(T)

CUMCO2Q3(T)

CUMCO2Q4(T)


TNPPFACQ1(T)

TNPPFACQ2(T)

TNPPFACQ3(T)

TNPPFACQ4(T)

TRESFACQ1(T)

TRESFACQ2(T)

TRESFACQ3(T)

TRESFACQ4(T)

Q10EFFQ1(T)

Q10EFFQ2(T)

Q10EFFQ3(T)

Q10EFFQ4(T)

Q10EFFQ5(T)

Q10EFFQ6(T)


CO2CONC1(T)

CO2CONC2(T)

CO2CONC3(T)

CO2CONC4(T)

CO2CONC5(T)

DCONCO2ENSOQ(T)

CONCH4Q(T)

CONN2OQ(T)

CONOHQ(T)


RADCO2Q(T)

RADCH4Q(T)

RADN2OQ(T)

RADHALTQ(T)

RADODSQ(T)

RADSO2Q(T)

RADAERINDQ(T)

RADCRBAERQ(T)

RADH2OSTRQ(T)


NORRADQ(T)

SCALEMEANTQ

SCALEMEANSQ

SQROOTQ(T)

RPCTT1Q(T)

RPCTT2Q(T)

RPCTA1Q(T)

RPCTA2Q(T)

SENSSCALEQ(T)


PCTRESQ(T,NPCR)

PCTQ(T)

T2MQ1(T)

T2MQ2(T)

T2MQ3(T)

T2MQ4(T)

T2MQ5(T)

T2MQ6(T)

T2MRT1990Q(T)


PCSRESQ(T,NPCR)

PCSQ(T)

PCSQ1990(T)


SLRGLQ1(T)            Sea level rise due to glaciers and small ice caps without glacier area contraction (for 1990)

SLRGLQ2P(T)           Sea level rise due to glaciers and small ice caps without glacier area contraction
*                     (for 1995 onward)

SLRGLQ2F(T)           Sea level rise due to glaciers and small ice caps without glacier area contraction

SLRGLQ3P(T)           Sea level rise due to glaciers and small ice caps with glacier area contraction (for 1995 onward)

SLRGLQ3F(T)           Sea level rise due to glaciers and small ice caps with glacier area contraction

SLRGLQ1990(T)         Sea level rise due to glaciers and small ice caps with glacier area contraction relative to 1990

SLRGRQ1(T)            Greenland ice sheet contribution to sea level rise (for 1990)

SLRGRQ2P(T)           Greenland ice sheet contribution to sea level rise (for 1995 onward)

SLRGRQ2F(T)           Greenland ice sheet contribution to sea level rise

SLRGRQ1990(T)         Greenland ice sheet contribution to sea level rise relative to 1990

SLRANQ1(T)            Antarctic ice sheet contribution to sea level rise (for 1990)

SLRANQ2P(T)           Antarctic ice sheet contribution to sea level rise (for 1995 onward)

SLRANQ2F(T)           Antarctic ice sheet contribution to sea level rise

SLRANQ1990(T)         Antarctic ice sheet contribution to sea level rise relative to 1990

SLRQP(T)              Total sea level rise

SLRQF(T)              Total sea level rise

SLRQ1990(T)           Total sea level rise relative to 1990


BETACONT(T)           Fix BETA as time-invariant

SENSCONT(T)           Fix SENS as time-invariant


KP1Q                  Land ocean heat transfer: Offset

KP2Q                  Land ocean heat transfer: Slope

LP1Q                  Land feedback strength: Offset

LP2Q                  Land feedback strength: Slope

OP1Q                  Ocean feedback strength: Offset

OP2Q                  Ocean feedback strength: Slope


RADSOLLQ(T)           Solar forcing over the land for past mode

RADSOLOQ(T)           Solar forcing over the ocean for past mode

RADTTLLQ(T)           Total radiative forcing over the land

RADTTLOQ(T)           Total radiative forcing over the ocean

RADNUMLQ(t)           Numerical implementation of forcing

RADNUMOQ(t)           Numerical implementation of forcing


CLIMP0

CLIMP1

CLIMP2

KHEAT


CONVOL10(I)

CONVOL1(I)

*CONVOL20(I)
*
*CONVOL2(I)
*
*CONVOL30(I)
*
*CONVOL3(I)

CONVOLS(I)

CLIMHLPS0(T)

CLIMHLPS(T)

CLIMHLPF0(T)

CLIMHLPF(T)


CLIMA11

CLIMA12

CLIMA21

CLIMA22

CLIMB11

CLIMB12

CLIMB21

CLIMB22

CLIMC11

CLIMC12

CLIMC21

CLIMC22


CLIMTAL(T)

CLIMTSS(T)

CLIMGMT(T)


FOQ(T)

HCIOQ(T)

HCMLQ(T)

HCALQ(T)

HCTLQ(T)

HCRDQ(T);


******* Equation definitions ===========================================================================================

*** Carbon partition between atmosphere and mixed layer ----------------------------------------------------------------

COC1ALLQ(T)      .. COC1ALL(T) =E= COC1ATMBAC+COC1MIXBAC+COC1ATM(T)+COC1MIX(T);

COC1BALQ(T)      .. COC1ATM(T)+COC1MIX(T) =E= COC(T,'1');

COC1ATMQ(T)      .. COC1ATMBAC+COC1ATM(T) =E= CONCO2INT(T)/CNVCO2;


COC1MIXQ(T) $(DYNCARBSWCH > 0.5)
                 .. COC1MIXBAC+COC1MIX(T) =E= CONCO2INT(T)/ALPHAS(T)*(1+AK1S(T)/H(T)+AK1S(T)*AK2S(T)/H(T)**2)*ALPHAC;
* Total DIC concentration is a sum of [CO2(aq)], [HCO3-], and [CO32-], each of which is expressed with the associated
* thermodynamic equilibrium constant. Unit conversion using ALPHAC between mol per Giga L and GtC

ALKLQ(T) $(DYNCARBSWCH > 0.5)
                 .. ALK =E= CONCO2INT(T)/ALPHAS(T)*(AK1S(T)/H(T)+2*AK1S(T)*AK2S(T)/H(T)**2)
                            +BTS/(1+H(T)/AKB(T))+AKW(T)/H(T)-H(T);
* The definition of total alkalinity. In contrast to equation COC1MIXQ,
* the unit conversion using ALPHAC is not necessary because ALK is defined in eq per Giga L.

BUFFERQ(T) $((DYNCARBSWCH < 0.5) and (Ord(T) > 1))
                 .. REVELLE(T)*((COC1MIX(T)-COC1MIX(T-1))/(COC1MIXBAC+(COC1MIX(T-1)+COC1MIX(T))/2))
                    =E= ((COC1ATM(T)-COC1ATM(T-1))/(COC1ATMBAC+(COC1ATM(T-1)+COC1ATM(T))/2));
* The Revelle factor can be calculated when the dynamic calculation for thermodynamic equilibria of carbonate species
* is activated (DYNCARBSWCH > 0.5). When it is not activated (DYNCARBSWCH < 0.5), the Revelle factor is prescribed
* (equation REVELLEQ below). Then, the carbon partition in the mixed layer and the atmosphere is determined.

REVELLEQ(T) $(DYNCARBSWCH < 0.5)
                 .. REVELLE(T) =E= REV1750+REVCOF*(CONCO2INT(T)/(COC1ATMBAC*CNVCO2)-1);
* Revelle and Munk (1977) cited in Mackenzie and Lerman (2006) p.265


* The following equations are written as outputs expressed in conventional units.
PHQ(T)           .. PH(T)         =E= -log10(H(T)/10**9/DENS);

DICQ(T)          .. DIC(T)        =E= (COC1MIXBAC+COC1MIX(T))/ALPHAC/10**9/DENS;

CO2AQQ(T)        .. CO2AQ(T)      =E= CONCO2INT(T)/ALPHAS(T)/10**9/DENS;

HCO3Q(T)         .. HCO3(T)       =E= CONCO2INT(T)/ALPHAS(T)*(AK1S(T)/H(T))/10**9/DENS;

CO32Q(T)         .. CO32(T)       =E= CONCO2INT(T)/ALPHAS(T)*(AK1S(T)*AK2S(T)/H(T)**2)/10**9/DENS;

FRCO2AQQ(T)      .. FRCO2AQ(T)    =E= CO2AQ(T)/DIC(T);

FRHCO3Q(T)       .. FRHCO3(T)     =E= HCO3(T)/DIC(T);

FRCO32Q(T)       .. FRCO32(T)     =E= CO32(T)/DIC(T);


*** CO2 fertilization effect and heterotrophic respiration -------------------------------------------------------------

NPPBACQ(T,NTR)   .. NPPBAC(T,NTR) =E= (BJ(NTR)*TAU(NTR))*FNPP*TNPPFAC(T);

CO2FRTQ(T,NTR)   .. CO2FRT(T,NTR) =E=
                    (BJ(NTR)*TAU(NTR))*(FNPP*BETA(T)*LOG((COC1ATM(T)+COC1ATMBAC)/COC1ATMBAC))*TNPPFAC(T);

HTRRESQ(T,NTR)   .. HTRRES(T,NTR) =E= CTR(T,NTR)/TAU(NTR)*TRESFAC(T);


*** Temperature feedback to atmosphere-ocean carbon cycle --------------------------------------------------------------

TMIXQ1(T) $((DYNCARBSWCH > 0.5) and (OCNCO2TFSWCH > 0.5) and (DOECLIMSWCH > 0.5))
                 .. TMIX(T) =E= TMIX1750+DMPTMIX*(TSS(T)-TSS1750);

TMIXQ2(T) $((DYNCARBSWCH > 0.5) and (OCNCO2TFSWCH > 0.5) and (DOECLIMSWCH < 0.5))
                 .. TMIX(T) =E= TMIX1750+DMPTMIX*(T2M(T)-T2M1750);

TMIXQ3(T) $(((DYNCARBSWCH > 0.5) and (OCNCO2TFSWCH < 0.5)) or (DYNCARBSWCH < 0.5))
                 .. TMIX(T) =E= TMIX1750;


* The following dissociation constants are defined with molar concentration units
* whereas the functional relationships in Millero literature are defined with molal concentration units.
* Thus, the correction terms using DENS are added to adjust the units.

AK1SQ(T) $(DYNCARBSWCH > 0.5) .. -log10(AK1S(T))-(-log10(AK1P(T))) =E=
                                 (AK1A+AK1B/(TMIX(T)+273.15)+AK1C*log((TMIX(T)+273.15)));
* The function for the dissociation constant of bicarbonate in actual seawater (Millero, 2006)
* cf. 9.952E2 in NICCS (Hooss, 2001)

AK1PQ(T) $(DYNCARBSWCH > 0.5) .. -log10(AK1P(T)) =E=
                                 (-126.34048+6320.813/(TMIX(T)+273.15)+19.568224*log((TMIX(T)+273.15)))-9-log10(DENS);

AK1AQ(T) $(DYNCARBSWCH > 0.5) .. AK1A =E= 13.4191*SQRT(SAL)+0.0331*SAL-0.0000533*SAL**2;

AK1BQ(T) $(DYNCARBSWCH > 0.5) .. AK1B =E= -530.123*SQRT(SAL)-6.103*SAL;

AK1CQ(T) $(DYNCARBSWCH > 0.5) .. AK1C =E= -2.06950*SQRT(SAL);


AK2SQ(T) $(DYNCARBSWCH > 0.5) .. -log10(AK2S(T))-(-log10(AK2P(T))) =E=
                                 (AK2A+AK2B/(TMIX(T)+273.15)+AK2C*log((TMIX(T)+273.15)));
* The function for the dissociation constant of bicarbonate in actual seawater (Millero, 2006)
* cf. 7.076E-1 in NICCS (Hooss, 2001)

AK2PQ(T) $(DYNCARBSWCH > 0.5) .. -log10(AK2P(T)) =E=
                                 (-90.18333+5143.692/(TMIX(T)+273.15)+14.613358*log((TMIX(T)+273.15)))-9-log10(DENS);

AK2AQ(T) $(DYNCARBSWCH > 0.5) .. AK2A =E= 21.0894*SQRT(SAL)+0.1248*SAL-0.0003687*SAL**2;

AK2BQ(T) $(DYNCARBSWCH > 0.5) .. AK2B =E= -772.483*SQRT(SAL)-20.051*SAL;

AK2CQ(T) $(DYNCARBSWCH > 0.5) .. AK2C =E= -3.3336*SQRT(SAL);


AKBQ(T) $(DYNCARBSWCH > 0.5)
        .. AKB(T) =E=
        exp[(-8966.90-2890.51*SQRT(SAL)-77.942*SAL+1.726*SAL*SQRT(SAL)-0.0993*POWER(SAL,2))/(TMIX(T)+273.15)
        +(148.0248+137.194*SQRT(SAL)+1.62247*SAL)+(-24.4344-25.085*SQRT(SAL)-0.2474*SAL)*log((TMIX(T)+273.15))
        +0.053105*SQRT(SAL)*(TMIX(T)+273.15)]*1000000000*DENS;
* Dissociation constant of boric acid (as a function of temperature, salinity and density)
* (Millero (1995) equation (52))
* cf. 1.881 in NICCS (Hooss, 2001)


AKWQ(T) $(DYNCARBSWCH > 0.5)
         .. AKW(T) =E= exp[148.9802-13847.26/(TMIX(T)+273.15)-23.6521*Log((TMIX(T)+273.15))
         +(-5.977+118.67/(TMIX(T)+273.15)+1.0495*Log((TMIX(T)+273.15)))*SQRT(SAL)-0.01615*SAL]
         *1000000000000000000*Power(DENS,2);
* Dissociation constant of water (as a function of temperature, salinity and density)
* (Millero (1995) equation (63))
* The original value assigned for this valuable in NICCS deviates from the estimate using Millero's equation
* by a factor of 5. The value calculated in ACC2 is in line with the typical value of Table 7.5 in Pilson (1996).
* cf. 6.463E3 in NICCS (Hooss, 2001)


ALPHASQ(T) $(DYNCARBSWCH > 0.5)
         .. ALPHAS(T) =E=
         [1/exp(-60.2409+93.4517*(100/(TMIX(T)+273.15))+23.3585*Log((TMIX(T)+273.15)/100)
         +SAL*(0.023517-0.023656*((TMIX(T)+273.15)/100)
         +0.0047036*Power((TMIX(T)+273.15)/100,2)))]*1000000/1000000000/DENS;
* Inverse of Henrys Law constant (as a function of temperature, salinity and density)
* (Millero (1995) equation (26))
* cf. 2.988E-2 in NICCS (Hooss, 2001)


*** Temperature feedback to land carbon cycle --------------------------------------------------------------------------

TNPPFACQ1(T) $((DOECLIMSWCH > 0.5) and (LNDCO2TFNPPSWCH > 0.5)) .. TNPPFAC(T) =E= 1+TNPPSENS*(TAL(T)-TAL1750);

TNPPFACQ2(T) $((DOECLIMSWCH < 0.5) and (LNDCO2TFNPPSWCH > 0.5)) .. TNPPFAC(T) =E= 1+TNPPSENS*(T2M(T)-T2M1750);

TNPPFACQ3(T) $((DOECLIMSWCH > 0.5) and (LNDCO2TFNPPSWCH < 0.5)) .. TNPPFAC(T) =E= 1;

TNPPFACQ4(T) $((DOECLIMSWCH < 0.5) and (LNDCO2TFNPPSWCH < 0.5)) .. TNPPFAC(T) =E= 1;


TRESFACQ1(T) $((DOECLIMSWCH > 0.5) and (LNDCO2TFRESSWCH > 0.5)) .. TRESFAC(T) =E= Q10EFF(T)**((TAL(T)-TAL1750)/10);

TRESFACQ2(T) $((DOECLIMSWCH < 0.5) and (LNDCO2TFRESSWCH > 0.5)) .. TRESFAC(T) =E= Q10EFF(T)**((T2M(T)-T2M1750)/10);

TRESFACQ3(T) $((DOECLIMSWCH > 0.5) and (LNDCO2TFRESSWCH < 0.5)) .. TRESFAC(T) =E= 1;

TRESFACQ4(T) $((DOECLIMSWCH < 0.5) and (LNDCO2TFRESSWCH < 0.5)) .. TRESFAC(T) =E= 1;


Q10EFFQ1(T)  $((DOECLIMSWCH > 0.5) and (LNDCO2TFRESSWCH > 0.5) and (LNDCO2TFQ10SWCH > 0.5))
                                                         .. Q10EFF(T) =E= Q10BAS-Q10TMP*(TAL(T)-TAL1750);

Q10EFFQ2(T)  $((DOECLIMSWCH < 0.5) and (LNDCO2TFRESSWCH > 0.5) and (LNDCO2TFQ10SWCH > 0.5))
                                                         .. Q10EFF(T) =E= Q10BAS-Q10TMP*(T2M(T)-T2M1750);

Q10EFFQ3(T)  $((DOECLIMSWCH > 0.5) and (LNDCO2TFRESSWCH > 0.5) and (LNDCO2TFQ10SWCH < 0.5))
                                                         .. Q10EFF(T) =E= Q10BAS;

Q10EFFQ4(T)  $((DOECLIMSWCH < 0.5) and (LNDCO2TFRESSWCH > 0.5) and (LNDCO2TFQ10SWCH < 0.5))
                                                         .. Q10EFF(T) =E= Q10BAS;

Q10EFFQ5(T)  $((DOECLIMSWCH > 0.5) and (LNDCO2TFRESSWCH < 0.5))
                                                         .. Q10EFF(T) =E= 1;

Q10EFFQ6(T)  $((DOECLIMSWCH < 0.5) and (LNDCO2TFRESSWCH < 0.5))
                                                         .. Q10EFF(T) =E= 1;


*** Ocean CO2 uptake ---------------------------------------------------------------------------------------------------

CO2OCNUPANTQ1(T) $((MODEFLAG > 0.5) and (Ord(T) = 1))
                                 .. CO2OCNUPANT(T) =E= 0;

CO2OCNUPANTQ2(T) $((MODEFLAG < 0.5) and (Ord(T) = 1))
                                 .. CO2OCNUPANT(T) =E= ((SUMCOC(T+1)-SUMCOC(T))+(SUMCOC(T)-SUMCOC1999))/2/DELTAT;

CO2OCNUPANTQ3(T) $((Ord(T) > 1) and (Ord(T) < Card(T)))
                                 .. CO2OCNUPANT(T) =E= ((SUMCOC(T+1)-SUMCOC(T))+(SUMCOC(T)-SUMCOC(T-1)))/2/DELTAT;

CO2OCNUPANTQ4(T) $(Ord(T) = Card(T))
                                 .. CO2OCNUPANT(T) =E= 0;


CO2OCNUPTTLQ1(T) $((MODEFLAG > 0.5) and (Ord(T) = 1))
                                 .. CO2OCNUPTTL(T) =E= 0;

CO2OCNUPTTLQ2(T) $((MODEFLAG < 0.5) and (Ord(T) = 1) and (CO2UPPRESWCH > 0.5))
                                 .. CO2OCNUPTTL(T) =E= CO2OCNUPPRE(T)+CO2OCNUPANT(T);

CO2OCNUPTTLQ3(T) $((MODEFLAG < 0.5) and (Ord(T) = 1) and (CO2UPPRESWCH < 0.5))
                                 .. CO2OCNUPTTL(T) =E= CO2OCNUPANT(T);

CO2OCNUPTTLQ4(T) $((Ord(T) > 1) and (Ord(T) < Card(T)) and (CO2UPPRESWCH > 0.5))
                                 .. CO2OCNUPTTL(T) =E= CO2OCNUPPRE(T)+CO2OCNUPANT(T);

CO2OCNUPTTLQ5(T) $((Ord(T) > 1) and (Ord(T) < Card(T)) and (CO2UPPRESWCH < 0.5))
                                 .. CO2OCNUPTTL(T) =E= CO2OCNUPANT(T);

CO2OCNUPTTLQ6(T) $(Ord(T) = Card(T))
                                 .. CO2OCNUPTTL(T) =E= 0;
* The CO2 uptake in the year 1750 and 2100 cannot be defined.
* The CO2 uptake in the year 2000 cannot be calculated in the past mode alone.
* The future mode utilizes the information in the carbon storage in 1999 (CO2OCNUPANTQ2 and CO2OCNUPTTLQ2)
* to calculate the CO2 uptake in the year 2000.


*** Land CO2 uptake ----------------------------------------------------------------------------------------------------

CO2LNDUPANTQ1(T) $((MODEFLAG > 0.5) and (Ord(T) = 1))
                                 .. CO2LNDUPANT(T) =E= 0;

CO2LNDUPANTQ2(T) $((MODEFLAG < 0.5) and (Ord(T) = 1))
                                 .. CO2LNDUPANT(T) =E= ((SUMCTR(T+1)-SUMCTR(T))+(SUMCTR(T)-SUMCTR1999))/2/DELTAT;

CO2LNDUPANTQ3(T) $((Ord(T) > 1) and (Ord(T) < Card(T)))
                                 .. CO2LNDUPANT(T) =E= ((SUMCTR(T+1)-SUMCTR(T))+(SUMCTR(T)-SUMCTR(T-1)))/2/DELTAT;

CO2LNDUPANTQ4(T) $(Ord(T) = Card(T))
                                 .. CO2LNDUPANT(T) =E= 0;


CO2LNDUPTTLQ1(T) $((MODEFLAG > 0.5) and (Ord(T) = 1))
                                 .. CO2LNDUPTTL(T) =E= 0;

CO2LNDUPTTLQ2(T) $((MODEFLAG < 0.5) and (Ord(T) = 1) and (CO2UPPRESWCH > 0.5))
                                 .. CO2LNDUPTTL(T) =E= CO2LNDUPPRE(T)+CO2LNDUPANT(T);

CO2LNDUPTTLQ3(T) $((MODEFLAG < 0.5) and (Ord(T) = 1) and (CO2UPPRESWCH < 0.5))
                                 .. CO2LNDUPTTL(T) =E= CO2LNDUPANT(T);

CO2LNDUPTTLQ4(T) $((Ord(T) > 1) and (Ord(T) < Card(T)) and (CO2UPPRESWCH > 0.5))
                                 .. CO2LNDUPTTL(T) =E= CO2LNDUPPRE(T)+CO2LNDUPANT(T);

CO2LNDUPTTLQ5(T) $((Ord(T) > 1) and (Ord(T) < Card(T)) and (CO2UPPRESWCH < 0.5))
                                 .. CO2LNDUPTTL(T) =E= CO2LNDUPANT(T);

CO2LNDUPTTLQ6(T) $(Ord(T) = Card(T))
                                 .. CO2LNDUPTTL(T) =E= 0;
* The CO2 uptake in the year 1750 and 2100 cannot be defined.
* The CO2 uptake in the year 2000 cannot be calculated in the past mode alone.
* The future mode utilizes the information in the carbon storage in 1999 (CO2LNDUPANTQ2 and CO2LNDUPTTLQ2)
* to calculate the CO2 uptake in the year 2000.


*** Preindustrial oceanic and land CO2 uptake --------------------------------------------------------------------------

CO2OCNUPPRECNSTQ(T) $(Ord(T) > 1) .. CO2OCNUPPRE(T) =E= CO2OCNUPPRE(T-1);

CO2LNDUPPRECNSTQ(T) $(Ord(T) > 1) .. CO2LNDUPPRE(T) =E= CO2LNDUPPRE(T-1);


*** Governing equations for carbon cycle box model (Based on Eq. 90 in Hooss et al., 1999) -----------------------------
*
* The equilibrium assumptions for the short-term land reservoirs (Hooss et al., 1999) are lifted
* because the time step in ACC2 has been shortened from five years to one year.

COCEQ1(T) $(Ord(T) > 1)       .. (COC(T,'1')-COC(T-1,'1'))/DELTAT =E= 0.5*(FLXCOC1(T)+FLXCOC1(T-1));

COCEQ2(T) $(Ord(T) > 1)       .. (COC(T,'2')-COC(T-1,'2'))/DELTAT =E= 0.5*(FLXCOC2(T)+FLXCOC2(T-1));

COCEQ3(T) $(Ord(T) > 1)       .. (COC(T,'3')-COC(T-1,'3'))/DELTAT =E= 0.5*(FLXCOC3(T)+FLXCOC3(T-1));

COCEQ4(T) $(Ord(T) > 1)       .. (COC(T,'4')-COC(T-1,'4'))/DELTAT =E= 0.5*(FLXCOC4(T)+FLXCOC4(T-1));

CTREQ(T,NTR) $(Ord(T) GT 1)   .. (CTR(T,NTR)-CTR(T-1,NTR))/DELTAT =E= 0.5*(FLXCTR(T,NTR)+FLXCTR(T-1,NTR));


FLXCOC1Q(T)                   .. FLXCOC1(T) =E= ANT(T,'CO2')-(CO2OCNUPPRE(T)+CO2LNDUPPRE(T))$(CO2UPPRESWCH > 0.5)
                                 -ETA('2')/HS*COC1MIX(T)+ETA('2')/HL('2')*COC(T,'2')-SUM(NTR,FLXCTR(T,NTR));

FLXCOC2Q(T)                   .. FLXCOC2(T) =E= ETA('2')/HS*COC1MIX(T)-(ETA('2')+ETA('3'))/HL('2')*COC(T,'2')+
                                 ETA('3')/HL('3')*COC(T,'3');

FLXCOC3Q(T)                   .. FLXCOC3(T) =E= ETA('3')/HL('2')*COC(T,'2')-(ETA('3')+ETA('4'))/HL('3')*COC(T,'3')+
                                 ETA('4')/HL('4')*COC(T,'4');

FLXCOC4Q(T)                   .. FLXCOC4(T) =E= ETA('4')/HL('3')*COC(T,'3')-ETA('4')/HL('4')*COC(T,'4');

FLXCTRQ(T,NTR)                .. FLXCTR(T,NTR) =E= NPPBAC(T,NTR)+CO2FRT(T,NTR)-HTRRES(T,NTR);


SUMCOCQ(T)                    .. SUMCOC(T) =E= SUM(NOC, COC(T,NOC))-COC1ATM(T);

SUMCTRQ(T)                    .. SUMCTR(T) =E= SUM(NTR, CTR(T,NTR));


CTRRELQ(T,NTR)                .. CTRREL(T,NTR) =E= CTR(T,NTR)-CTR1750(NTR);


CUMCO2Q1(T) $[(MODEFLAG > 0.5) and (Ord(T) = 1)]
                              .. CUMCO2(T) =E= ANT(T,'CO2');

CUMCO2Q2(T) $[(MODEFLAG > 0.5) and (Ord(T) > 1)]
                              .. CUMCO2(T) =E= CUMCO2(T-1)+ANT(T,'CO2');

CUMCO2Q3(T) $[(MODEFLAG < 0.5) and (Ord(T) = 1)]
                              .. CUMCO2(T) =E= CUMCO21999+ANT(T,'CO2');

CUMCO2Q4(T) $[(MODEFLAG < 0.5) and (Ord(T) > 1)]
                              .. CUMCO2(T) =E= CUMCO2(T-1)+ANT(T,'CO2');


*** Concentrations -----------------------------------------------------------------------------------------------------
*
* Below is a technical note on two different variables representing the atmospheric CO2 concentration
* (due to the treatment of ENSO influence on carbon cycle).
* In the past mode, CON(T,'CO2') includes the ENSO influence to the atmospheric CO2 concentration
* whereas CONCO2INT(T) does not include it.
* In the future mode, as the ENSO influence cannot be considered, these two valuables have the same values over time.
* Thus, CON(T,'CO2') in the future mode in 2000 is different from that in the past mode in 2000
* by the ENSO contribution to the CO2 concentration.
* It also follows that CONCO2INT('2000') rather than CON('2000','CO2') is transferred
* from the past mode to the future mode (see ACC2_init_future.inc).
* In the past mode, the ENSO effect is added in the final atmospheric CO2 concentration using the regression line
* (see ACC2_settings_past.inc) after all the calculations involving carbon allocations are done.

CO2CONC1(T)$((ENSOCO2SWCH > 0.5) and (MODEFLAG > 0.5)) .. CON(T,'CO2') =E= CONCO2INT(T)+DCONCO2ENSO(T)*CNVCO2;
* Atmospheric carbon dioxide concentration (with ENSO contribution in the past mode)

CO2CONC2(T)$((ENSOCO2SWCH < 0.5) and (MODEFLAG > 0.5)) .. CON(T,'CO2') =E= CONCO2INT(T);
* Atmospheric carbon dioxide concentration (without ENSO contribution in the past mode)

CO2CONC3(T)$((ENSOCO2SWCH > 0.5) and (MODEFLAG < 0.5) and (Ord (T) = 1))
                                                       .. CON(T,'CO2') =E= CONCO2INT(T)+DCONCO2ENSO2000*CNVCO2;
* Atmospheric carbon dioxide concentration in 2000 (with ENSO contribution in the future mode)

CO2CONC4(T)$((ENSOCO2SWCH > 0.5) and (MODEFLAG < 0.5) and (Ord (T) > 1))
                                                       .. CON(T,'CO2') =E= CONCO2INT(T);
* Atmospheric carbon dioxide concentration after 2000 (with ENSO contribution in the future mode)
* CON(T,'CO2') is different between the past and future mode because ENSO effect is removed in the future mode.

CO2CONC5(T)$((ENSOCO2SWCH < 0.5) and (MODEFLAG < 0.5)) .. CON(T,'CO2') =E= CONCO2INT(T);
* Atmospheric carbon dioxide concentration (without ENSO contribution in the past mode)

DCONCO2ENSOQ(T) $(MODEFLAG < 0.5)      .. DCONCO2ENSO(T) =E= 0;
* To fix DCONCO2ENSO at zero not to consider the ENSO effect in the future mode


CONCH4Q(T) $(Ord(T) GT 1)     .. (CON(T,'CH4')-CON(T-1,'CH4'))/DELTAT =E= 0.5*
                                 ((1/CNVCH4*(ANT(T,'CH4')+ NATCH4)-
                                 (1/(1/(CON(T,'rOH')/TAUCH4OH+1/TAUCH4SS)))*CON(T,'CH4'))+
                                 (1/CNVCH4*(ANT(T-1,'CH4')+ NATCH4)-
                                 (1/(1/(CON(T-1,'rOH')/TAUCH4OH+1/TAUCH4SS)))*CON(T-1,'CH4')));

CONN2OQ(T) $(Ord(T) GT 1)     .. (CON(T,'N2O')-CON(T-1,'N2O'))/DELTAT =E= 0.5*
                                 ((1/CNVN2O*(ANT(T,'N2O')+ NATN2O)-(1/TAUN2O*(CON(T,'N2O')/CONN2O2000R)
                                 **(-SENTAUN2O)*CON(T,'N2O')))+
                                 (1/CNVN2O*(ANT(T-1,'N2O')+ NATN2O)-(1/TAUN2O*(CON(T-1,'N2O')/CONN2O2000R)
                                 **(-SENTAUN2O)*CON(T-1,'N2O'))));

CONOHQ(T)                     .. CON(T,'rOH') =E= exp(-0.32*(LOG(CON(T,'CH4'))-LOG(CON('2000','CH4')))
                                 +0.0042*(ANT(T,'NOx')-ANT('2000','NOx'))-0.000105*(ANT(T,'CO')-ANT('2000','CO'))
                                 -0.000315*(ANT(T,'VOC')-ANT('2000','VOC')));
* Note: LOG is natural logarithm in GAMS.


*** Radiative forcings -------------------------------------------------------------------------------------------------

RADCO2Q(T)       .. RAD(T,'CO2') =E= RAD2XCO2/LOG(2.0)*LOG((CON(T,'CO2')/CNVCO2)/COC1ATMBAC);

RADCH4Q(T)       .. RAD(T,'CH4') =E= RHOCH4*(SQRT(CON(T,'CH4'))-SQRT(CONCH4PRE))
                                     -OVERLFAC1*LOG(1+OVERLFAC2*(CON(T,'CH4')*CONN2OPRE)**OVERLEXP1+
                                     OVERLFAC3*CON(T,'CH4')*(CON(T,'CH4')*CONN2OPRE)**OVERLEXP2)
                                     +OVERLFAC1*LOG(1+OVERLFAC2*(CONCH4PRE*CONN2OPRE)**OVERLEXP1+
                                     OVERLFAC3*CONCH4PRE*(CONCH4PRE*CONN2OPRE)**OVERLEXP2);

RADN2OQ(T)       .. RAD(T,'N2O') =E= RHON2O*(SQRT(CON(T,'N2O'))-SQRT(CONN2OPRE))
                                     -OVERLFAC1*LOG(1+OVERLFAC2*(CONCH4PRE*CON(T,'N2O'))**OVERLEXP1+
                                     OVERLFAC3*CONCH4PRE*(CONCH4PRE*CON(T,'N2O'))**OVERLEXP2)
                                     +OVERLFAC1*LOG(1+OVERLFAC2*(CONCH4PRE*CONN2OPRE)**OVERLEXP1+
                                     OVERLFAC3*CONCH4PRE*(CONCH4PRE*CONN2OPRE)**OVERLEXP2);

RADHALTQ(T)      .. RADHAL(T) =E= SUM(HALR,RAD(T,HALR));

RADODSQ(T)       .. RADODS(T) =E= SUM(ODSR,RAD(T,ODSR));

RADSO2Q(T)       .. RAD(T,'SO2') =E=
                     [SCLRAD*RADSO2DR2000*ANT(T,'SO2')/ANT('2000','SO2')]$(MODEFLAG > 0.5)
                    +[SCLRAD*RADSO2DR2000*ANT(T,'SO2')/ANT('2000','SO2')]$((MODEFLAG < 0.5) and (MISTRNSWCH < 0.5))
                    +[SCLRAD*RADSO2DR2000*ANT(T,'SO2')/ANT('2000','SO2')*(1+RADMISAVE/RADAERAVE)]
                     $((MODEFLAG < 0.5) and (MISTRNSWCH > 0.5));
*RADSO2Q(T)       .. RAD(T,'SO2') =E= 0;
* Direct effect of sulfur aerosol
* Note a difference in the SO2 emission estimates in 2000 between RIVM (used for past run)
* and IPCC SRES (used for future run).

RADAERINDQ(T)    .. RAD(T,'AERIND') =E=
                     [SCLRAD*(RADSO2ID2000*LOG((NATSO2+ANT(T,'SO2'))/NATSO2)/LOG((NATSO2+ANT('2000','SO2'))/NATSO2))]
                     $(MODEFLAG > 0.5)
                    +[SCLRAD*(RADSO2ID2000*LOG((NATSO2+ANT(T,'SO2'))/NATSO2)/LOG((NATSO2+ANT('2000','SO2'))/NATSO2))]
                     $((MODEFLAG < 0.5) and (MISTRNSWCH < 0.5))
                    +[SCLRAD*(RADSO2ID2000*LOG((NATSO2+ANT(T,'SO2'))/NATSO2)/LOG((NATSO2+ANT('2000','SO2'))/NATSO2))
                     *(1+RADMISAVE/RADAERAVE)]$((MODEFLAG < 0.5) and (MISTRNSWCH > 0.5));
*RADAERINDQ(T)    .. RAD(T,'AERIND') =E= 0;
* (First and second) indirect effects of all aerosols

RADCRBAERQ(T)    .. RAD(T,'CRBAER') =E= [SCLRAD*RADCRBAER2000*ANT(T,'CO')/ANT('2000','CO')]$(MODEFLAG > 0.5)
                                       +[SCLRAD*RADCRBAER2000*ANT(T,'CO')/ANT('2000','CO')]
                                         $((MODEFLAG < 0.5) and (MISTRNSWCH < 0.5))
                                       +[SCLRAD*RADCRBAER2000*ANT(T,'CO')/ANT('2000','CO')*(1+RADMISAVE/RADAERAVE)]
                                         $((MODEFLAG < 0.5) and (MISTRNSWCH > 0.5));
*RADCRBAERQ(T)    .. RAD(T,'CRBAER') =E= 0;
* Direct effect of carbonaceous (BC+OC) aerosol forcing

RADH2OSTRQ(T)    .. RAD(T,'H2OSTR') =E= ALPHCH4H2O*RHOCH4*(SQRT(CON(T,'CH4'))-SQRT(CONCH4PRE));


*** Principal components of Empirical Orthogonal Function --------------------------------------------------------------
*
* Original PC values of Hooss et al. (2001) are scaled by a factor depending on the mean value
* captured in the EOF pattern. The following scale factors are used to obtain an EOF mean value of 1.

NORRADQ(T)       .. NORRAD(T) =E= RADTTL(T)/RAD2XCO2;

SCALEMEANTQ      .. SCALEMEANT =E= 4.67384/40.4741;

SCALEMEANSQ      .. SCALEMEANS =E= 1.52056/10.1464;


SQROOTQ(T)       .. SQROOT(T) =E= SQRT(33.9031639 + 8.36722152 * SENS(T) + SENS(T)**2);

RPCTT1Q(T)       .. RPCTT(T,'1') =E= 91.0141709*SENS(T)/(5.8226761+SENS(T)-SQROOT(T));

RPCTT2Q(T)       .. RPCTT(T,'2') =E= 91.0141709*SENS(T)/(5.8226761+SENS(T)+SQROOT(T));

RPCTA1Q(T)       .. RPCTA(T,'1') =E= SENS(T)/SQROOT(T)*(-2.9113173 + 0.5*(SQROOT(T)+SENS(T)));

RPCTA2Q(T)       .. RPCTA(T,'2') =E= SENS(T)/SQROOT(T)*(2.9113173 + 0.5*(SQROOT(T)-SENS(T)));


SENSSCALEQ(T)    .. SENSSCALE(T) =E= SENS(T)/(SUM(NPCRSEN,RPCTAP(NPCRSEN))*SCALEMEANT);


*** Temperature change calculated by IRF -------------------------------------------------------------------------------
*
* Below is a technical note on two different variables representing the surface air temperature change
* (due to the treatment of ENSO influence on temperature).
* In the past mode, the ENSO effect is added in the final temperature change (T2M(T))
* after all the energy balance model or IRF calculation are done.
* In the future mode, the ENSO influence is not added in T2M(T).
* Thus, T2M('2000') in the future mode is different from that in the past mode by the ENSO contribution
* to the temperature.
* On the other hand, GMT('2000') and PCT('2000') are consistent in the past and future modes
* because they are internally calculated in DEOCLIM and IRF and are not influenced by ENSO.
* Thus, GMT('2000') and PCT('2000') rather than T2M('2000') is transferred in the initial across the operation modes.

PCTRESQ(T,NPCR) $(Ord(T) > 1) .. (PCTRES(T,NPCR)-PCTRES(T-1,NPCR))/DELTAT =E= 0.5/RPCTT(T,NPCR)*
                                 ((RPCTA(T,NPCR)*NORRAD(T)-PCTRES(T,NPCR))+
                                 (RPCTA(T-1,NPCR)*NORRAD(T-1)-PCTRES(T-1,NPCR)));

PCTQ(T)                       .. PCT(T) =E= SUM(NPCR,PCTRES(T,NPCR));

T2MQ1(T) $((DOECLIMSWCH > 0.5) and (ENSOT2MSWCH > 0.5) and (MODEFLAG > 0.5)) .. T2M(T) =E= GMT(T)+DT2MENSO(T);

T2MQ2(T) $((DOECLIMSWCH < 0.5) and (ENSOT2MSWCH > 0.5) and (MODEFLAG > 0.5)) .. T2M(T) =E= PCT(T)+DT2MENSO(T);

T2MQ3(T) $((DOECLIMSWCH > 0.5) and (ENSOT2MSWCH < 0.5) and (MODEFLAG > 0.5)) .. T2M(T) =E= GMT(T);

T2MQ4(T) $((DOECLIMSWCH < 0.5) and (ENSOT2MSWCH < 0.5) and (MODEFLAG > 0.5)) .. T2M(T) =E= PCT(T);

T2MQ5(T) $((DOECLIMSWCH > 0.5)                         and (MODEFLAG < 0.5)) .. T2M(T) =E= GMT(T);

T2MQ6(T) $((DOECLIMSWCH < 0.5)                         and (MODEFLAG < 0.5)) .. T2M(T) =E= PCT(T);

T2MRT1990Q(T)                 .. T2MRT1990(T) =E= T2M(T)-T2M1990;


*** Sea level rise calculated by IRF (due to thermal expansion of oceans) ----------------------------------------------
*
* Evolution of first principle component of sea level rise

PCSRESQ(T,NPCR) $(Ord(T) > 1) .. (PCSRES(T,NPCR)-PCSRES(T-1,NPCR))/DELTAT =E= 0.5/RPCST(NPCR)*
                                 ((RPCSA(NPCR)*SENSSCALE(T)*SCALEMEANS*NORRAD(T)-PCSRES(T,NPCR))+
                                 (RPCSA(NPCR)*SENSSCALE(T-1)*SCALEMEANS*NORRAD(T-1)-PCSRES(T-1,NPCR)));

PCSQ(T)                       .. PCS(T) =E= SUM(NPCR,PCSRES(T,NPCR));

PCSQ1990(T) $(MODEFLAG < 0.5) .. PCSRT1990(T) =E= PCS(T)-PCS1990;


*** Sea level rise calculated by IRF (due to factors other than thermal expansion) -------------------------------------
*
* Equations with the declarations ending with F(=Future) are equivalent to the corresponding equations with P(=Past)

SLRGLQ1(T) $[(MODEFLAG GT 0.5) AND [(Ord(T) EQ 241) AND (DELTAT = 1)]] .. SLRGLWT(T) =E= SLRGLWT1990;

SLRGLQ2P(T) $[(MODEFLAG GT 0.5) AND [(Ord(T) GT 241) AND (DELTAT = 1)]]
                              .. (SLRGLWT(T)-SLRGLWT(T-1))/DELTAT =E= 0.5*((T2M1990AOGCM+T2M19ST+T2M(T)-T2M1990)*SENSGL
                                 +(T2M1990AOGCM+T2M19ST+T2M(T-1)-T2M1990)*SENSGL);

SLRGLQ2F(T) $[(MODEFLAG LT 0.5) AND (Ord(T) GT 1)]
                              .. (SLRGLWT(T)-SLRGLWT(T-1))/DELTAT =E= 0.5*((T2M1990AOGCM+T2M19ST+T2M(T)-T2M1990)*SENSGL
                                 +(T2M1990AOGCM+T2M19ST+T2M(T-1)-T2M1990)*SENSGL);

SLRGLQ3P(T) $[(MODEFLAG GT 0.5) AND [(Ord(T) GT 240) AND (DELTAT = 1)]]
                              .. SLRGLWC(T) =E= COEFGL1*SLRGLWT(T)+COEFGL2*Power(SLRGLWT(T),2);
* Power function is used for the case when SLRGLWT(T) becomes negative.
* Operation '**' cannot be used for a negative number.

SLRGLQ3F(T) $[(MODEFLAG LT 0.5) AND (Ord(T) GT 1)]
                              .. SLRGLWC(T) =E= COEFGL1*SLRGLWT(T)+COEFGL2*Power(SLRGLWT(T),2);


SLRGLQ1990(T) $(MODEFLAG LT 0.5) .. SLRGLWCRT1990(T) =E= SLRGLWC(T)-SLRGLWC1990;


SLRGRQ1(T) $[(MODEFLAG GT 0.5) AND [(Ord(T) EQ 241) AND (DELTAT = 1)]] .. SLRGR(T) =E= SLRGR1990;

SLRGRQ2P(T) $[(MODEFLAG GT 0.5) AND [(Ord(T) GT 241) AND (DELTAT = 1)]]
                              .. (SLRGR(T)-SLRGR(T-1))/DELTAT =E= 0.5*((T2M1990AOGCM+T2M(T)-T2M1990)*SENSGR
                                 +(T2M1990AOGCM+T2M(T-1)-T2M1990)*SENSGR);

SLRGRQ2F(T) $[(MODEFLAG LT 0.5) AND (Ord(T) GT 1)]
                              .. (SLRGR(T)-SLRGR(T-1))/DELTAT =E= 0.5*((T2M1990AOGCM+T2M(T)-T2M1990)*SENSGR
                                 +(T2M1990AOGCM+T2M(T-1)-T2M1990)*SENSGR);

SLRGRQ1990(T) $(MODEFLAG LT 0.5) .. SLRGRRT1990(T) =E= SLRGR(T)-SLRGR1990;


SLRANQ1(T) $[(MODEFLAG GT 0.5) AND [(Ord(T) EQ 241) AND (DELTAT = 1)]] .. SLRAN(T) =E= SLRAN1990;

SLRANQ2P(T) $[(MODEFLAG GT 0.5) AND [(Ord(T) GT 241) AND (DELTAT = 1)]]
                              .. (SLRAN(T)-SLRAN(T-1))/DELTAT =E= 0.5*((T2M1990AOGCM+T2M(T)-T2M1990)*SENSAN
                                 +(T2M1990AOGCM+T2M(T-1)-T2M1990)*SENSAN);

SLRANQ2F(T) $[(MODEFLAG LT 0.5) AND (Ord(T) GT 1)]
                              .. (SLRAN(T)-SLRAN(T-1))/DELTAT =E= 0.5*((T2M1990AOGCM+T2M(T)-T2M1990)*SENSAN
                                 +(T2M1990AOGCM+T2M(T-1)-T2M1990)*SENSAN);

SLRANQ1990(T) $(MODEFLAG LT 0.5) .. SLRANRT1990(T) =E= SLRAN(T)-SLRAN1990;


SLRQP(T) $[(MODEFLAG GT 0.5) AND [(Ord(T) GT 240) AND (DELTAT = 1)]] .. SLRRT1750(T) =E= PCS(T)+SLRGLWC(T)
                                 +SLRGR(T)+SLRAN(T);

SLRQF(T) $(MODEFLAG LT 0.5) .. SLRRT1750(T) =E= PCS(T)+SLRGLWC(T)+SLRGR(T)+SLRAN(T);

SLRQ1990(T) $(MODEFLAG LT 0.5) .. SLRRT1990(T) =E= SLRRT1750(T)-SLR1990+(SLRI+SLRP+SLRS)*((Ord(T)-1)*DELTAT+10);
* +10 takes account of the sea level rise due to the minor components (SLRI, SLRP, and SLRS) between 1990 and 2000.
* Note that Ord('2000') returns 1.


*** Uncertain parameters -----------------------------------------------------------------------------------------------
*
* Parameters BETA and SENS are time-invariant by default as fixed below.
* The time dimension is defined for these parameters only for experiments.

BETACONT(T) $(Ord(T) > 1) .. BETA(T) =E= BETA(T-1);

SENSCONT(T) $(Ord(T) > 1) .. SENS(T) =E= SENS(T-1);
* The Impulse Response Function rather than DOECLIM (DOECLIMSWCH = 0) has to be used to assume time-dependent
* climate sensitivity.


* ----------------------------------------------------------------------------------------------------------------------
* Note: The equations below are taken from DOECLIM (Kriegler, 2005)
* ----------------------------------------------------------------------------------------------------------------------

*** Heat transfer coefficients -----------------------------------------------------------------------------------------
* (equations A.19-A.21 in Kriegler, 2005)

*                   Land ocean heat transfer: Offset
KP1Q             .. KP1 =E= BK * RLAM * FLA / (RLAM*FLA-AK*(RLAM-BSI));

*                   Land ocean heat transfer: Slope
KP2Q             .. KP2 =E= - AK * FLA * (RLAM*FLA+(1-FLA)*BSI) / (RLAM*FLA-AK*(RLAM-BSI)) * RAD2XCO2 ;

*                   Land feedback strength: Offset
LP1Q             .. LP1 =E= -BK * (RLAM-BSI) /  (RLAM*FLA-AK*(RLAM-BSI)) ;

*                   Land feedback strength: Slope
LP2Q             .. LP2 =E= FLA * (RLAM*FLA+(1-FLA)*BSI) / (RLAM*FLA-AK*(RLAM-BSI)) * RAD2XCO2 ;

*                   Ocean feedback strength: Offset
OP1Q             .. OP1 =E= BK * RLAM*FLA / (1-FLA) * (RLAM-BSI) /  (RLAM*FLA-AK*(RLAM-BSI)) ;

*                   Ocean feedback strength: Slope
OP2Q             .. OP2 =E= (RLAM*FLA - AK*(RLAM-BSI)/(1-FLA)) *
                            (RLAM*FLA+(1-FLA)*BSI) / (RLAM*FLA-AK*(RLAM-BSI)) * RAD2XCO2 ;


RADSOLLQ(T)      .. RADSOLL(T) =E= DQL*RADSOL(T)*1/[(1-(FLA*ALBL+(1-FLA)*ALBO))/4]*(1-ALBL)/4;

RADSOLOQ(T)      .. RADSOLO(T) =E= DQO*RADSOL(T)*1/[(1-(FLA*ALBL+(1-FLA)*ALBO))/4]*(1-ALBO)/4;

RADTTLLQ(T)      .. RADTTLL(T) =E= RADANTL(T)+RADVOL(T)+RADSOLL(T)+RADMIS(T);

RADTTLOQ(T)      .. RADTTLO(T) =E= RADANTO(T)+RADVOL(T)+RADSOLO(T)+RADMIS(T);


*** Climate model parameters -------------------------------------------------------------------------------------------
* (see equations A.19-A.21 in Kriegler, 2005)

CLIMP0(T)        .. KSL*SENS(T) =E= KP1*SENS(T) + KP2 ;

CLIMP1(T)        .. LLAM*SENS(T) =E= LP1*SENS(T) + LP2 ;

CLIMP2(T)        .. OLAM*SENS(T) =E= OP1*SENS(T) + OP2 ;

KHEAT            .. KDIF =E= CKD * KAPPA ;


*** Convolution kernel of integration ----------------------------------------------------------------------------------

* Equals a_{i-j}-a_{i-j+1} in Equation A.26 in Kriegler, 2005
* See definition of a_{i-j} on top of page XIII.
*
* Infinite series can be included up to third order
*
* KER1(I): Solution for infinitely deep ocean and first order bottom correction term
* KER2(I): Second order bottom correction term
* KER3(I): Third order bottom correction term
*
* ATTENTION! ERRORF IN GAMS IS NOT IDENTICAL TO STANDARD DEFINITION OF ERF
* IT IS: ERF(x) = 2*ERRORF(SQRT(2)*x)-1 !!
* The following convolution equations were modified by Katsumasa Tanaka to fit past and future modes.

CONVOL10(I) $(ORD(I) = 1)     .. KER1(I) =E= 2*(1-2*Exp(-BOT**2/KDIF/DELTAT)) +
                                  4 * BOT*SQRT(PI/KDIF/DELTAT) *
                                  2 * (1 - Errorf(BOT*SQRT(2/KDIF/DELTAT)) ) ;

CONVOL1(I) $(ORD(I) > 1)      .. KER1(I) =E= 2*SQRT(ORD(I))*(1-2*Exp(-BOT**2/KDIF/DELTAT/ORD(I))) -
                                  2*SQRT(ORD(I)-1)*(1-2*Exp(-BOT**2/KDIF/DELTAT/(ORD(I)-1))) +
                                  4 * BOT*SQRT(PI/KDIF/DELTAT) *
                                  2 * ( Errorf(BOT*SQRT(2/KDIF/DELTAT/(ORD(I)-1))) -
                                  Errorf(BOT*SQRT(2/KDIF/DELTAT/ORD(I))) ) ;

*CONVOL20(I) $(ORD(I) = 1)     .. KER2(I) =E= 4*Exp(-4*BOT**2/KDIF/DELTAT) -
*                                  8 * BOT*SQRT(PI/KDIF/DELTAT) *
*                                  2 * (1 - Errorf(2*BOT*SQRT(2/KDIF/DELTAT)) ) ;
*
*CONVOL2(I) $(ORD(I) > 1)      .. KER2(I) =E= 2*SQRT(ORD(I))*2*Exp(-4*BOT**2/KDIF/DELTAT/ORD(I)) -
*                                  2*SQRT(ORD(I)-1)*2*Exp(-4*BOT**2/KDIF/DELTAT/(ORD(I)-1)) -
*                                  8 * BOT*SQRT(PI/KDIF/DELTAT) *
*                                  2 * ( Errorf(2*BOT*SQRT(2/KDIF/DELTAT/(ORD(I)-1))) -
*                                  Errorf(2*BOT*SQRT(2/KDIF/DELTAT/ORD(I))) ) ;
*
*CONVOL30(I) $(ORD(I) = 1)     .. KER3(I) =E= -4*Exp(-9*BOT**2/KDIF/DELTAT) +
*                                  12 * BOT*SQRT(PI/KDIF/DELTAT) *
*                                  2 * (1 - Errorf(3*BOT*SQRT(2/KDIF/DELTAT)) ) ;
*
*CONVOL3(I) $(ORD(I) > 1)      .. KER3(I) =E= -2*SQRT(ORD(I))*2*Exp(-9*BOT**2/KDIF/DELTAT/ORD(I)) +
*                                  2*SQRT(ORD(I)-1)*2*Exp(-9*BOT**2/KDIF/DELTAT/(ORD(I)-1)) +
*                                  12 * BOT*SQRT(PI/KDIF/DELTAT) *
*                                  2 * ( Errorf(3*BOT*SQRT(2/KDIF/DELTAT/(ORD(I)-1))) -
*                                  Errorf(3*BOT*SQRT(2/KDIF/DELTAT/ORD(I))) ) ;
*
* Choose number of bottom correction terms
* to 1st order: KER(I) =E=  KER1(I) ;
* to 2nd order: KER(I) =E=  KER1(I) + KER2(I) ;
* to 3rd order: KER(I) =E=  KER1(I) + KER2(I) + KER3(I) ;

CONVOLS(I)                    .. KER(I) =E= KER1(I) ;


*** Temperature memory -------------------------------------------------------------------------------------------------


CLIMHLPS0(T) $((MODEFLAG > 0.5) AND (ORD(T) LE 2))
                              .. TSUM(T) =E= 0 ;

CLIMHLPS(T)  $((MODEFLAG > 0.5) AND (ORD(T) GT 2))
                              .. TSUM(T) =E= SUM( I$(ORD(I) LT ORD(T)-1), TSS(T-ORD(I)) * (KER(I+1)-KER(I+2))) ;

CLIMHLPF0(T) $((MODEFLAG < 0.5) AND (ORD(T) EQ 1))
                              .. TSUM(T) =E= SUM( I$(Ord(I) LT 250), TSSMM(I+1)*(KER(I+1)-KER(I+2)) );

CLIMHLPF(T)  $((MODEFLAG < 0.5) AND (ORD(T) GT 1))
                              .. TSUM(T) =E= SUM( I$(Ord(I) LT Ord(T)-1), TSS(T-Ord(I))*(KER(I+1)-KER(I+2)))
                                             + SUM( I$(Ord(I) LT 251), TSSMM(I)*(KER(I+(Ord(T)-1))-KER(I+Ord(T))) );


*** Matrix entries for implicit numerical integration ------------------------------------------------------------------
*
* (Equation A.27 in Kriegler, 2005)

CLIMA11                       .. A11 =E= 1 - 0.5*DELTAT * (LLAM/CAL+KSL/CAL/FLA) ;

CLIMA12                       .. A12 =E= 0.5*DELTAT * BSI * KSL/CAL/FLA ;

CLIMA21                       .. A21 =E= 0.5*DELTAT * KSL/CML/(1-FLA) ;

CLIMA22                       .. A22 =E= 1 - 0.5*DELTAT * (OLAM/CML + BSI*KSL/CML/(1-FLA)) +
                                 (KER("1")-KER("2")) * RML * SHC/CML*SQRT(DELTAT*KDIF/PI) ;

CLIMB11                       .. B11 =E= 1 + 0.5*DELTAT * (LLAM/CAL + KSL/CAL/FLA) ;

CLIMB12                       .. B12 =E= - 0.5*DELTAT * BSI * KSL/CAL/FLA ;

CLIMB21                       .. B21 =E= - 0.5*DELTAT * KSL/CML/(1-FLA) ;

CLIMB22                       .. B22 =E= 1 + 0.5*DELTAT * (OLAM/CML + BSI*KSL/CML/(1-FLA)) +
                                 KER("1") * RML * SHC/CML*SQRT(DELTAT*KDIF/PI) ;

* Two different ways of numerical integration are provided for DOECLIM.
* Runge-Kutta 2nd order is less computationally demanding but would produce oscillations mainly in the land temperature
* when a large change in the radiative forcing is imposed. The oscillations are due to a short turnover time
* for land biosphere (approximately 1.5 months), which cannot be resolved on the time scale of this model.
* However, except for the oscillations, these two numerical integration methods produce identical temperature
* projections in a long run (Kriegler, personal communication). The use of Runge-Kutta 4th order is necessary
* for the past mode, which is characterized by large shifts in the radiative forcing (stemming from volcanic forcing).
* Either the 4th order or the 2nd order can be used for the future mode as long as the system is operated
* with a smooth emission scenario. When Runge-Kutta 2nd order is used, comment out the 4th order equations
* and uncomment the 2nd order equations (both CLIMC** equations and RADNUM** equations).

*# Vanishing correction terms for implicit Runge-Kutta 2nd order

* CLIMC11                       .. C11 =E= 0 ;

* CLIMC12                       .. C12 =E= 0 ;

* CLIMC21                       .. C21 =E= 0 ;

* CLIMC22                       .. C22 =E= 0 ;

*# Non-vanishing correction terms for implicit Runge-Kutta 4th order
*#(added by Elmar Kriegler, 02.02.2006)

CLIMC11                      .. C11 =E= DELTAT**2/12 * ( LLAM*LLAM/CAL**2 + KSL**2/(CAL*FLA)**2 +
                                        2*LLAM*KSL/CAL**2/FLA + BSI*KSL**2/(CAL*CML*FLA*(1-FLA)) );

CLIMC12                      .. C12 =E= DELTAT**2/12 * ( BSI*LLAM*KSL/CAL**2/FLA + BSI*KSL**2/(CAL*FLA)**2 +
                                        BSI*OLAM*KSL/(CAL*CML)/FLA + (BSI*KSL)**2/(CAL*CML*FLA*(1-FLA)) );

CLIMC21                      .. C21 =E= DELTAT**2/12 * ( OLAM*KSL/CML**2/(1-FLA) + BSI*KSL**2/(CML*(1-FLA))**2 +
                                        LLAM*KSL/(CAL*CML)/(1-FLA) + KSL**2/(CAL*CML*FLA*(1-FLA)) );

CLIMC22                      .. C22 =E= DELTAT**2/12 * ( OLAM*OLAM/CML**2 + (BSI*KSL)**2/(CML*(1-FLA))**2 +
                                        2*BSI*OLAM*KSL/CML**2/(1-FLA) + BSI*KSL**2/(CAL*CML*FLA*(1-FLA)) );


*** Model integration --------------------------------------------------------------------------------------------------
*
* See the explanation on Runge-Kutta 2nd order and 4th order above.
* See equation A.27 in Kriegler, 2005 (For Runge-Kutta 2nd order)

*# Implicit Runge-Kutta 2nd order
*
* RADNUMLQ(T) $(ORD(T) < CARD(T)) .. RADNUML(T) =E= 0.5*DELTAT * (RADTTLL(T) + RADTTLL(T+1))/CAL ;
*
* RADNUMOQ(T) $(ORD(T) < CARD(T)) .. RADNUMO(T) =E= 0.5*DELTAT * (RADTTLO(T) + RADTTLO(T+1))/CML ;

*# Implicit Runge-Kutta 4th order (added by Elmar Kriegler, 02.02.2006)

RADNUMLQ(T) $(ORD(T) < CARD(T)) .. RADNUML(T) =E= 0.5*DELTAT * (RADTTLL(T) + RADTTLL(T+1))/CAL
                                    + DELTAT**2/12 * ( (RADTTLL(T+1) - RADTTLL(T))/CAL**2 * (LLAM + KSL/FLA) -
                                    BSI*(RADTTLO(T+1) - RADTTLO(T))/(CML*CAL)/FLA*KSL ) ;

RADNUMOQ(T) $(ORD(T) < CARD(T)) .. RADNUMO(T) =E= 0.5*DELTAT * (RADTTLO(T) + RADTTLO(T+1))/CML
                                    + DELTAT**2/12 * ( (RADTTLO(T+1) - RADTTLO(T))/CML**2*(OLAM + BSI*KSL/(1-FLA)) -
                                    (RADTTLL(T+1) - RADTTLL(T))/(CML*CAL)/(1-FLA)*KSL ) ;

*** Land surface temperature -------------------------------------------------------------------------------------------

CLIMTAL(T) $(ORD(T) > 1)      .. (B11+C11)* TAL(T) + (B12-C12) * TSS(T) =E=
                                 RADNUML(T-1) + (A11+C11) * TAL(T-1) + (A12-C12) * TSS(t-1) ;


*** Sea surface temperature --------------------------------------------------------------------------------------------

CLIMTSS(T) $(ORD(T) > 1)      .. (B22+C22) * TSS(T) + (B21-C21) * TAL(T) =E=
                                 RADNUMO(t-1) + (A21-C21) * TAL(t-1) + (A22+C22) * TSS(t-1) +
                                 RML * SHC/CML*SQRT(DELTAT*KDIF/PI) * TSUM(t-1) ;


*** Global mean surface temperature ------------------------------------------------------------------------------------

CLIMGMT(T)                    .. GMT(T) =E= FLA*TAL(T) + (1-FLA)*TSS(T) ;
* The original corresponding equation in DOECLIM is as follows: GMT(T) =E= FLA*TAL(T) + (1-FLA)*BSI*TSS(T)
* In the framework of ACC2, ocean surface air warming enhancement factor (BSI) should be removed from the equation,
* which was suggested by Elmar Kriegler (personal communication on March 10, 2006).
* This correction gives less enhanced temperature change, resulting in a higher climate sensitivity in the inversion.
* Refer to Appendix p.7 of Krieger (2005) and Section 4.2. of Krieger (2005) for further explanations.


*** Heat budget calculation --------------------------------------------------------------------------------------------
* Ongoing work. The calculation results are not yet in line with Box 3.1, Figure 1 of IPCC AR5 WG1.

FOQ(T) $(ORD(T) > 1)          .. FO(T) =E= (RADTTLO(T)+RADTTLO(T-1))/2/CML - (TSS(T)+TSS(T-1))/(2*CML/OLAM)
                                           - (BSI*(TSS(T)+TSS(T-1)) - (TAL(T)+TAL(T-1)))/(2*(1-FLA)*CML/KSL)
                                           - (TSS(T)-TSS(T-1))/DELTAT;
* Heat flux (W per m2) to the interior ocean from the mixed layer
* Based on equation 2.3.22 of Tanaka (2008, PhD thesis)

HCIOQ(T) $(ORD(T) > 1)        .. HCIO(T) =E= FO(T)*GSA*(1-FLA)*RML*(10**6)*(10**6)*60*60*24*365.25/(10**21);
* Heat flux (ZJ: J*10**21) to the interior ocean from the mixed layer

HCMLQ(T) $(ORD(T) > 1)        .. HCML(T) =E= CML*(TSS(T)-TSS(T-1))*GSA*(1-FLA)*(10**6)*(10**6)*60*60*24*365.25/(10**21);
* Change in the heat content in the mixed-layer ocean-atmosphere system

HCALQ(T) $(ORD(T) > 1)        .. HCAL(T) =E= CAL*(TAL(T)-TAL(T-1))*GSA*FLA*(10**6)*(10**6)*60*60*24*365.25/(10**21);
* Change in the heat content in the land-atmosphere system

HCTLQ(T) $(ORD(T) > 1)        .. HCTL(T) =E= HCIO(T)+HCML(T)+HCAL(T);
* Change in the heat content in the earth system

HCRDQ(T) $(ORD(T) > 1)        .. HCRD(T) =E= (RADTTLL(T)+RADTTLL(T-1))*GSA*FLA*(10**6)*(10**6)*60*60*24*365.25/(10**21)
                                 + (RADTTLO(T)+RADTTLO(T-1))*GSA*(1-FLA)*(10**6)*(10**6)*60*60*24*365.25/(10**21);
* Heat input from the radiative forcing to the system


****** Equation Scaling ================================================================================================

ALKLQ.SCALE(T) = 100000;

AK1PQ.SCALE(T) = 100;

AK1AQ.SCALE(T) = 10;

AK1BQ.SCALE(T) = 1000;

AK1CQ.SCALE(T) = 10;

AK2PQ.SCALE(T) = 100;

AK2AQ.SCALE(T) = 100;

AK2BQ.SCALE(T) = 1000;

AK2CQ.SCALE(T) = 10;

RADH2OSTRQ.SCALE(T) = 0.001;

NORRADQ.SCALE(T) = 0.0001;

PCTRESQ.SCALE(T,NPCR) = 0.01;

SLRGLQ2P.SCALE(T) = 0.0001;

SLRGRQ2P.SCALE(T) = 0.0001;

SLRANQ2P.SCALE(T) = 0.0001;

RADANTLQ.SCALE(T) = 0.01;

RADANTOQ.SCALE(T) = 0.01;

RADSOLLQ.SCALE(T) = 0.0001;

RADSOLOQ.SCALE(T) = 0.0001;

RADTTLLQ.SCALE(T) = 0.01;

RADTTLOQ.SCALE(T) = 0.01;

RADNUMLQ.SCALE(T) = 0.0001;

RADNUMOQ.SCALE(T) = 0.0001;

KHEAT.SCALE = 1000;

CONVOL1.SCALE(I) = 0.1;
