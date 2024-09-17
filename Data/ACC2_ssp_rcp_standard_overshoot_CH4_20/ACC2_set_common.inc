********** ACC2 version 4.2: ACC2_set_common.inc ***********************************************************************
********** Set declarations and assignments (past and future modes) ====================================================

SET FOA        Forcing agents and ozone depleting substances group A
               /CO2,
                CH4, N2O,
                rOH,
                CF4, C2F6, C4F10, HFC23, HFC32, HFC4310mee, HFC125, HFC134a, HFC143a, HFC152a, HFC227ea, HFC236fa,
                 HFC245ca,
                CFC11, CFC12, CFC113, CFC114, CFC115, CCl4, CH3CCl3, Halon1211, Halon1301, Halon2402, HCFC22, HCFC141b,
                 HCFC142b, HCFC123,
                CH3Cl, CH3Br,
                EESC,
                SF6,
                NOx, CO, VOC,
                O3TRP,
                O3STR,
                SO2,
                CRBAER, AERIND,
                H2OSTR
                /;

SET FOB(FOA)   Forcing agents and ozone depleting substances group B
               /CO2,
                CH4, N2O,
                rOH,
                CF4, C2F6, C4F10, HFC23, HFC32, HFC4310mee, HFC125, HFC134a, HFC143a, HFC152a, HFC227ea, HFC236fa,
                 HFC245ca,
                CFC11, CFC12, CFC113, CFC114, CFC115, CCl4, CH3CCl3, Halon1211,Halon1301, Halon2402, HCFC22, HCFC141b,
                 HCFC142b, HCFC123,
                CH3Cl, CH3Br,
                EESC
                SF6,
                NOx, CO, VOC,
                O3TRP,
                SO2/;

SET FOC(FOB)   Forcing agents and ozone depleting substances group C
               /CH4, N2O,
                CF4, C2F6, C4F10, HFC23, HFC32, HFC4310mee, HFC125, HFC134a, HFC143a, HFC152a, HFC227ea, HFC236fa,
                 HFC245ca,
                CFC11, CFC12, CFC113, CFC114, CFC115, CCl4, CH3CCl3, Halon1211, Halon1301, Halon2402, HCFC22, HCFC141b,
                 HCFC142b, HCFC123,
                CH3Cl, CH3Br,
                SF6,
                NOx, CO, VOC,
                SO2/;

SET FOD(FOB)   Forcing agents group D
               /CO2,
                CH4, N2O,
                NOx, CO, VOC,
                SO2/;

SET FOE(FOD)   Forcing agents group E
               /CH4, N2O,
                NOx, CO, VOC,
                SO2/;

SET FOF(FOB)   Forcing agents group F
               /CO2,
                CH4, N2O,
                rOH/;

SET FOG(FOE)   Forcing agents group G
               /CH4, N2O/;

SET FOH(FOB)   Forcing agents group H
               /CO2, CH4, N2O/;

SET HAL(FOC)   Halocarbons
               /CF4, C2F6, C4F10, HFC23, HFC32, HFC4310mee, HFC125, HFC134a, HFC143a, HFC152a, HFC227ea, HFC236fa,
                 HFC245ca,
                CFC11, CFC12, CFC113, CFC114, CFC115, CCl4, CH3CCl3, Halon1211, Halon1301, Halon2402, HCFC22, HCFC141b,
                 HCFC142b, HCFC123,
                CH3Cl, CH3Br/;

SET HALR(FOC)  Halocarbons except for CH3Cl and CH3Br
               /CF4, C2F6, C4F10, HFC23, HFC32, HFC4310mee, HFC125, HFC134a, HFC143a, HFC152a, HFC227ea, HFC236fa,
                 HFC245ca,
                CFC11, CFC12, CFC113, CFC114, CFC115, CCl4, CH3CCl3, Halon1211, Halon1301, Halon2402, HCFC22, HCFC141b,
                 HCFC142b, HCFC123/;

SET ODS(FOC)   Ozone depleting substances
               /CFC11, CFC12, CFC113, CFC114, CFC115, CCl4, CH3CCl3, Halon1211, Halon1301, Halon2402, HCFC22, HCFC141b,
                 HCFC142b, HCFC123,
                CH3Cl, CH3Br/;

SET ODSR(FOC)  Ozone depleting substances except for CH3Cl and CH3Br
               /CFC11, CFC12, CFC113, CFC114, CFC115, CCl4, CH3CCl3, Halon1211, Halon1301, Halon2402, HCFC22, HCFC141b,
                 HCFC142b, HCFC123/;

SET I          Time step index from 1750 to 3000 for integral convolution for DOECLIM
               /1*1251/;
* If you run the model beyond the year 3000, the member for set I needs to be expanded accordingly.

SET NOC        Designation of reservoirs in ocean carbon cycle model
               /1*4/;

SET NTR        Designation of reservoirs in land carbon cycle model
               /1*4/;

SET NPCR       Designation of principal components
               /1,2/;

SET LAGCO2     Time lag for natural atmospheric CO2 concentration change behind the NINO3 index
               /0*10/;

SET LAGT2M     Time lag for surface air temperature behind the SOI index
               /0*10/;

SET ENS        Ensemble members of CESM and CCSM runs
               /1*8/;


ALIAS (NPCR,NPCRSEN);
* The above alias statement declares that NPCR and NPCRSEN contain the identical set members.
