********** ACC2 version 4.2: ACC2_switchboard_common.inc ***************************************************************
****** Switchboard (past and future modes) =============================================================================
*
*
* If you wish to perform the standard inverse calculation,
* under the parameters and data sections, assign '0' to Q10TMPSWCH, TNPPSENSSWCH, RAD2XCO2SWCH, SCLRADSWCH, KAPPASWCH,
* RADMISAR1PROSWCH, and T2MAR1PROSWCH, and set all the other switches to '1'.
* Furthermore, set LNDCO2TFNPPSWCH, LNDCO2TFQ10SWCH, RADMISAR1SWCH, and T2MAR1SWCH to be '0' and the rest to be '1'.
* If you wish to perform the forward calculation (without any inverse estimation),
* assign '0' to all the switches under PARAMETERS and DATA sections.
*
* Assignment of '1' to a switch under the Parameters section indicates the associated residual term enters
* in the cost function. When a switch is set to '0', the corresponding parameter is fixed at the prior mean value.
* Similarly, the assignment of '1' for a switch under the Data section indicates that the associated
* residual term appears in the cost function. When a switch is set to '0', the corresponding measurements are
* not used in the inverse calculation.
*
* Upon any change in the switch assignments, one should be aware of the degrees of freedom (indicated in brackets).
* The degree of freedom for parameters must not exceed that for data.
*
*
*** Parameters ---------------------------------------------------------------------------------------------------------
*
         CO2FFSWCH =        1;
*        df(251): CO2 emission due to fossil fuel combustion between 1750 and 2000
*
         CO2LUSWCH =        1;
*        df(251): CO2 emission due to land use change between 1750 and 2000
*
         ANTCH4N2OSWCH =    1;
*        df(502): Anthropogenic emissions of CH4 and N2O between 1750 and 2000
*
         RADMISSWCH =       0;
*        df(251): Missing forcing between 1750 and 2000
*
         TMIX1750SWCH =     1;
*        df(1): Preindustrial mixed layer temperature
*
         DMPTMIXSWCH =      1;
*        df(1): Atmosphere-mixed layer temperature scaling factor
*
         BETASWCH =         0;
*        df(1): Beta factor for CO2 fertilization
*
         Q10BASSWCH =       0;
*        df(1): Q10 for heterotrophic respiration
*
         Q10TMPSWCH =       0;
*        df(1): Temperature sensitivity to Q10
*        '0' is assigned by default. See the discussion on LNDCO2TFQ10SWCH below.
*
         TNPPSENSSWCH =     0;
*        df(1): Temperature sensitivity to NPP
*        '0' is assigned by default. Further discussion is in Tanaka (2008) p.38-39.
*
         CO2OCNUPPRESWCH =  1;
*        df(1): Preindustrial ocean CO2 uptake between 1750 and 2000 (constant flux)
*
         CO2LNDUPPRESWCH =  1;
*        df(1): Preindustrial land CO2 uptake between 1750 and 2000 (constant flux)
*
         NATCH4N2OSWCH =    1;
*        df(2): Natural emissions of CH4 and N2O between 1750 and 2000 (constant flux)
*
         TAUCH4N2OSWCH =    1;
*        df(2): CH4 Lifetime with respect to OH depletion and N2O lifetime
*
         RAD2XCO2SWCH =     0;
*        df(1): CO2 doubling forcing
*        '0' is assigned by default as CO2 doubling forcing cannot be tightly constrained
*        together with climate sensitivity in the ACC2 inversion.
*
         SCLRADSWCH =       1;
*        df(1): Forcing scaling factor
*        '0' is assigned by default. With this switch set to '1' and RADMISSWCH set to '0', only the time evolution of
*        the total aerosol forcing is adjusted in the inversion (forcing scaling experiment).
*
         SENSSWCH =         0;
*        df(1): Climate sensitivity
*
         KAPPASWCH =        0;
*        df(1): Vertical ocean diffusivity
*        Kappa is by default fixed at Kriegler's best estimate (KAPPA =0.55 cm^2/s) in the inversion
*        because Kappa cannot be simultaneously constrained with climate sensitivity
*        in the absence of the ocean heat distribution data in the inversion.
*        If the switch assignment is '1,', the ocean vertical diffusivity can be treated as a parameter
*        in the inverse calculation (DOECLIMSWCH must be set to '1' here).
*        See Tanaka (2008) p.109-110.
*
         RADMISAR1PROSWCH = 0;
*        df(1): Autoregressive propagator for missing forcing residuals
*        '0' is assigned by default. See the discussion for RADMISAR1SWCH below.
*
         T2MAR1PROSWCH =    0;
*        df(1): Autoregressive propagator for global-mean surface air temperature change residuals
*        '0' is assigned by default. See the discussion for T2MAR1SWCH below.
*
*
*** Data ---------------------------------------------------------------------------------------------------------------
*
         CO2UPANTSWCH =     1;
*        df(498): Ocean and land CO2 uptake between 1751 and 1999
*
         CONCO2SWCH =       1;
*        df(250): Atmospheric CO2 concentration between 1751 and 2000
*
         CONCH4N2OSWCH =    1;
*        df(500): Atmospheric CH4 and N2O concentrations between 1751 and 2000
*
         T2MSWCH =          1;
*        df(251): Global-mean surface air temperature change between 1750 and 2000
*
*
*** Choice of climate model --------------------------------------------------------------------------------------------
*
         DOECLIMSWCH =      1;
*        If the switch assignment is '1,' DOECLIM is used for calculation from radiative forcing to surface temperature.
*        If the switch assignment is '0,' physical interpretation of the impulse response functions
*        (Hooss, 2001; Kriegler, unpublished) is used instead.
*        The sea level rise projection is calculated from the physical interpretation of the impulse response functions
*        regardless of the assignment for this switch.
*
*
*** Feedbacks ----------------------------------------------------------------------------------------------------------
*
         DYNCARBSWCH =      1;
*        A switch assignment of '1' indicates that the thermodynamic equilibria for carbonate species are calculated,
*        indicating that the saturation of the ocean CO2 uptake under rising atmospheric CO2 concentration can be
*        dynamically reproduced.
*        When the switch assignment is '0,' the thermodynamic equilibria are calculated based on a prescribed Revelle
*        factor (see equation REVELLEQ).
*
         OCNCO2TFSWCH =     0;
*        A switch assignment of '1' indicates that the ocean CO2 uptake is sensitive to the mixed-layer temperature.
*        The associated parameters (TMIX1750 and DMPTMIX) can be estimated in the inverse calculation by assigning the
*        respective switches (TMIX1750SWCH and DMPTMIXSWCH) to '1'.
*        When OCNCO2TFSWCH is '0', thermodynamic constants are fixed at the values obtained for the preindustrial
*        mixed-layer temperature. Thus, ocean CO2 uptake is not sensitive to temperature change.
*
         LNDCO2TFNPPSWCH =  0;
*        A switch assignment of '1' indicates that terrestrial NPP is sensitive to the surface air temperature change
*        over land (when DOECLIM is used) or the global-mean surface temperature change (when impulse response functions
*        are used). The associated parameter (TNPPSENS) can be estimated in the inverse calculation with TNPPSENSSWCH.
*        If LNDCO2TFNPPSWCH is assigned to be '0', terrestrial NPP is kept at its preindustrial estimate.
*        '0' is assigned by default. Further discussion is in Tanaka (2008) p.38-39.
*
         LNDCO2TFRESSWCH =  1;
*        A switch assignment of '1' indicates that the heterotrophic respiration is sensitive to the surface air
*        temperature change over land (when DOECLIM is used) or the global-mean surface temperature change (when impulse
*        response functions are used). The associated parameter (Q10BAS (and Q10TMP)) can be estimated in the inverse
*        calculation by assigning the respective switches (Q10BASSWCH and Q10TMPSWCH) to '1'.
*        When LNDCO2TFRESSWCH is assigned to be '0', heterotrophic respiration is not sensitive to the temperature
*        change. Heterotrophic respiration is then reproduced as in the impulse response function of Joos et al. (1996).
*
         LNDCO2TFQ10SWCH =  0;
*        A switch assignment of '1' indicates that Q10 itself is sensitive to the surface air temperature change over
*        land (when DOECLIM is used) or the global-mean surface temperature change (when impulse response functions are
*        used). The associated sensitivity parameter (Q10TMP) can be estimated in the inverse calculation by assigning
*        the corresponding switch Q10TMPSWCH to '1'.
*        When LNDCO2TFQ10SWCH is assigned to be '0', Q10 is kept constant irrespective of the temperature.
*        '0' is assigned by default as the introduction of the temperature sensitivity for Q10 would lead to an
*        unusually low Q10 in the latter half of the 21st century.
*        Further discussion is in Tanaka (2008) p.38-39.
*
*
*** Assumptions for past mode ------------------------------------------------------------------------------------------
*
         CO2UPPRESWCH =     1;
*        A switch assignment of '1' indicates that the preindustrial ocean and land CO2 uptake between 1750 and 2000
*        (constant fluxes) is added to the carbon cycle. The magnitudes of such fluxes can be estimated in the inverse
*        calculation with respective switch assignments (CO2OCNUPPRESWCH and CO2LNDUPPRESWCH).
*
         CO2LUSGMSWCH =     1;
*        A switch assignment of '1' indicates that the uncertainty range of the CO2 emission due to land use change
*        is assumed to be larger than the literature estimate (Houghton, 2003) by a factor designated in parameter
*        'CO2LUSGMFAC' (= 2 by default). Rationales for such a modification are stated in Tanaka (2008) p.106-107.
*        If the switch is assigned to be '0', there is no modification for the Houghton's uncertainty range.
*
         CONCO2SGMSWCH =    1;
*        A switch assignment of '1' indicates that the uncertainty range of the atmospheric CO2 concentration
*        is assumed to be larger after 1959 than in the literature (Keeling, 2005) by a factor designated in
*        parameter 'CONCO2SGMFAC' (= 4 by default). Otherwise, the switch is assigned to be '0.'
*        See Tanaka (2008) p.104-105.
*
         ANTCH4N2OSGMSWCH = 1;
*        A switch assignment of '1' indicates that the prior uncertainty ranges of the emissions of CH4 and N2O are
*        assumed to be smaller than in the literature (van Aardenne et al, 2001) by a factor designated in parameter
*        'ANTCH4N2OSGMFAC' (= 0.5 by default). Such a modification of the uncertainty ranges reflects the fact
*        that the prior uncertainties in van Aardenne et al. (2001) are assessed only at a qualitative level
*        and also the fact that these literature ranges are disproportionally larger than those for the CO2 emission.
*        Otherwise, the switch is assigned to be '0.'
*
         CONCH4SGMSWCH =    1;
*        Treatment for the prior uncertainty range of CH4 concentration similar to CONCO2SGMSWCH (see above).
*        Parameter 'CONCH4SGMFAC' is set to 4 by default.
*
         CONN2OSGMSWCH =    1;
*        Treatment for the prior uncertainty range of N2O concentration similar to CONCO2SGMSWCH (see above).
*        Parameter 'CONN2OSGMFAC' is set to 4 by default.
*
         ENSOCO2SWCH =      1;
*        A switch assignment of '1' indicates that the inversion takes account of the ENSO influence to the variation
*        in the atmospheric CO2 concentration. In other words, the ENSO effect is added to the
*        model-calculated natural CO2 concentration changes between 1960 and 2000 by using a regression line as a
*        function of NINO3 (with 7-month time lag) computed dynamically in the inversion.
*        Otherwise, the switch is assigned to be '0'. See Tanaka (2008) Section 3.5.1.
*        Note that, when the optimization period is set to end before 1985 (that is, TENDINV < 196), the ENSO effect
*        on the carbon cycle is not considered regardless of the switch assignment because of the lack of data points
*        to derive a regression line.
*
         ENSOT2MSWCH =      1;
*        A switch assignment of '1' indicates that the ENSO influence to the surface temperature change
*        is included in the model-calculated temperature change between 1930 and 2000. In contrast to the case of
*        ENSOCO2SWCH, a fixed linear relationship between SOI and the surface temperature change is prescribed.
*        Otherwise, the switch is assigned to be '0'. See Tanaka (2008) Section 3.5.1.
*
         RADMISAR1SWCH =    0;
*        A switch assignment of '1' indicates that the residuals of missing forcing are corrected
*        by using a first-order autoregressive model. The autoregressive model is a stochastic model exhibiting a red
*        power spectrum, which explains the high frequency variability in the missing forcing.
*        Otherwise, the switch is assigned to be '0'.
*        This switch is '0' in the standard inverse calculation because it is not clear how an autocorrelation model
*        fit into the probability densities -- the autocorrelation model is not compatible with the probabilistic
*        inverse estimation theory (Tanaka, 2008, p.97-101).
*
         T2MAR1SWCH =       0;
*        A switch assignment of '1' indicates that the residuals of global-mean surface air temperature are corrected
*        by using a first-order autoregressive model. The autoregressive model is a stochastic model exhibiting a red
*        power spectrum, which explains the high frequency variability in the temperature measurements.
*        The implementation follows Kriegler (2005).
*        Otherwise, the switch is assigned to be '0'.
*        This switch is '0' in the standard inverse calculation because it is not clear how an autocorrelation model
*        fit into the probability densities -- the autocorrelation model is not compatible with the probabilistic
*        inverse estimation theory (Tanaka, 2008, p.97-101).
*
         VOLSWCH =          1;
*        A switch assignment of '1' indicates that the uncertainty ranges for atmospheric CO2, CH4, and N2O
*        concentrations, missing radiative forcing, and surface temperature change are assumed to be larger than
*        original estimates by a factor designated in parameter 'VOLSGMFAC' (= 4 by default) when the volcanic forcing
*        is stronger than -0.5 W/m2.
*        See Tanaka (2008) Section 3.5.2.
*
*
*** Assumptions for future mode ----------------------------------------------------------------------------------------
*
         MISTRNSWCH =       0;
*        A switch assignment of '1' indicates that the average ratio of the posterior missing forcing to the total
*        aerosol forcing during 1950-2000 is used to scale the total aerosol forcing of the future mode. The averages
*        of the missing forcing and the total aerosol forcing during 1950-2000 are transferred from the past mode to
*        the future mode with the variables RADMISAVE and RADAERAVE, respectively. The missing forcing during 1950-2000
*        can be interpreted as the aerosol forcing uncertainty. The scaling of the future aerosol forcing aims at
*        ensuring consistency in the estimates of uncertainties in the past and future modes.
*
         SOLSWCH =          1;
*        A switch assignment of '1' indicates that solar forcing is added by making the following assumptions:
*        1) The average solar forcing over 20th century is used between 2031 and 2100 (or thereafter).
*        2) The solar forcing prior to 2030 are linearly interpolated with the solar forcing in 2000 of the past mode.
