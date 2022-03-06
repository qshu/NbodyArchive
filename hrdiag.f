      SUBROUTINE HRDIAG( MASS, AJ, TM, TN, LUMS, TSCLS, R, TE, KW)
*
*
*       H-R diagram for population I stars.
*       -----------------------------------
*
*       Computes the new mass, radius, effective temperature & stellar type.
*       Input (MASS, AJ, TM, TN, LUMS & TSCLS) supplied by routine STAR.
*       Ref: P.P. Eggleton, M.J. Fitchett & C.A. Tout (1989) Ap.J. 347, 998.
*
      REAL  MASS, TSCLS(6), LUMS(6), TM, TN, AM, CL1, CL2, TE, AH, AJ,
     &        AY, VMF, R, FF, R1, R2, R3, RZ, RT, RG, LUM, TGS, AG, LG
*
*
*       ---------------------------------------------------------------------
*       MASS    Stellar mass in solar units (input: old, output: new value).
*       AJ      Current age in Myr.
*       TM      Main sequence time.
*       TN      Nuclear burning time.
*       LUMS    Characteristic luminosity.
*       TSCLS   Time scale for different stages.
*       R       Stellar radius in solar units.
*       TE      Effective temperature (suppressed).
*       KW      Classification type (1 - 8).
*       ---------------------------------------------------------------------
*
*
      TGS = TSCLS(5)
      LG = LUMS(3)
      AM = ALOG10( MASS )
      AH = TM + TSCLS(1)
      AG = AH + TSCLS(2)
      AY = AG + TSCLS(3)
*
      IF (MASS.GT.8.0) THEN
          VMF = 5.0
      ELSE
          VMF = MIN(1.4, 0.4 + MASS/8.0)
      END IF
*
      IF (AJ.GT.TN) THEN
          LUM = 40.0 / (AJ - TN)**1.4
          R = .01
          KW = 7
          IF ( MASS .GT.8.0) KW = 8
*       White dwarf / neutron star.
          MASS = MIN(VMF, 1.4)
*
      ELSE IF ( MASS .GT.25.0.AND.AJ.GT.TM) THEN
          LUM = 1.0E+05
          R = 5.0
          MASS = 10.0
          KW = 6
*       Wolf-Rayet star.
*
      ELSE IF (AJ.GE.AY) THEN
          LUM = LG*(TGS / (TN + TSCLS(6) - AJ))**1.17
          MASS = MASS - (MASS - VMF)* (LUM / LUMS(6)) **0.8
          R = (0.25* LUM **0.4 + 0.8* LUM **0.67) / MASS **0.27
          KW = 5
*       Asymptotic red supergiant.
*
      ELSE IF (AJ.GE.AG) THEN
          LUM = LUMS(5)
          MASS = MASS - (MASS - VMF)* (LUMS(4) / LUMS(6)) **0.8
          R = (0.25* LUM **0.4 + 0.8* LUM **0.67) / MASS **0.27
          FF = (AJ - AG) / TSCLS(3)
          IF (R.GE.25.0) THEN
              R = R*(25.0/R)**FF
          ELSE
              R = R*(1.0 - 0.1*FF)
          END IF
          KW = 4
*       Core helium burning star.
*
      ELSE IF (AJ.LT.AH) THEN
          IF ( MASS .GT.1.334) THEN
              R1 = .1509 + .1709*AM
              R2 = .06656 - .4805*AM
              R3 = .007395 + .5083*AM
              RZ = (.7388* MASS **1.679 - 1.968* MASS **2.887) /
     &                                       (1.0 - 1.821* MASS **2.337)
          ELSE
              R1 = .08353 + .0565*AM
              R2 = .01291 + .2226*AM
              R3 = .1151 + .06267*AM
              RZ =  MASS **1.25*(.1148 + .8604* MASS * MASS ) /
     &                                           (.04651 + MASS * MASS )
          END IF
*
          IF (AJ.LT.TM) THEN
              IF ( MASS .GT.1.334) THEN
                  CL1 = .092088 + .059338*AM
                  CL2 = -.05713 + AM*(.3756 - .1744*AM)
              ELSE
                  CL1 = .2594 + .1348*AM
                  CL2 = .144 - .8333*AM
              END IF
              FF = AJ / TM
              LUM = LUMS(1) *10**(FF*(FF*CL2 + CL1))
              R = RZ*10**(FF*(FF*(FF*R3 + R2) + R1))
              KW = 1
*       Main sequence star.
*
          ELSE
              RT = RZ*10**(R1 + R2 + R3)
              RG = (0.25*LG**0.4 + 0.8*LG**0.67) / MASS **0.27
              FF = (AJ - TM) / TSCLS(1)
              R = RT*(RG/RT)**FF
              LUM = LUMS(2)*(LG/LUMS(2))**FF
              KW = 2
*       Hertzsprung-gap star.
          END IF
*
      ELSE
          LUM = LG*(TGS / (AH + TGS - AJ))**1.17
          MASS = MASS - (MASS - VMF)* (LUM / LUMS(6)) **0.8
          R = (0.25* LUM **0.4 + 0.8* LUM **0.67) / MASS **0.27
          KW = 3
*       Red giant.
      END IF
*
*     TE = (1130* LUM / (R*R))**0.25
      TE = 0.0
*
      RETURN
*
      END
