      SUBROUTINE KSPERI(IPAIR)
!VAMPINST include
*
*
*       Pericentre KS variables.
*       ------------------------
*
      INCLUDE 'common6.h'
*
*
*       Set c.m. index & semi-major axis.
!VAMPINST subroutine_start
      ICM = N + IPAIR
      SEMI = -0.5D0*BODY(ICM)/H(IPAIR)
*
*       Define auxiliary quantities and obtain the eccentricity.
      ZETA = 1.0 - R(IPAIR)/SEMI
      PSI = TDOT2(IPAIR)/SQRT(BODY(ICM))
      ECC = SQRT(ZETA**2 + PSI**2/SEMI)
*
*       Avoid nearly circular orbits (undefined pericentre).
      IF (ECC.LT.0.0001) GO TO 20
*
*       Distinguish between near-collision, elliptic & hyperbolic case.
      Q = (PSI/ECC)**2
      IF (ZETA.GT.0.0D0.AND.ABS(Q/SEMI).LT.0.1) THEN
*       Expansion of THETA/(SIN(THETA)*COS(THETA)) - 1 for Q = SIN(THETA)**2.
          SN1 = 1.0
          S = 0.0D0
          Q = Q/SEMI
          DO 1 IT = 1,10
              SN1 = Q*SN1*FLOAT(2*IT)/FLOAT(2*IT + 1)
              S = S + SN1
    1     CONTINUE
          S = S + SN1*Q*0.9D0/(1.0 - Q)
          DT = (R(IPAIR)*ZETA - PSI**2 + ZETA*S*SEMI)*PSI/
     &                                          (ECC**2*SQRT(BODY(ICM)))
      ELSE IF (SEMI.GT.0.0D0) THEN
*       Determine the eccentric anomaly with respect to pericentre (0,PI).
          THETA = ATAN2(ABS(PSI)/SQRT(SEMI),ZETA)
*
*       Obtain pericentre time interval from Kepler's equation.
          DT = SEMI*SQRT(SEMI/BODY(ICM))*(THETA - ABS(PSI)/SQRT(SEMI))
      ELSE IF (SEMI.LT.0.0D0) THEN
*       Hyperbolic case.
          A1 = PSI/(ECC*SQRT(ABS(SEMI)))
*       Use EXP(F) = SINH(F) + COSH(F) to obtain the eccentric anomaly THETA.
          A2 = ABS(A1) + SQRT(A1**2 + 1.0D0)
          THETA = LOG(A2)
          IF (A1.LT.0.0D0) THETA = -THETA
          A0 = ABS(SEMI)
          DT = A0*SQRT(A0/BODY(ICM))*(ABS(PSI)/SQRT(A0) - THETA)
      END IF
*
*       Specify transformation coefficients (Mikkola's procedure).
      XC = SQRT(0.5D0 + 0.5D0*ZETA/ECC)
      YS = PSI/(ECC*XC*SQRT(BODY(ICM)))
      ZZ = BODY(ICM)/(4.0*SEMI)
      R(IPAIR) = 0.0D0
      TDOT2(IPAIR) = 0.0D0
*
*       Generate analytical solutions for U & UDOT using old U0 & UDOT.
      DO 10 K = 1,4
          U(K,IPAIR) = U0(K,IPAIR)*XC - UDOT(K,IPAIR)*YS
          UDOT(K,IPAIR) = U0(K,IPAIR)*YS*ZZ + UDOT(K,IPAIR)*XC
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
          TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0D0*U(K,IPAIR)*UDOT(K,IPAIR)
   10 CONTINUE
*
*       Do not allow R' < 0 for repeated pericentre in routine KSINT.
      IF (TDOT2(IPAIR).LT.0.0D0) THEN
          TDOT2(IPAIR) = 0.0D0
      END IF
*
*       Set pericentre time and predict c.m. coordinates & velocities.
      TIME = TIME - DT
      CALL XVPRED(ICM,0)
*
!VAMPINST return
20        RETURN
      RETURN
*
!VAMPINST subroutine_end
      END
