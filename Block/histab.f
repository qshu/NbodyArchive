      SUBROUTINE HISTAB(IPAIR,J,PMIN,RSTAB)
*
*
*       Hierarchical stability condition.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'kscmn6.h'
      REAL*8  XX(3,3),VV(3,3)
      LOGICAL SWAPPED
*
*
*       Define KS & c.m. indices and save IPAIR for RETURN.
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      I = N + IPAIR
      IP = IPAIR
*
*       Evaluate semi-major axis & eccentricity of inner binary.
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      ECC2 = (1.0D0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      ECC = SQRT(ECC2)
*
*       Form stability ratio (Eggleton & Kiseleva, Ap.J. 455, 640, 1995).
*     Q = BODY(I)/BODY(J)
*     Q1 = MAX(BODY(I2)/BODY(I1),BODY(I1)/BODY(I2))
*     Q3 = Q**0.33333
*     Q13 = Q1**0.33333
*       Note sign error in third term of Eq. (2).
*     AR = 1.0 + 3.7/Q3 - 2.2/(1.0 + Q3) + 1.4/Q13*(Q3 - 1.0)/(Q3 + 1.0)
*
*       Save hierarchical stability separation in COMMON variable.
*     PCRIT = AR*SEMI*(1.0D0 + ECC)
*
*       Determine the outer eccentricity.
      RIJ2 = 0.0
      VIJ2 = 0.0
      RDOT = 0.0
      DO 5 K = 1,3
          RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
          VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,J))**2
          RDOT = RDOT + (X(K,I) - X(K,J))*(XDOT(K,I) - XDOT(K,J))
    5 CONTINUE
      RIJ = SQRT(RIJ2)
      A1 = 2.0/RIJ - VIJ2/(BODY(I) + BODY(J))
*
*       Exit on hyperbolic orbit.
      IF (A1.LT.0.0) THEN
          PMIN = 1.0
          RSTAB = 0.0
          GO TO 20
      END IF
*
      SEMI1 = 1.0/A1
      ECC2 = (1.0 - RIJ/SEMI1)**2 + RDOT**2/(SEMI1*(BODY(I) + BODY(J)))
      ECC1 = SQRT(ECC2)
      PMIN = SEMI1*(1.0 - ECC1)
*
*       Evaluate the basic stability condition without fudge factor.
      Q = BODY(J)/BODY(I)
      IF (ECC1.LT.1.0) THEN
          XFAC = (1.0 + Q)*(1.0 + ECC1)/SQRT(1.0 - ECC1)
      ELSE
          XFAC = 40.0*(1.0 + Q)
      END IF
      PCRIT1 = 2.8*XFAC**0.4*SEMI
*
*       Exit if stability value falls outside practical limits.
      IF ((PCRIT1.GT.1.5*PMIN.OR.PCRIT1.LT.0.5*PMIN).AND.J.LE.N) THEN
          RSTAB = PCRIT1
          GO TO 20
      END IF
*
*       Choose the most active triple in case of two binaries.
      JJ = J
      SWAPPED = .FALSE.
      IF (J.GT.N) THEN
          JP = J - N
          IF (LKSINT(JP)) THEN
              HJ = HX(JP)
              RJ = RX(JP)
              TDOT2J = TDOT2X(JP)
          ELSE
              HJ = H(JP)
              RJ = R(JP)
              TDOT2J = TDOT2(JP)
          END IF
          SEMI2 = -0.5D0*BODY(J)/HJ
*       Adopt 10% fudge factor with linear dependence on smallest ratio.
          YFAC = 1.0 + 0.1*MIN(SEMI2/SEMI,SEMI/SEMI2)
          IF (SEMI2.GT.SEMI) THEN
              ECC2 = (1.0 - RJ/SEMI2)**2 + TDOT2J**2/(BODY(J)*SEMI2)
              ECC = SQRT(ECC2)
              SEMIZ = SEMI2
              SEMI2 = SEMI
              SEMI = SEMIZ
              IP = JP
              I1 = 2*IP - 1
              I2 = I1 + 1
              JJ = I
              I = J
              IF (LKSINT(IPAIR).AND.LKSINT(IP)) SWAPPED = .TRUE.
          END IF
      ELSE
          YFAC = 1.0
      END IF
*
*       Resolve weakly perturbed binary (prevent X(K,I1) = X(K,I2)).
      IF (SWAPPED) THEN
          IF (GAMMAX(IP).LT.GMIN.OR.XKS(1,I1).EQ.XKS(1,I2)) THEN
              CALL RESOLVX(IP)
          END IF
      ELSE
          IF (GAMMA(IP).LT.GMIN.OR.X(1,I1).EQ.X(1,I2)) THEN
              CALL RESOLV(IP,1)
          END IF
      END IF
*
*       Determine inclination for triple configuration (NB! radians).
      DO 10 K = 1,3
          IF (SWAPPED) THEN
              XX(K,1) = XKS(K,I1)
              XX(K,2) = XKS(K,I2)
              VV(K,1) = XDOTKS(K,I1)
              VV(K,2) = XDOTKS(K,I2)
          ELSE
              XX(K,1) = X(K,I1)
              XX(K,2) = X(K,I2)
              VV(K,1) = XDOT(K,I1)
              VV(K,2) = XDOT(K,I2)
          END IF
          XX(K,3) = X(K,JJ)
          VV(K,3) = XDOT(K,JJ)
  10  CONTINUE
      CALL INCLIN(XX,VV,X(1,I),XDOT(1,I),ANGLE)
*
*       Employ the improved stability criterion for doubtful cases.
      RSTAB = QSTAB(ECC,ECC1,ANGLE,BODY(I1),BODY(I2),BODY(JJ))*SEMI
*     RSTAB = YFAC*RSTAB
*
   20 RETURN
*
      END
