      SUBROUTINE HISTAB(IPAIR,J,PMIN,RSTAB)
*
*
*       Hierarchical stability condition.
*       ---------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Define KS & c.m. indices.
      I1 = 2*IPAIR - 1
*     I2 = I1 + 1
      I = N + IPAIR
*
*       Set semi-major axis & eccentricity of inner binary.
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
      SEMI1 = 1.0/A1
      ECC2 = (1.0 - RIJ/SEMI1)**2 + RDOT**2/(SEMI1*(BODY(I) + BODY(J)))
      ECC1 = SQRT(ECC2)
      PMIN = SEMI1*(1.0 - ECC1)
*
*       Evaluate the stability formula without fudge factor.
      Q = BODY(J)/BODY(I)
      IF (ECC1.LT.1.0) THEN
*       ZFAC = (1.0 + ECC1)/(1.0 - ECC1)**0.182*(1.0 + Q)/(1.0 + ECC)**3
          XFAC = (1.0 + Q)*(1.0 + ECC1)/SQRT(1.0 - ECC1)
      ELSE
          XFAC = 1.0 + Q
      END IF
*     PCRIT = 2.6*ZFAC**0.355*SEMI*(1.0 + ECC)
      FE  = 1.0
      PCRIT = 2.8*FE*XFAC**0.4*SEMI
      RSTAB = PCRIT
*
      RETURN
*
      END
