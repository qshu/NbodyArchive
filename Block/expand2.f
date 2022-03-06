      SUBROUTINE EXPAND2(IPAIR,SEMI0)
*
*
*       Expansion (contraction) of KS orbit.
*       ------------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Evaluate square regularized velocity.
      V20 = 0.0
      DO 10 K = 1,4
          V20 = V20 + UDOT(K,IPAIR)**2
   10 CONTINUE
*
*       Form KS coordinate & velocity scaling factors (general point is OK).
      I = N + IPAIR
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      C2 = SQRT(SEMI/SEMI0)
      V2 = 0.5*(BODY(I) + H(IPAIR)*R(IPAIR)*(SEMI/SEMI0))
      IF (V2.LT.0.0) V2 = V20
      C1 = SQRT(V2/V20)
*
*       Re-scale KS variables to new energy (H < 0: constant eccentricity).
      R(IPAIR) = 0.0D0
      TDOT2(IPAIR) = 0.0
      DO 20 K = 1,4
          U(K,IPAIR) = C2*U(K,IPAIR)
          UDOT(K,IPAIR) = C1*UDOT(K,IPAIR)
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
          TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0*U(K,IPAIR)*UDOT(K,IPAIR)
   20 CONTINUE
*
      RETURN
*
      END
