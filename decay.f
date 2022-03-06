      SUBROUTINE DECAY(IPAIR,SEMI)
*
*
*       Re-scaling of KS variables to low eccentricity.
*       -----------------------------------------------
*
      INCLUDE 'common4.h'
*
*
*       Define new apo or pericentre for small eccentricity.
      IF (R(IPAIR).GT.SEMI) THEN
          RX = 1.002*SEMI
      ELSE
          RX = 0.998*SEMI
      END IF
*
*       Form KS coordinate & velocity scaling factors assuming J = const.
      C2 = SQRT(RX/R(IPAIR))
      C1 = 1.0/C2
*
*       Re-scale KS variables to eccentricity 0.002 at apo or pericentre.
      R(IPAIR) = 0.0D0
*     TDOT2(IPAIR) = 0.0D0
      DO 20 K = 1,4
          U(K,IPAIR) = C2*U(K,IPAIR)
          UDOT(K,IPAIR) = C1*UDOT(K,IPAIR)
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
*         TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0*U(K,IPAIR)*UDOT(K,IPAIR)
   20 CONTINUE
*
      RETURN
*
      END
