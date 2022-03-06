      SUBROUTINE KSLIST(IPAIR)
*
*
*       KS perturber selection.
*       -----------------------
*
      INCLUDE 'common6.h'
      REAL  W1,W2,W3,RSEP2
*
*
*       Set component & c.m. index and semi-major axis.
      I1 = 2*IPAIR - 1
      I = N + IPAIR
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
*
*       Prepare standard perturber selection using apocentre distance.
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      RAP = SEMI*(1.0 + SQRT(ECC2))
      RCRIT2 = 2.0*RAP**2/BODY(I)
      RCRIT3 = RCRIT2*RAP/GMIN
*       Base fast search on maximum binary mass (2*BODY1).
      RCRIT2 = 2.0*RCRIT2*BODY1*CMSEP2
*
*       Select new perturbers from c.m. neighbour list.
      NNB1 = 1
      NNB2 = LIST(1,I) + 1
      DO 10 L = 2,NNB2
          J = LIST(L,I)
          W1 = X(1,J) - X(1,I)
          W2 = X(2,J) - X(2,I)
          W3 = X(3,J) - X(3,I)
          RSEP2 = W1*W1 + W2*W2 + W3*W3
*       Also include any merged c.m. or chain c.m. bodies in the fast test.
          IF (RSEP2.LT.RCRIT2.OR.NAME(J).LE.0) THEN
              RIJ3 = RSEP2*SQRT(RSEP2)
*       Estimate unperturbed distance from tidal limit approximation.
              IF (RIJ3.LT.BODY(J)*RCRIT3) THEN
                  NNB1 = NNB1 + 1
                  LIST(NNB1,I1) = J
              END IF
          END IF
   10 CONTINUE
*
*       Check case of no perturbers (dual purpose).
      IF (NNB1.EQ.1) THEN
*       Retain the previous perturber after partial reflection.
          IF (TDOT2(IPAIR).GT.0.0D0) THEN
              NNB1 = 2
          ELSE IF (KZ(27).EQ.0) THEN
*       Specify one unperturbed period at apocentre (standard case).
              STEP(I1) = TWOPI*SEMI*SQRT(SEMI/BODY(I))
          ELSE
*       Maintain perturbed motion during tidal capture process.
              RP = SEMI*(1.0D0 - SQRT(ECC2))
              RT = 2.7*MAX(RADIUS(I1),RADIUS(I1+1))
              IF (SEMI.GT.RSYNC.AND.RP.LT.RT*(1.0 + GAMMA(IPAIR))) THEN
                  NNB1 = 2
              ELSE
                  STEP(I1) = TWOPI*SEMI*SQRT(SEMI/BODY(I))
              END IF
          END IF
      END IF
*
*       Save perturber membership.
      LIST(1,I1) = NNB1 - 1
*
      RETURN
*
      END
