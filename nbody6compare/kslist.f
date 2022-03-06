      SUBROUTINE KSLIST(IPAIR)
*
*
*       KS perturber selection.
*       -----------------------
*
      INCLUDE 'common6.h'
*
*
*       Set component & c.m. index and form semi-major axis & eccentricity.
      I1 = 2*IPAIR - 1
      I = N + IPAIR
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      EB = -0.5*BODY(I1)*BODY(I1+1)/SEMI
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      ECC = SQRT(ECC2)
*
*       Use semi-major axis and/or RMIN for perturber selection.
      IF (EB.LT.EBH) THEN
          RAP = SEMI*(1.0 + ECC)
      ELSE
*       Include a tapered criterion depending on energy for soft binaries.
          IF (EB.LT.0.0) THEN
              ZFAC = 1.0 + ABS(EB - EBH)/ABS(EBH)
          ELSE
              ZFAC = 1.0
          END IF
*       Adopt actual apocentre for perturber selection if R > SEMI.
          IF (SEMI.GT.0.0D0.AND.R(IPAIR).GT.SEMI) THEN
              RAP = SEMI*(1.0 + ECC)
          ELSE
              RAP = MAX(ZFAC*RMIN,R(IPAIR))
          END IF
*       Ensure extra perturbers at new regularizations (R may be small).
          IF (IPHASE.GT.0.AND.SEMI.GT.0.0) THEN
              RAP = SEMI*(1.0 + ECC)
          END IF
      END IF
*
      RCRIT2 = 2.0*RAP**2/BODY(I)
      RCRIT3 = RCRIT2*RAP/GMIN
*       Base fast search on maximum binary mass (2*BODY1).
      RCRIT2 = 2.0*RCRIT2*BODY1*CMSEP2
*
*       Select new perturbers from the neighbour list.
      NNB1 = 1
      NNB2 = LIST(1,I) + 1
      DO 10 L = 2,NNB2
          J = LIST(L,I)
          W1 = X(1,J) - X(1,I)
          W2 = X(2,J) - X(2,I)
          W3 = X(3,J) - X(3,I)
          RSEP2 = W1*W1 + W2*W2 + W3*W3
*       Include any merged c.m. or chain c.m. bodies in the fast test.
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
      IF (NNB1.EQ.1.AND.SEMI.GT.0.0) THEN
*       Retain the previous perturber after partial reflection.
          IF (TDOT2(IPAIR).GT.0.0D0.AND.KZ(25).GT.0) THEN
              NNB1 = 2
          ELSE IF (KZ(27).LE.1) THEN
*       Specify one unperturbed period at apocentre (NB! check STEP(I)).
              STEP(I1) = TWOPI*SEMI*SQRT(SEMI/BODY(I))
              STEP(I1) = MIN(STEP(I1),STEP(I))
*       Set GAMMA < GMIN since CMFIRR uses resolved KS otherwise (bug 4/00).
              GAMMA(IPAIR) = 0.5*GMIN
          ELSE
*       Maintain perturbed motion during circularization (cf. KSTIDE test).
              RP = SEMI*(1.0D0 - ECC)
              RT = 4.0*MAX(RADIUS(I1),RADIUS(I1+1))
              IF (RP.LT.0.99*RT.AND.KSTAR(I).NE.20) THEN
                  IF (LIST(1,I1).GT.0) THEN
                      NNB1 = 2
                  ELSE
*       Adopt orbital period since current value could be small many times.
                      STEP(I1) = TWOPI*SEMI*SQRT(SEMI/BODY(I))
                      STEP(I1) = MIN(STEP(I1),STEP(I))
                  END IF
              ELSE
                  STEP(I1) = TWOPI*SEMI*SQRT(SEMI/BODY(I))
                  STEP(I1) = MIN(STEP(I1),STEP(I))
                  GAMMA(IPAIR) = 0.5*GMIN
              END IF
          END IF
*       Copy all neighbours (< LMAX-1) for soft binary in pericentre region.
          IF (R(IPAIR).LT.SEMI.AND.EB.GT.EBH) THEN
              DO 15 L = 2,NNB2
                  IF (NNB1.LT.LMAX-3) THEN
                      NNB1 = NNB1 + 1
                      LIST(NNB1,I1) = LIST(L,I)
                  END IF
   15         CONTINUE
          END IF
      END IF
*
*       Save perturber membership.
      LIST(1,I1) = NNB1 - 1
*
      RETURN
*
      END
