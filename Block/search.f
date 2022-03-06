      SUBROUTINE SEARCH(I,IKS)
*
*
*       Close encounter search.
*       -----------------------

      INCLUDE 'common6.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
*
*
*       Increase counter for regularization attempts and set critical step.
      NKSTRY = NKSTRY + 1
      RJMIN2 = 1.0
      DTS = MAX(SMIN,8.0*STEP(I))
*
      FMAX = 0.0
      NCLOSE = 0
*       Find dominant neighbours by selecting all STEP(J) <= 4*STEP(I).
      L = LIST(1,I) + 2
    2 L = L - 1
      IF (L.LT.2) GO TO 30
      IF (LIST(L,I).LE.N) GO TO 4
*
*       Check first whether any c.m. with small step is dominant.
      J = LIST(L,I)
*       Include mass condition (STEP may be large).
      IF (STEP(J).GT.DTS.AND.BODY(J).LT.10.0*BODY(I)) GO TO 2
      A1 = X(1,J) - X(1,I)
      A2 = X(2,J) - X(2,I)
      A3 = X(3,J) - X(3,I)
      RIJ2 = A1*A1 + A2*A2 + A3*A3
*       Include more distant massive perturbers from the neighbour list.
      IF (RIJ2.GT.9.0*RMIN2) GO TO 2
*
      FIJ = (BODY(I) + BODY(J))/RIJ2
      IF (FMAX.LT.FIJ) FMAX = FIJ
*       Abandon further search if c.m. force exceeds half total force.
      IF (FMAX**2.LT.F(1,I)**2 + F(2,I)**2 + F(3,I)**2) THEN
          NCLOSE = NCLOSE + 1
          JLIST(NCLOSE) = J
          GO TO 2
      ELSE
          GO TO 30
      END IF
*
*       Continue searching single particles with current value of FMAX.
    4 JCOMP = 0
      RX2 = 16.0*RMIN22
      DO 6 K = L,2,-1
          J = LIST(K,I)
*         IF (STEP(J).GT.DTS) GO TO 6
          A1 = X(1,J) - X(1,I)
          A2 = X(2,J) - X(2,I)
          A3 = X(3,J) - X(3,I)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*         IF (RIJ2.LT.4.0*RMIN2) THEN
          IF (RIJ2.LT.RX2) THEN
              NCLOSE = NCLOSE + 1
              JLIST(NCLOSE) = J
*       Record index of every single body with small step inside 2*RMIN.
               FIJ = (BODY(I) + BODY(J))/RIJ2
              IF (FIJ.GT.FMAX) THEN
                  FMAX = FIJ
*       Save square distance and global index of dominant body.
                  RJMIN2 = RIJ2
                  JCOMP = J
              END IF
          END IF
    6 CONTINUE
*
*       See whether dominant component is a single particle inside RMIN.
      IF (JCOMP.LT.IFIRST.OR.JCOMP.GT.N) GO TO 30
*       Accept one single candidate inside 2*RMIN (which makes PERT = 0).
      IF (RJMIN2.GT.RMIN22) GO TO 30
*
      RDOT = (X(1,I) - X(1,JCOMP))*(XDOT(1,I) - XDOT(1,JCOMP)) +
     &       (X(2,I) - X(2,JCOMP))*(XDOT(2,I) - XDOT(2,JCOMP)) +
     &       (X(3,I) - X(3,JCOMP))*(XDOT(3,I) - XDOT(3,JCOMP))
*
      FMAX0 = FMAX
*       Ensure one perturber (& JCOMP) in FPERT estimate.
      IF (NCLOSE.EQ.-1) THEN
          NNB1 = LIST(1,I) + 1
          NCLOSE = 1
          JLIST(NCLOSE) = JCOMP 
          FMAX = 0.0
          DO 20 L = 2,NNB1
              J = LIST(L,I)
              IF (J.EQ.JCOMP) GO TO 20
              RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2
     &                                    + (X(3,I) - X(3,J))**2
              FIJ = (BODY(I) + BODY(J))/RIJ2
              IF (FIJ.GT.MAX(FMAX,FMAX0)) THEN
                  FMAX = FIJ
                  NCLOSE = NCLOSE + 1
                  JLIST(NCLOSE) = J
                  JCOMP = J
              END IF
   20     CONTINUE
      END IF
*
*       Only select approaching particles (include nearly circular case).
      RIJMIN = SQRT(RJMIN2)
      IF (RDOT.GT.0.02*SQRT((BODY(I) + BODY(JCOMP))*RIJMIN)) GO TO 30
*
*       Evaluate vectorial perturbation due to NCLOSE bodies.
      CALL FPERT(I,JCOMP,NCLOSE,JLIST,PERT)
*
*       Accept #I & JCOMP if the relative motion is dominant (GI < 0.01).
      BCM = BODY(I) + BODY(JCOMP)
      GI = PERT*RJMIN2/BCM
      IF (GI.GT.0.01) THEN
*         IF (KZ(4).GT.0.AND.TIME-TLASTT.GT.4.44*TCR/FLOAT(N))
*    &                                             CALL EVOLVE(JCOMP,0)
          GO TO 30
      END IF
*
*       Exclude any c.m. body of compact subsystem (TRIPLE, QUAD or CHAIN).
      DO 8 ISUB = 1,NSUB
          NAMEI = NAMES(1,ISUB)
          IF (NAMEI.EQ.NAME(I).OR.NAMEI.EQ.NAME(JCOMP)) GO TO 30
    8 CONTINUE
*
*       Also check possible c.m. body of chain regularization (NAME = 0).
      IF (NCH.GT.0) THEN
          IF (NAME(I).EQ.0.OR.NAME(JCOMP).EQ.0) GO TO 30
      END IF
*
*       Save index and increase indicator to denote new regularization.
      ICOMP = I
      IKS = IKS + 1
*     WRITE (6,33)  NCLOSE, NAME(ICOMP), NAME(JCOMP), SQRT(RX2), GI
*  33 FORMAT (' FINAL    NCL NM RX GI ',I4,2I6,1P,5E10.2)
*
   30 RETURN
*
      END
