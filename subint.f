      SUBROUTINE SUBINT(IQ,I10,STEP0)
*
*
*       Decision-making for subsystems.
*       -------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      INTEGER  KSLIST(KMAX)
*
*       See whether to advance KS solutions at first new block time.
    1 IF (NPAIRS.GT.0) THEN
*       Obtain list of all KS pairs due in interval DTB.
          IF (TBLIST.LE.TBLOCK.OR.NPAIRS.NE.NBPREV) THEN
              IF (DTB.EQ.0.0D0) THEN
                  DTB = MAX(DTMIN,TBLOCK - TPREV)
              END IF
    2         TBLIST = TPREV + DTB
              NNTB = 0
              DO 3 JPAIR = 1,NPAIRS
                  J1 = 2*JPAIR - 1
                  IF (T0(J1) + STEP(J1).LE.TBLIST) THEN
                      NNTB = NNTB + 1
                      KBLIST(NNTB) = J1
                  END IF
    3         CONTINUE
*       Increase interval on zero membership.
              IF (NNTB.EQ.0) THEN
                  DTB = 2.0*DTB
                  GO TO 2
              END IF
*       Stabilize interval on membership of NPAIRS/4.
              RATIO = 0.25*FLOAT(NPAIRS+3)/FLOAT(NNTB)
              DTB = SQRT(RATIO)*DTB
              DTB = MAX(2.0*(TBLOCK - TPREV),DTB)
*       Do not exceed next output boundary because of KS escaper removal.
              TBLIST = MIN(TADJ,TBLIST)
          END IF
*
*       Form list of any KS pairs due in the current block-step.
          KSB = 0
          DO 4 L = 1,NNTB
              J1 = KBLIST(L)
              IF (T0(J1) + STEP(J1).LE.TBLOCK) THEN
                  IF(LIST(1,J1).GT.0)THEN
                     IBINP = IBINP + 1
                  ELSE
                     IBINU = IBINU + 1
                  END IF
                  KSB = KSB + 1
                  KSLIST(KSB) = J1
              END IF
    4     CONTINUE
*
*       Search members of KSLIST for smallest value of T0 + STEP.
    5     TKSMIN = 1.0E+20
          DO 10 K = 1,KSB
              J1 = KSLIST(K)
              TKS = T0(J1) + STEP(J1)
              IF (TKS.LE.TKSMIN) THEN
                  I1 = J1
                  TKSMIN = TKS
              END IF
   10     CONTINUE
*
*       See whether the smallest KS time falls due before next block step.
          IF (TKSMIN.LE.TBLOCK) THEN
              TIME = TKSMIN
              CALL KSINT(I1)
*
*       Set KS indicator on termination, multiple regularization or merger.
              IF (IPHASE.NE.0) THEN
                  TBLIST = TPREV
                  IF (IQ.EQ.0.OR.IPHASE.LT.0) THEN
                      IQ = IPHASE
*       Do not overwrite new location #I1 in case of collision (IQ < 0).
                      IF (IQ.GT.0) THEN
                          I10 = I1
                          STEP0 = STEP(I1)
                          STEP(I1) = 1.0
                      END IF
                  END IF
*       Reset decision indicator and make new list on change of KS sequence.
                  IF (IPHASE.GT.0) THEN
                      IPHASE = 0
                  ELSE
                      IPHASE = 0
                      GO TO 1
                  END IF
              END IF
              GO TO 5
          END IF
*       Copy original block time at end of KS treatment.
          TIME = TBLOCK
          NBPREV = NPAIRS
      END IF
*
*       Check time for advancing any triple, quad or chain regularization.
      IF (NSUB.GT.0) THEN
   30     TSUB = 1.0D+10
          DO 40 L = 1,NSUB
              IF (TS(L).LT.TSUB) THEN
                  ISUB = L
                  TSUB = TS(L)
              END IF
   40     CONTINUE
*
          IF (TSUB.LE.TBLOCK) THEN
              TIME = TSUB
*       Decide between triple, quad or chain.
              IF (ISYS(ISUB).EQ.1) THEN
*       Update unperturbed size of subsystem and copy c.m. step.
                  CALL EXTEND(ISUB)
                  CALL TRIPLE(ISUB)
              ELSE IF (ISYS(ISUB).EQ.2) THEN
                  CALL EXTEND(ISUB)
                  CALL QUAD(ISUB)
              ELSE
                  IF (STEPS(ISUB).LT.0.0D0) THEN
                      STEPS(ISUB) = 1.0D-10
                      GO TO 50
                  END IF
                  CALL CHAIN(ISUB)
                  IF (STEPS(ISUB).LT.0.0D0) THEN
                      STEPS(ISUB) = 1.0D-10
                      GO TO 50
                  END IF
              END IF
*
*       Check for termination (set TPREV < TIME and set IQ < 0).
              IF (ISUB.LT.0.OR.IPHASE.LT.0) THEN
                  TPREV = TIME - STEP(NTOT)
                  IQ = -1
              END IF
              GO TO 30
          END IF
   50     TIME = TBLOCK
      END IF
*
      RETURN
*
      END
