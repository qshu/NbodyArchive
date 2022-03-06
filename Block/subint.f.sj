      SUBROUTINE SUBINT(IQ,I10)
*
*
*       Decision-making for subsystems.
*       -------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'kscmn6.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CONNECT/  TIME_CH
      INTEGER LISTKS(KMAX),NEXTL(KMAX)
      PARAMETER (NLMAX=3000)
      SAVE LY,NK,IT
      DATA LY,IT /0,0/
*     SAVE COMP,TT1
*     DATA COMP /0.0D0/
*
*
*     TT1 = DBLTIM()  ! Needs wtime.o from GPU2 (also in Makefile)
*
*       Initialize NEW CHAIN delay time after start/restart.
      IF (IT.EQ.0) THEN
          TIME_CH = TTOT
          READ (5,*) CLIGHT
          IT = 1
      END IF
*
*       See whether to advance any KS solutions at start of block-step.
      IF (NPAIRS.GT.0) THEN
*       Make a list of all members due up to end of block-step.
          I1 = -1
          NK = 0
          DO 2 IP = 1,NPAIRS
            I1 = I1 + 2
            IF (T0(I1) + STEP(I1).LE.TBLOCK) THEN
              NK = NK + 1
              LISTKS(NK) = I1
            END IF
    2     CONTINUE
          IF (NK.EQ.0) GO TO 30
*
*       Determine quantized value of TMIN (next small block-step).
    5     TMIN = 1.0D+10                ! Bug fix 3/19.
          DO 6 L = 1,NK
            I1 = LISTKS(L)
            TMIN = MIN(T0(I1) + STEP(I1),TMIN)
    6     CONTINUE
*
*       Update time and form list of binaries due at TMIN.
          TIME = TMIN
          LENGTH = 0
          DO 8 L = 1,NK
            I1 = LISTKS(L)
            IF (T0(I1) + STEP(I1).EQ.TMIN) THEN
              LENGTH = LENGTH + 1
              NEXTL(LENGTH) = I1
            END IF
    8     CONTINUE
*
*       Check possible chain integration on zero LENGTH.
          IF (LENGTH.EQ.0) THEN
            GO TO 30
          END IF
*
*       Clean current length of ISTAT.
          DO 9 L = 1,LENGTH
            ISTAT(L) = 0
    9     CONTINUE
*
*       Begin KS integration (sequential or parallel).
          IF (LENGTH.LE.NLMAX) THEN
            DO 10 LI = 1,LENGTH
              I1 = NEXTL(LI)
              CALL XINIT(I1)
              CALL KSINT(I1,LI)
              CALL XUPDATE(I1)
   10       CONTINUE
            NBPRED = NBPRED + LENGTH
          ELSE
*       Perform the parallel execution loop.
            DO 20 L = 1,LENGTH
              I1 = NEXTL(L)
              CALL XINIT(I1)
   20       CONTINUE
!$omp parallel do private(LI,I1) schedule(guided)
            DO 21 LI = 1,LENGTH
              I1 = NEXTL(LI)
              CALL KSINTP(I1,LI)
   21       CONTINUE
!$omp end parallel do
            DO 22 LI = 1,LENGTH
              I1 = NEXTL(LI)
              CALL XUPDATE(I1)
   22       CONTINUE
*       Note that unperturbed KS are now included in parallel step counter.
            NSTEPU = NSTEPU + LENGTH
          END IF
*
*       Search non-zero flags (ISTAT < 0 means collision).
          IF (IQ.EQ.0) THEN
            DO 25 LI = 1,LENGTH
              IF (ISTAT(LI).NE.0) THEN
                IQ = ISTAT(LI)
                I10 = NEXTL(LI)
                IPAIR = KVEC(I10)
*       Note IPHASE must be defined for chain (not needed for collision).
                IF (IQ.GE.4.AND.IQ.NE.7) THEN
!       TODO remove in final version
                  if (IQ/=IPHASEX(IPAIR)) write(6,*) 'WARN SUBINT#1',IQ,
     &                                              IPHASEX(IPAIR),IPAIR
                  if (IPAIR/=KSPAIRX(IPAIR)) write(6,*) 'WARN SUBINT#2',
     &                                    IQ,IPHASEX(IPAIR),IPAIR,KSPAIR
                  IPHASEX(IPAIR) = IQ
                  CALL DELAY(0,IPAIR)
                END IF
                GO TO 28
              END IF 
   25       CONTINUE
          END IF
*
*       Treat collision explicitly before quitting (suppressed in KSINT).
   28     IF (IQ.LT.0) THEN
            IP = KVEC(I10)
            KSPAIR = IP
            QPERI = R(IP)
            IQCOLL = -2
            CALL CMBODY(QPERI,2)
*       Ensure block-step time for normal continuation.
            TIME = TBLOCK
*       Note rest of the block will be done next time (need GO TO 999 first).
            GO TO 60
          END IF
*
*       Complete all block-steps before special procedure (except coll).
*       Check more block-step levels unless collision (final exit at TBLOCK).
          IF (TMIN.LT.TBLOCK.AND.IQ.EQ.0) GO TO 5
*       Copy original block time at end of KS treatment.
          TIME = TBLOCK
      END IF
*
*       Check time for advancing any triple, quad or chain regularization.
   30 IF (NSUB.GT.0) THEN
   40     TSUB = 1.0D+10
          DO 45 L = 1,NSUB
              IF (TS(L).LT.TSUB) THEN
                  ISUB = L
                  TSUB = TS(L)
              END IF
   45     CONTINUE
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
                  IF (ISUB.GT.0) THEN
                      IF (STEPS(ISUB).LT.0.0D0) THEN
                          STEPS(ISUB) = 1.0D-10
                          GO TO 50
                      END IF
                  END IF
              END IF
*
*       Check for termination (set TPREV < TIME and set IQ < 0).
              IF (ISUB.LT.0.OR.IPHASE.LT.0) THEN
                  TPREV = TIME - STEP(NTOT)
                  IQ = -1
              END IF
              GO TO 40
          END IF
   50     TIME = TBLOCK
      END IF
*
*     TT2 = DBLTIM()
*     COMP = COMP + (TT2 - TT1)
*     IF (ABS(TIME-TCRIT).LT.1.0D-05) WRITE (6,70)  COMP
*  70 FORMAT (' SUBINT   ',1P,E12.3)
*
   60 RETURN
*
      END
