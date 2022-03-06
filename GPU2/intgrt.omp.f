      SUBROUTINE INTGRT
*
*
*       N-body integrator flow control.
*       -------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/PREDICT/ TPRED(NMAX)
      PARAMETER (NIMAX=1024,NPMAX=16)
      REAL*8   H2I(NIMAX),XI(3,NIMAX),VI(3,NIMAX),GPUACC(3,NIMAX),
     &         GPUJRK(3,NIMAX),GPUPHI(NIMAX),GF(3,NMAX),GFD(3,NMAX)
      INTEGER  NXTLST(NMAX),LISTQ(NMAX),NL(20)
      INTEGER  IR(NMAX),LISTGP(LMAX,NIMAX)
      LOGICAL LOOP,LSTEPM
      SAVE IQ,ICALL,NQ,LQ,LOOP,LSTEPM,STEPM,ISAVE,JSAVE,ISTART,NNPREV
      DATA IQ,ICALL,LQ,LOOP,LSTEPM,STEPM /0,2,11,.TRUE.,.FALSE.,0.03125/
      DATA ISAVE,JSAVE,ISTART /0,0,0/
*     SAVE CPRED,CPRED2,CPRED3
*     DATA CPRED,CPRED2,CPRED3 /0.0D0,0.0D0,0.0D0/
*
*
*       Enforce level search on return, except new and terminated KS.
      IF (IPHASE.NE.1.AND.IPHASE.NE.2) LOOP = .TRUE.
*
*       Update quantized value of STEPM for large N (first time only).
      IF (.NOT.LSTEPM.AND.NZERO.GT.1024) THEN
          K = (FLOAT(NZERO)/1024.0)**0.333333
          STEPM = 0.03125D0/2**(K-1)
          LSTEPM = .TRUE.
      END IF
*
*       Open the GPU libraries on each new run (note nnbmax = NN is printed).
      IF (ISTART.EQ.0) THEN
          NN = N
          NNPREV = NN
          CALL GPUNB_OPEN(NN)
*       Set larger value for GPUIRR (but note further increase of NTOT).
          NNN = NTOT + 10
          CALL GPUIRR_OPEN(NNN,LMAX)
          ISTART = 1
      END IF
*
*       Search for high velocities after escape or KS/chain termination.
  999 IF (KZ(37).GT.0.AND.(IPHASE.EQ.-1.OR.IPHASE.GE.2)) THEN
          CALL HIVEL(0)
      END IF
*
*       Update equence on the GPU after major changes.
!$omp parallel do private(I)
      DO 995 I = IFIRST,NTOT
          CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),F(1,I),FDOT(1,I),
     &                                            BODY(I),T0(I))
          CALL GPUIRR_SET_LIST(I,LIST(1,I))
  995 CONTINUE
!$omp end parallel do
*
*       Reset control & regularization indicators.
      IPHASE = 0
      IKS = 0
      DTM = 1.0
      TPREV = TIME
*       Initialize end-point of integration times and set DTM & TPRED = TIME.
      DO 1000 I = IFIRST,NTOT
          TNEW(I) = T0(I) + STEP(I)
          DTM = MIN(DTM,STEP(I))
          TPRED(I) = -1.0
          CALL JPRED(I)
 1000 CONTINUE
*
*       Determine level for the smallest step (ignore extreme values).
      LQS = 20
      DO 1001 L = 6,20
          IF (DTM.EQ.DTK(L)) THEN
              LQS = L
          END IF
 1001 CONTINUE
*
*       Specify upper level for optimized membership.
      LQB = LQS - 4
      IF (IQ.LT.0) ICALL = 0
      IQ = 0
*       Enforce new block step search initially and on significant change.
      TLISTQ = TIME
*
*       Check updating new list of block steps with T0 + STEP =< TLISTQ.
    1 ICALL = ICALL + 1
*     TT1 = DBLTIM()
*       Reset TMIN second & third time after change to catch new chain step.
      IF (TIME.GE.TLISTQ.OR.ICALL.LE.3) THEN
*       Update interval by optimization at major times (sqrt of N-NPAIRS).
          IF (DMOD(TLISTQ,2.0D0).EQ.0.0D0.OR.LOOP) THEN
              LOOP = .FALSE.
              DO 10 L = 1,20
                  NL(L) = 0
   10         CONTINUE
!$omp parallel do private(I)
              DO 14 I = IFIRST,NTOT
*       Count steps at five different levels for the smallest values.
                  DO 12 L = LQB,LQS
                      IF (STEP(I).LT.DTK(L)) NL(L) = NL(L) + 1
   12             CONTINUE
   14         CONTINUE
!$omp end parallel do
              NLSUM = 0
*       Determine interval by summing smallest steps until near sqrt(N-N_b).
              NSQ = SQRT(FLOAT(N - NPAIRS))
              LQ = LQS
              DO 15 L = LQS,LQB,-1
                  NLSUM = NLSUM + NL(L)
                  IF (NLSUM.LE.NSQ) LQ = L
   15         CONTINUE
*             WRITE (6,16)  TIME+TOFF,NQ,NLSUM,LQ,(NL(K),K=LQB,LQS)
*  16         FORMAT (' LEVEL CHECK:    T NQ NLSUM LQ NL  ',
*    &                                  F9.3,3I5,2X,7I4)
          END IF
*
*       Increase interval by optimized value.
          NQ = 0
          TMIN = 1.0D+10
   18     TLISTQ = TLISTQ + DTK(LQ)
          DO 20 I = IFIRST,NTOT
              IF (TNEW(I).LE.TLISTQ) THEN
                  NQ = NQ + 1
                  LISTQ(NQ) = I
                  TMIN = MIN(TNEW(I),TMIN)
              END IF
   20     CONTINUE
*       Increase interval in rare case of zero membership.
          IF (NQ.EQ.0) GO TO 18
*       Make a slight adjustment for high levels and small membership.
          IF (LQ.LE.15.AND.NQ.LE.2) LQ = LQ - 1
      END IF
*
*       Find all particles in next block (TNEW = TMIN).
      CALL INEXT(NQ,LISTQ,TMIN,NXTLEN,NXTLST)
*
*       Set new time and save block time (for regularization terminations).
      I = NXTLST(1)
      TIME = T0(I) + STEP(I)
      TBLOCK = TIME
*
*     WRITE (6,22)  I, NXTLEN, NSTEPU, NSTEPI, TIME, STEP(I), STEPR(I)
*  22 FORMAT (' INTGRT   I LEN #U #I T S SR  ',2I6,2I11,F9.4,1P,2E10.2)
*     IF (STEP(I).LT.1.0D-08) STOP
*     CALL FLUSH(6)
*
*     TT2 = DBLTIM()
*     CPRED3 = CPRED3 + (TT2 - TT1)
*       Re-determine list if current time exceeds boundary.
      IF (TIME.GT.TLISTQ) GO TO 1
*
*       Check option for advancing interstellar clouds.
      IF (KZ(13).GT.0) THEN
          CALL CLINT
      END IF
*
*       Check optional integration of cluster guiding centre.
      IF (KZ(14).EQ.3.OR.KZ(14).EQ.4) THEN
          IF (KZ(14).EQ.3.AND.DMOD(TIME,STEPX).EQ.0.0D0) THEN
              CALL GCINT
          END IF
*       Include mass loss by gas expulsion (Kroupa et al. MN 321, 699).
          IF (MPDOT.GT.0.0D0.AND.TIME + TOFF.GT.TDELAY) THEN
              MP = MP0/(1.0 + MPDOT*(TIME + TOFF - TDELAY))
          END IF
      END IF
*
*       Include commensurability test (may be suppressed if no problems).
*     IF (DMOD(TIME,STEP(I)).NE.0.0D0) THEN
*         WRITE (6,25)  I, NAME(I), NSTEPI, TIME, STEP(I), TIME/STEP(I)
*  25     FORMAT (' DANGER!   I NM # TIME STEP T/DT ',
*    &                        2I6,I11,F12.5,1P,E9.1,0P,F16.4)
*         STOP
*     END IF
*
*       Check for new regularization at end of previous block.
      IF (IKS.GT.0) THEN
          IPHASE = 1
*       Copy the saved KS component indices and time.
          ICOMP = ISAVE
          JCOMP = JSAVE
          TIME = TSAVE
          GO TO 100
      END IF
*
*     TT1 = DBLTIM()
*       Check next adjust time before beginning a new block.
      IF (TIME.GT.TADJ) THEN
          TIME = TADJ
          IPHASE = 3
          GO TO 100
      END IF
*
*       Check output time in case DTADJ & DELTAT not commensurate.
      IF (TIME.GT.TNEXT) THEN
          TIME = TNEXT
          CALL OUTPUT
          GO TO 1
      END IF
*
*       See whether to advance ARchain or KS at first new time.
      IF (TIME.GT.TPREV) THEN
          CALL SUBINT(IQ,I10)
          IF (IQ.LT.0) GO TO 999
      END IF
*
*       Define array for regular force condition (STEPR is commensurate).
      NFR = 0
      DO 28 L = 1,NXTLEN
          J = NXTLST(L)
          IF (TNEW(J).GE.T0R(J) + STEPR(J)) THEN
              NFR = NFR + 1
              IR(NFR) = J
          END IF
   28 CONTINUE
*
*       Decide between merging of neighbour lists or full N prediction.
      IF (NXTLEN.LE.75.AND.NFR.EQ.0) THEN
*
*       Predict active particles in the library.
          CALL GPUIRR_PRED_ACT(NXTLEN,NXTLST,TIME)
      ELSE
          CALL GPUIRR_PRED_ALL(IFIRST,NTOT,TIME)
          NNPRED = NNPRED + 1
*!$omp parallel do private(S, S1, S2)
*         DO 40 J = IFIRST,NTOT
*             IF (BODY(J).EQ.0.0D0) GO TO 40
*             S = TIME - T0(J)
*             S1 = 1.5*S
*             S2 = 2.0*S
*             X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S +X0(1,J)
*             X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S +X0(2,J)
*             X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S +X0(3,J)
*             XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
*             XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
*             XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
*  40     CONTINUE
*!$omp end parallel do
*         CALL CXVPRED(IFIRST,NTOT,TIME,T0,X0,X0DOT,F,FDOT,X,XDOT,TPRED)
      END IF
*
*       Resolve all active KS binaries.
      DO 45 L = 1,NXTLEN
          J = NXTLST(L)
          IF (J.GT.N) THEN
              JPAIR = J - N
              IF (LIST(1,2*JPAIR-1).GT.0) THEN
                  CALL JPRED(J)
              END IF
          END IF
   45 CONTINUE
*
*       Save new time (output time at TIME > TADJ) and increase # blocks.
      TPREV = TIME
      NBLOCK = NBLOCK + 1
      TMIN = 1.0D+10
      IKS0 = IKS
*
*       Predict chain variables and perturber cordinates at new block-time.
      IF (NCH.GT.0) THEN
          IF (TNEW(ICH).GT.TBLOCK) THEN
              CALL XVPRED(ICH,0)
              TPRED(ICH) = TIME
          END IF
          CALL XCPRED(2)
      END IF
*
*       Evaluate all irregular forces in the library.
*     DO 46 II = 1,NXTLEN
*         I = NXTLST(II)
*         CALL GPUIRR_FIRR(I,GF(1,II),GFD(1,II))
*  46 CONTINUE
      CALL GPUIRR_FIRR_VEC(NXTLEN,NXTLST,GF,GFD)
*
      IF (NXTLEN.LE.NPMAX) THEN
*
*       Obtain all irregular forces in one loop.
      DO 50 II = 1,NXTLEN
*       Define pointer array and scalar for regular force condition.
          I = NXTLST(II)
          IF (TNEW(I).GE.T0R(I) + STEPR(I)) THEN
              IR1 = 1
          ELSE
              IR1 = 0
          END IF
*
*       Advance the irregular step and copy corrector to the GPU.
          CALL NBINT(I,IKS,IR1,GF(1,II),GFD(1,II))
          CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),F(1,I),FDOT(1,I),
     &                                            BODY(I),T0(I))
*
*       Save indices and TIME of first KS candidates in the block.
          IF (IKS0.EQ.0.AND.IKS.GT.0) THEN
              ISAVE = ICOMP
              JSAVE = JCOMP
              TSAVE = TIME
          END IF
   50 CONTINUE
*
      ELSE
*       Obtain all irregular forces in parallel.
!$omp parallel do private(II, I, IR1)
      DO II = 1,NXTLEN
*       Define pointer array and scalar for regular force condition.
          I = NXTLST(II)
          IF (TNEW(I).GE.T0R(I) + STEPR(I)) THEN
              IR1 = 1
          ELSE
              IR1 = 0
          END IF
*
*       Advance the irregular steps (parallel version).
*         CALL GPUIRR_FIRR(I,GF(1,II),GFD(1,II))
          CALL NBINTP(I,IR1,GF(1,II),GFD(1,II))
          CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),F(1,I),FDOT(1,I),
     &                                            BODY(I),T0(I))
      END DO
!$omp end parallel do
          NSTEPI = NSTEPI + NXTLEN
      END IF
*
*       Predict all particles and send to library (general case NFR > 0).
      IF (NFR.GT.0) THEN
      NN = NTOT - IFIRST + 1
*     TT3 = DBLTIM()
      CALL CXVPRED(IFIRST,NTOT,TIME,T0,X0,X0DOT,F,FDOT,X,XDOT,TPRED)
*     TT4 = DBLTIM()
*     CPRED = CPRED + (TT4 - TT3)
      CALL GPUNB_SEND(NN,BODY(IFIRST),X(1,IFIRST),XDOT(1,IFIRST))
*     CALL GPUNB_REGF(NI,H2I,XI,VI,GPUACC,GPUJRK,LMAX,NNBMAX,LISTGP)
*
*       Perform regular force loop.
      NOFL(1) = 0
      JNEXT = 0
      DO 55 II = 1,NFR,NIMAX
  550     NI = MIN(NFR-JNEXT,NIMAX)
*       Copy neighbour radius and state vector for each block.
!$omp parallel do private(I, K)
          DO 52 LL = 1,NI
              I = IR(JNEXT+LL)
              H2I(LL) = RS(I)**2
              DO 51 K = 1,3
                  XI(K,LL) = X(K,I)
                  VI(K,LL) = XDOT(K,I)
   51         CONTINUE
   52     CONTINUE
!$omp end parallel do
*
          CALL GPUNB_REGF(NI,H2I,XI,VI,GPUACC,GPUJRK,GPUPHI,LMAX,
     &                                               NNBMAX,LISTGP)
*       Copy neighbour list after overflow check.
          DO 54 LL = 1, NI
              I = IR(JNEXT + LL)
              NNB = LISTGP(1,LL)
*       Repeat the last block with reduced RS(I) on NNB < 0.
              IF (NNB.LT.0) THEN
                  WRITE (41,556)  NSTEPR, NAME(I), LIST(1,I), RS(I)
  556             FORMAT (' OVERFLOW!   #R NAME NB0 RS ',I11,I6,I4,F8.2)
                  CALL FLUSH(41)
                  NOFL(1) = NOFL(1) + 1
                  RS(I) = 0.9*RS(I)
                  GO TO 550
              END IF
   54     CONTINUE
!$omp parallel do private(I, ITEMP, NNB, L1, L)
          DO 56 LL=1, NI
              I = IR(JNEXT + LL)
              NNB = LISTGP(1,LL)
              L1 = 1
              DO 53 L = 2,NNB+1
*       Note GPU address starts from 0 (hence add IFIRST to neighbour list).
                  ITEMP = LISTGP(L,LL) + IFIRST
                  IF (ITEMP.NE.I) THEN
                      L1 = L1 + 1
                      LISTGP(L1,LL) = ITEMP
                  END IF
   53         CONTINUE
              LISTGP(1,LL) = L1 - 1
              CALL GPUIRR_SET_LIST(I,LISTGP(1,LL))
   56     CONTINUE
!$omp end parallel do
*
          CALL GPUIRR_FIRR_VEC(NI,IR(II),GF(1,1),GFD(1,1))
!$omp parallel do private(I)
          DO 57 LL = 1, NI
              I = IR(JNEXT+LL)
*       Obtain new irregular force and perform correction.
*             CALL GPUIRR_FIRR(I,GF(1,LL),GFD(1,LL))
              CALL GPUCOR(I,XI(1,LL),VI(1,LL),GPUACC(1,LL),GPUJRK(1,LL),
     &                               GF(1,LL),GFD(1,LL),LISTGP(1,LL))
*       Update neighbour lists on GPU (mostly changed).
              CALL GPUIRR_SET_LIST(I,LIST(1,I))
*             POT = POT + BODY(I)*GPUPHI(LL)
   57     CONTINUE
!$omp end parallel do
          JNEXT = JNEXT + NI
          NSTEPR = NSTEPR + NI
   55 CONTINUE
*
*       Accumulate the sum of overflows (NOFL(1) holds current number).
      NOFL(2) = NOFL(2) + NOFL(1)
*
      END IF
*
*       Determine next block time (note STEP may shrink in GPUCOR).
      DO 350 L = 1,NXTLEN
*       Determine next block time (note STEP may shrink in GPUCOR).
          I = NXTLST(L)
          TMIN = MIN(TNEW(I),TMIN)
  350 CONTINUE
*
*       Copy current coordinates & velocities from corrected values.
!$omp parallel do private(I, L, K)
      DO 360 L = 1,NXTLEN
          I = NXTLST(L)
*         TMIN = MIN(TNEW(I),TMIN)
          DO 58 K = 1,3
              X(K,I) = X0(K,I)
              XDOT(K,I) = X0DOT(K,I)
   58     CONTINUE
  360 CONTINUE
!$omp end parallel do
*
*        Send active particle data to GPU.
!$omp parallel do private(I, L, K)
      DO 60 L = 1,NXTLEN
          I = NXTLST(L)
          CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),F(1,I),FDOT(1,I),
     &                                            BODY(I),T0(I))
   60 CONTINUE
!$omp end parallel do
*
*     TT2 = DBLTIM()
*     CPRED2 = CPRED2 + (TT2 - TT1)
*     IF (TIME.EQ.TCRIT) WRITE (6,340)  CPRED2, TT2, CPRED, CPRED3
* 340 FORMAT (' TIMING   CPRED2 TT2 CPRED CPRED3 ',1P,4E12.2)
*
*       Check integration of tidal tail members.
      IF (NTAIL.GT.0) THEN
*       Allow large quantized interval with internal iteration.
          IF (DMOD(TIME,0.25D0).EQ.0.0D0) THEN
              DO 65 J = ITAIL0,NTTOT
                  IF (TNEW(J).LE.TIME) THEN
                      CALL NTINT(J)
                  END IF
   65         CONTINUE
          END IF
      END IF
*
*       Exit on KS termination, new multiple regularization or merger.
      IF (IQ.NE.0) THEN
          NBPREV = 0
          IF (IQ.GE.4.AND.IQ.NE.7) THEN
              CALL DELAY(IQ,-1)
          ELSE
*       Ensure correct KS index (KSPAIR may denote second termination).
              KSPAIR = KVEC(I10)
              IPHASE = IQ
          END IF
          GO TO 100
      END IF
*
*       Perform optional check on high-velocity particles at major times.
      IF (KZ(37).GT.0.AND.LISTV(1).GT.0) THEN
          IF (DMOD(TIME,STEPM).EQ.0.0D0) THEN
              CALL SHRINK(TMIN)
              IF (LISTV(1).GT.0) THEN
                  CALL HIVEL(-1)
              END IF
          END IF
      END IF
*
*       Check optional mass loss time at end of block-step.
      IF (KZ(19).GT.0) THEN
*       Delay until time commensurate with 100-year step (new polynomials).
          IF (TIME.GT.TMDOT.AND.DMOD(TIME,STEPX).EQ.0.0D0) THEN
              IF (KZ(19).GE.3) THEN
                  CALL MDOT
              ELSE
                  CALL MLOSS
              END IF
              IF (IPHASE.LT.0) GO TO 999
          END IF
      END IF
*
*       Advance counters and check timer & optional COMMON save (NSUB = 0).
      NTIMER = NTIMER + 1
      IF (NTIMER.LT.NMAX) GO TO 1

      NTIMER = 0
      NSTEPS = NSTEPS + NMAX
*
      IF (NSTEPS.GE.100*NMAX.AND.NSUB.EQ.0) THEN
          NSTEPS = 0
          IF (KZ(1).EQ.3) CALL MYDUMP(1,1)
      END IF
*
*       Check option for general binary search.
      IF (KZ(4).GT.0.AND.TIME - TLASTS.GT.DELTAS) THEN  
          CALL EVOLVE(0,0)
      END IF
*
*       Include facility for termination of run (create dummy file STOP).
      OPEN (99,FILE='STOP',STATUS='OLD',FORM='FORMATTED',IOSTAT=IO)
      IF (IO.EQ.0) THEN
          CLOSE (99)
          IF (NSUB.EQ.0)  WRITE (6,70)
   70     FORMAT  (/,9X,'TERMINATION BY MANUAL INTERVENTION')
          CPU = 0.0
      END IF
*
*       Repeat cycle until elapsed computing time exceeds the limit.
      CALL CPUTIM(TCOMP)
      IF (TCOMP.LT.CPU) GO TO 1
*
*       Do not terminate during triple, quad or chain regularization.
      IF (NSUB.GT.0) THEN
*       Specify zero step to enforce termination.
          DO 75 L = 1,NSUB
              STEPS(L) = 0.0D0
   75     CONTINUE
          NTIMER = NMAX
          GO TO 1
      END IF
*
*       Obtain elapsed wall-clock time (hours, minutes & seconds).
      CALL WTIME(IHOUR,IMIN,ISEC)
      SECS = 3600.0*IHOUR + 60.0*IMIN + ISEC
      WTOT = WTOT + SECS - WTOT0
      WTOT0 = SECS
*
*       Terminate run with optional COMMON save.
      IF (KZ(1).GT.0) THEN
          CPUTOT = CPUTOT + TCOMP - CPU0
          WT = WTOT/3600.0
          CALL MYDUMP(1,1)
          WRITE (6,80)  TIME+TOFF, TCOMP, CPUTOT/60.0, ERRTOT, DETOT, WT
   80     FORMAT (/,9X,'COMMON SAVED AT TIME =',F8.2,'  TCOMP =',F7.1,
     &                 '  CPUTOT =',F6.1,'  ERRTOT =',F10.6,
     &                 '  DETOT =',F10.6,'  WTOT =',F7.1)
      END IF
*
      CALL GPUNB_CLOSE
      CALL GPUIRR_CLOSE
      STOP
*
*       Set current global time.
  100 TTOT = TIME + TOFF
*
*       Close the GPU library at termination (use modified NCRIT).
      IF (TTOT.GE.TCRIT.OR.N.LE.NCRIT+10) THEN
          CALL GPUNB_CLOSE
          CALL GPUIRR_CLOSE
          IF (N.LE.NCRIT+10) NCRIT = N
      END IF
      NNPREV = NN
*
      RETURN
*
      END
