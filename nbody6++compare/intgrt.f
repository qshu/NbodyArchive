# 1 "intgrt.F"
# 1 "<built-in>"
# 1 "<command line>"
# 1 "intgrt.F"
      SUBROUTINE INTGRT
*
*
* N-body integrator flow control.
* -------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CLUMP/ BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     & NAMES(NCMAX,5),ISYS(5)
      INTEGER JHIST,JHISTR
      LOGICAL LSHRNK,LSTEPM
      EXTERNAL SHORT
      COMMON/BLKLVL/JHIST(0:NMAX),JHISTR(0:NMAX)
      INTEGER NXTLST(NMAX),IREG(NMAX),NBLIST(NMAX),IBL(LMAX)
      COMMON/STSTAT/ TINIT,NIR,NIB,NRGL,NKS
*




      INTEGER IMPI(LMAX,NMAX),JMPI(11,maxpe)
*
      SAVE IQ,ICALL,LSTEPM,STEPM
      DATA IQ,ICALL,LSTEPM,STEPM /0,2,.FALSE.,0.03125/
*
* Update quantized value of STEPM for large N (first time only).
      IF (.NOT.LSTEPM.AND.NZERO.GT.1024) THEN
          K = (FLOAT(NZERO)/1024.0)**0.333333
          STEPM = 0.03125D0/2**(K-1)
          LSTEPM = .TRUE.
      END IF
*
* Search for high velocities after escape or KS/chain termination.
  999 IF (KZ(18).GT.0.AND.(IPHASE.EQ.-1.OR.IPHASE.GE.2)) THEN
          CALL HIVEL(0)
      END IF
*
* Reset control & regularization indicators.
      IPHASE = 0
      IKS = 0
* Initialize end-point of integration times and set TMIN.
      TMIN = 1.0D+10
* Initialize end-point of integration times and set TMIN and DTM.
      DO 1000 I = IFIRST,NTOT
         TIMENW(I) = T0(I) + STEP(I)
          IF(TIMENW(I).LT.TMIN)THEN
             TMIN = TIMENW(I)
             IMIN = I
          END IF
 1000 CONTINUE
*
      IF (IQ.LT.0) ICALL = 0
      IQ = 0
*
* Find all particles due at next block time.
    1 CONTINUE
*
      NXTLEN = 0
*
* Reset TSMALL second time after main change to catch new small steps.
      ICALL = ICALL + 1
      IF (ICALL.EQ.2) GO TO 999
*
* determine next block particles without assuming
* sorted time step list (R.Sp.)
*
      DO 5 J = IFIRST, NTOT
         IF(DABS(TIMENW(J)-TMIN).LT.DTK(40)) THEN
            NXTLEN = NXTLEN + 1
            NXTLST(NXTLEN) = J
         END IF
  5 CONTINUE
*
* if(ixxxx.gt.0)then
* if(time.gt.2.0D0)
* *print*,' nxtlen,block=',nxtlen,(name(nxtlst(k)),k=1,nxtlen)
* if(ixxxx.gt.3)ixxxx=0
* end if
*
* Update short timestep list for regularization search.
      CALL SHORT(NXTLEN,NXTLST)
*
* Set new time and save block time (for regularization terminations).
      TIME = TMIN
      TBLOCK = TIME
*
* Check option for advancing interstellar clouds.
      IF (KZ(13).GT.0) THEN
          CALL CLINT
      END IF
*
* Include commensurability test (may be suppressed if no problems).
* IF (STEP(IMIN).LT.1.0E-15.OR.DMOD(TIME,STEP(IMIN)).NE.0.0D0) THEN
* WRITE (6,1005) IMIN, NAME(IMIN), NSTEPI, TIME, STEP(IMIN),
* & TIME/STEP(IMIN)
*1005 FORMAT (' DANGER!   I NM # TIME STEP T/DT ',
* & 2I5,I11,F12.5,1P,E9.1,0P,F16.4)
* CALL ABORT
* END IF
*
* Check for new regularization at end of block.
      IF (IKS.GT.0) THEN
          TIME = TPREV
          IPHASE = 1
          GO TO 100
      END IF
*
* Check next adjust time before beginning a new block.
      IF (TIME.GT.TADJ) THEN
          TIME = TPREV
          IPHASE = 3
          GO TO 100
      END IF
*
* Also check output time in case DTADJ & DELTAT not commensurate.
      IF (TIME.GT.TNEXT) THEN
          TIME = TPREV
          CALL OUTPUT
          GO TO 1
      END IF
*
* See whether to advance any close encounters at first new time.
      IF (TIME.GT.TPREV) THEN
            call cputim(tt5)
          CALL SUBINT(IQ,I10)
            call cputim(tt6)
            ttks = ttks + (tt6-tt5)*60.
*
          IF (IQ.LT.0) GO TO 999
      END IF
*
* Check regularization criterion for single particles.
      IKS = 0
      ISMIN = 0
      DSMIN = DTMIN
* Search only in prepared list of short-step particles. (R.Sp.)
      ISHORT = LSHORT(1)
      DO 50 L = 2,ISHORT+1
          I = LSHORT(L)
* Search for minimum timestep candidate for not ordered steps (R.Sp.)
* Beware that members of LSHORT may be members of KS pair (R.Sp.)
          IF (STEP(I).LT.DTMIN.AND.STEP(I).LT.DSMIN.AND.I.LE.N.AND.
     & I.GE.IFIRST) THEN
              DSMIN = STEP(I)
              ISMIN = I
          END IF
   50 CONTINUE
*
* See whether dominant body can be regularized.
      IF(ISMIN.GT.0) THEN
          CALL SEARCH(ISMIN,IKS)
*
* Include close encounter search for low-eccentric massive binaries.
      IF (IKS.EQ.0.AND.STEP(ISMIN).LT.4.0*DTMIN) THEN
* Consider massive single bodies in absence of subsystems.
          IF (ISMIN.LE.N.AND.BODY(I).GT.2.0*BODYM.AND.NSUB.EQ.0) THEN
*
* Obtain two-body elements and relative perturbation.
              JMIN = 0
              CALL ORBIT(ISMIN,JMIN,SEMI,ECC,GI)
*
              EB = -0.5*BODY(ISMIN)*BODY(JMIN)/SEMI
              IF (EB.LT.EBH.AND.GI.LT.0.25.AND.JMIN.GE.IFIRST) THEN
                  APO = SEMI*(1.0 + ECC)
* Check eccentricity (cf. max perturbation) and neighbour radius.
                  IF (ECC.LT.0.5.AND.APO.LT.0.02*RS(ISMIN)) THEN
* PRINT*, ' KS TRY: NAM E A EB ',
* * NAME(ISMIN), NAME(JMIN), ECC, SEMI, EB
* CALL FLUSH(6)
                      IKS = IKS + 1
                      ICOMP = ISMIN
                      JCOMP = JMIN
                  END IF
              END IF
          END IF
      END IF
      END IF
*
* Check regular force condition for small block memberships.
      IR = 0
      IF (NXTLEN.LE.10) THEN
          DO 28 L = 1,NXTLEN
              J = NXTLST(L)
              IF (TIMENW(J).GE.T0R(J) + STEPR(J)) THEN
                  IR = IR + 1
              END IF
   28 CONTINUE
      END IF
*
* Choose between predicting all neighbours or full N.
* Warning do not distribute prediction on PEs for consistency
          call cputim(tt1)
      IF (NXTLEN.LT.10.AND.IR.EQ.0) THEN
*
* Initialize pointers for neighbour lists.
          DO 30 L = 1,NXTLEN
              IBL(L) = NXTLST(L)
   30 CONTINUE
*
* Merge all neighbour lists (with absent members of NXTLST added).
          CALL NBSORT(NXTLEN,IBL,NNB,NBLIST)
*
* Predict coordinates & velocities of neighbours and #I to order FDOT.
          NBPRED = NBPRED + NNB
          NBFLAG = 1
          IPRED = 0
          DO 35 L = 1,NNB
              J = NBLIST(L)
              S = TIME - T0(J)
              S1 = 1.5*S
              S2 = 2.0*S
              X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S +X0(1,J)
              X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S +X0(2,J)
              X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S +X0(3,J)
              XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
              XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
              XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
   35 CONTINUE
      ELSE
          NNPRED = NNPRED + 1
          NBFLAG = 1
          IPRED = 1
          DO 40 J = IFIRST,NTOT
              S = TIME - T0(J)
              S1 = 1.5*S
              S2 = 2.0*S
              X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S +X0(1,J)
              X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S +X0(2,J)
              X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S +X0(3,J)
              XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
              XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
              XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
 40 CONTINUE
      END IF
*
* Resolve any KS coordinates & velocities using most recent c.m.
      IF (NPAIRS.GT.0) THEN
* Resolve perturbed KS pairs with c.m. prediction after NBSORT.
      JJ = -1
      DO 45 JPAIR = 1,NPAIRS
      JJ = JJ + 2
      IF (LIST(1,JJ).GT.0) THEN
* Ignore c.m. prediction after full N loop (all active KS needed).
          IF (IPRED.EQ.0) THEN
              J = N + JPAIR
              S = TIME - T0(J)
              S1 = 1.5*S
              S2 = 2.0*S
              X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S +X0(1,J)
              X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S +X0(2,J)
              X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S +X0(3,J)
              XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
              XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
              XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
          END IF
          IZZ = -1
* Predict ALL binaries even unperturbed ones for parallel code (R.Sp.)
          ZZ = 0.0
          CALL KSRES2(JPAIR,J1,J2,ZZ,IZZ)
      END IF
   45 CONTINUE
*
      END IF
* Resolve Chain if it is in block or in neighbour lists (RS Nov. 03).
      IF(NCH.GT.0) THEN
      ICHPR = 0
* First check whether chain c.m. is in current block (ICHPR = 2).
          DO 46 L = 1,NXTLEN
              I = NXTLST(L)
              IF(I.EQ.ICH) ICHPR = 2
 46 CONTINUE
* Second check whether neighbour lists contain chain c.m. (ICHPR = 1).
          IF (ICHPR.EQ.0) THEN
              DO 47 L = 1,NXTLEN
                  I = NXTLST(L)
                  NNB1 = LIST(1,I) + 1
* Second check whether neighbour lists contain chain c.m. (ICHPR = 1).
                  DO 48 K = 2,NNB1
                      J = LIST(K,I)
                      IF (J.GT.ICH) GO TO 47
                      IF (J.EQ.ICH) ICHPR = 1
 48 CONTINUE
 47 CONTINUE
          END IF
          IF (ICHPR.GT.1) CALL CHLIST(ICH)
          IF (ICHPR.GT.0) CALL XCPRED(0)
      END IF
*
* Save new time (output time at TIME> TADJ) and increase # of blocks.
      TPREV = TIME
      NBLOCK = NBLOCK + 1
      TMIN = 1.0D+10
*
                  call cputim(tt2)
          ttnbp = ttnbp + (tt2-tt1)*60.
*
* Initialize counters for irregular & regular integrations.
      NREG = 0
*
* Advance the irregular step for all particles in the current block.
* Block-Step Level Diagnostics (R.Sp. 29.Apr. 1993)
          IF(KZ(33).GT.0)JHIST(NXTLEN) = JHIST(NXTLEN) + 1
*



*
      DO 701 L = 1,NXTLEN
*
          I = NXTLST(L)
*
      CALL NBINT(I,NBFLAG)
*
 701 CONTINUE
*
          call cputim(tt3)
          ttirr = ttirr + (tt3-tt2)*60.
*
# 443 "intgrt.F"
*
       DO 71 L = 1,NXTLEN
          I = NXTLST(L)
* Save new block step and update T0 & next time
          T0(I) = TIME
          TIMENW(I) = T0(I) + STEP(I)
*
* Set non-zero indicator for new regular force.
          IF (T0R(I) + STEPR(I).LE.TIME) THEN
              NREG = NREG + 1
              IREG(NREG) = I
          ELSE
* Extrapolate regular force & first derivatives to obtain F & FDOT.
              DTR = TIME - T0R(I)
              DO 65 K = 1,3
                  F(K,I) = 0.5*(FRDOT(K,I)*DTR + FR(K,I) + FI(K,I))
                  FDOT(K,I) = ONE6*(FRDOT(K,I) + FIDOT(K,I))
* Higher order extrapolation?
* F(K,I) = FI(K,I) + FR(K,I) + DTR*(FRDOT(K,I)
* * + DTR*(D2R(K,I)/2.D0 + DTR*D3R(K,I)/6.D0))
* FDOT(K,I) = FIDOT(K,I) + FRDOT(K,I)
* * + DTR*(D2R(K,I) + DTR*D3R(K,I)/2.D0)
* F(K,I) = F(K,I)/2.D0
* FDOT(K,I) = FDOT(K,I)/6.D0
   65 CONTINUE
          END IF
*
              DO 67 K = 1,3
                  X0(K,I) = XN(K,I)
                  X0DOT(K,I) = XNDOT(K,I)
                  D0(K,I) = FI(K,I)
                  D1(K,I) = FIDOT(K,I)
   67 CONTINUE
*
 71 CONTINUE
*
* CALL nemo_savestate(n,3,time,body,x,xdot)
*
* See whether any KS candidates are in the same block.
      IF (IKS.GT.0) THEN
* Accept same time, otherwise reduce STEP(ICOMP) and/or delay.
          IF (T0(JCOMP).EQ.T0(ICOMP)) THEN
              I = ICOMP
              ICOMP = MIN(ICOMP,JCOMP)
              JCOMP = MAX(I,JCOMP)
          ELSE IF (T0(JCOMP) + STEP(JCOMP).LT.T0(ICOMP)) THEN
              STEP(ICOMP) = 0.5D0*STEP(ICOMP)
              TIMENW(ICOMP) = T0(ICOMP) + STEP(ICOMP)
              IKS = 0
          ELSE
              IKS = 0
          END IF
      END IF
*
      NSTEPI = NSTEPI + NXTLEN
*
* Obtain total force for all particles due in the current block.
*
      IF(NREG.GT.0)THEN
*
           call cputim(tt7)
              DO 811 J=IFIRST,NTOT
              S = TIME - T0(J)
              S1 = 1.5*S
              S2 = 2.0*S
              X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S +X0(1,J)
              X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S +X0(2,J)
              X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S +X0(3,J)
              XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
              XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
              XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
 811 CONTINUE
* Resolve any KS coordinates & velocities using most recent c.m.
      IF (NPAIRS.GT.0) THEN
          JJ = -1
          DO 86 JPAIR = 1,NPAIRS
          JJ = JJ + 2
          IF (LIST(1,JJ).GT.0) THEN
              ZZ = 1.0
              IZZ = -2
* Distinguish between low and high-order prediction of U & UDOT.
              IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
              CALL KSRES2(JPAIR,J1,J2,ZZ,IZZ)
          END IF
   86 CONTINUE
      END IF
*
          call cputim(tt8)
          ttpre = ttpre + (tt8-tt7)*60.
*
* Block-Step Level Diagnostics (R.Sp. 29.Apr. 1993)
          IF(KZ(33).GT.1)JHISTR(NREG) = JHISTR(NREG) + 1
*
          call cputim(tt1)
*



*
      DO 801 L = 1,NREG
          I = IREG(L)
*
          DO 655 K = 1,LMAX
 655 IMPI(K,L) = LIST(K,I)
*
          NBSUM = 0
*
       CALL REGINT(I,IMPI(1,L))
*
         DO 615 K = 1,3
              F(K,I) = 0.5D0*(FI(K,I) + FR(K,I))
              FDOT(K,I) = ONE6*(FIDOT(K,I) + FRDOT(K,I))
  615 CONTINUE
  801 CONTINUE
*
      call cputim(tt2)
      ttreg = ttreg + (tt2-tt1)*60.
# 768 "intgrt.F"
*
      NSTEPR = NSTEPR + NREG
      NBLCKR = NBLCKR + 1
      LSHRNK = .FALSE.
*
      DO 81 L = 1,NREG
          I = IREG(L)
*
              DO 816 K = 1,LMAX
 816 LIST(K,I) = IMPI(K,L)
*
* Check minimum neighbor sphere since last output
              IF(LIST(1,I).GT.0)RSMIN = MIN(RSMIN,RS(I))
*
 81 CONTINUE
*
* OPEN(98,STATUS='OLD',ERR=123)
* print*,' last reg block t=',time,' length=',nreg
* print*,' first 10 =',(name(ireg(l)),l=1,min(nreg,10))
* call flush(6)
* CLOSE(98)
*123 CONTINUE
*
      END IF
*
* Copy all corrected coordinates & velocities (NB! only at the end).
      DO 85 L = 1,NXTLEN
          I = NXTLST(L)
*
      IF (I.GT.N) THEN
          IPAIR = I - N
          IF (LIST(1,2*IPAIR-1).GT.0) NSTEPB = NSTEPB + 1
      END IF
*
          IF(TIMENW(I).LT.TMIN)THEN
              TMIN = TIMENW(I)
              IMIN = I
          END IF
*
          DO 82 K = 1,3
              X0(K,I) = XN(K,I)
              X0DOT(K,I) = XNDOT(K,I)
              X(K,I) = XN(K,I)
              XDOT(K,I) = XNDOT(K,I)
   82 CONTINUE
   85 CONTINUE
*
* Exit on KS termination, new multiple regularization or merger.
      IF (IQ.GT.0) THEN
          NBPREV = 0
          IF (IQ.GE.4.AND.IQ.NE.7) THEN
              CALL DELAY(IQ,-1)
          ELSE
* Ensure correct KS index (KSPAIR may denote second termination).
              KSPAIR = KVEC(I10)
              IPHASE = IQ
          END IF
          GO TO 100
      END IF
*
* Perform optional check on high-velocity particles at major times.
      IF (KZ(18).GT.0.AND.LISTV(1).GT.0) THEN
          IF (DMOD(TIME,STEPM).EQ.0.0D0) THEN
              CALL SHRINK
              IF (LISTV(1).GT.0) THEN
                  CALL HIVEL(-1)
              END IF
          END IF
      END IF
*
* Check optional mass loss time.
      IF (KZ(19).GT.0) THEN
* Delay until time commensurate with 1000-year step (new polynomials).
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
* Advance counters and check timer & optional COMMON save (NSUB = 0).
      NTIMER = NTIMER + NXTLEN
      IF (NTIMER.LT.NMAX) GO TO 1
      NTIMER = 0
      NSTEPS = NSTEPS + NMAX
*
      IF (NSTEPS.GE.100*NMAX.AND.NSUB.EQ.0) THEN
          NSTEPS = 0
          IF (KZ(1).GT.1) CALL MYDUMP(1,1)
      END IF
*
* Check option for general binary search.
      IF (KZ(4).NE.0.AND.TIME - TLASTS.GT.DELTAS) THEN
          CALL EVOLVE(0,0)
      END IF
*
* Include facility for termination of run (create dummy file STOP).
      IF(rank.EQ.0)THEN
      OPEN (99,FILE='STOP',STATUS='OLD',FORM='FORMATTED',IOSTAT=IO)
      IF (IO.EQ.0) THEN
          CLOSE (99)
          IF (NSUM.EQ.0.and.rank.eq.0) WRITE (6,90)
   90 FORMAT (/,9X,'TERMINATION BY MANUAL INTERVENTION')
          CPU = 0.0
      END IF
      END IF
*
* Repeat cycle until elapsed computing time exceeds the limit.
      CALL CPUTIM(TCOMP)
          TCOMP = (TCOMP-TTOTA)*60.
*
      IF (TCOMP.LT.CPU) GO TO 1
*
* Do not terminate during triple, quad or chain regularization.
      IF (NSUB.GT.0) THEN
* Specify zero step to enforce termination.
          DO 95 L = 1,NSUB
              STEPS(L) = 0.0D0
   95 CONTINUE
          NTIMER = NMAX
          GO TO 1
      END IF
*
* Terminate run with optional COMMON save.
      IF (KZ(1).NE.0) THEN
          CPUTOT = CPUTOT + TCOMP - CPU0
          CALL MYDUMP(1,1)
          if(rank.eq.0)
     & WRITE (6,98) TOFF, TIME, TIME+TOFF, TCOMP, CPUTOT/60.0,
     & ERRTOT, DETOT
   98 FORMAT (
     & '  TCOMP =',F7.1,'  CPUTOT =',F6.1,
     & '  ERRTOT =',F10.6,'  DETOT =',F10.6)
      END IF
*
* Determine time interval and step numbers per time unit
      TIMINT = TIME + TOFF - TINIT
*



      WRITE (6,195) rank,TIMINT,NSTEPI-NIR,NSTEPB-NIB,NSTEPR-NRGL,
     & NSTEPU-NKS
  195 FORMAT (
     & ' NIRRB=',I11,' NREG=',I11,' NKS=',I11)
      WRITE (6,196) (NSTEPI-NIR)/TIMINT,(NSTEPB-NIB)/TIMINT,
     & (NSTEPR-NRGL)/TIMINT,(NSTEPU-NKS)/TIMINT
  196 FORMAT (
     & D12.5,' NREG=',D12.5,' NKS=',D12.5)







          STOP
*
 100 CONTINUE
*
* Set current global time.
          TTOT = TIME + TOFF
* Full prediction at end of intgrt to preserve consistency after
* parallel execution
          DO 400 J = IFIRST,NTOT
              S = TIME - T0(J)
              S1 = 1.5*S
              S2 = 2.0*S
              X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S +X0(1,J)
              X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S +X0(2,J)
              X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S +X0(3,J)
              XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
              XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
              XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
 400 CONTINUE
* Resolve any KS coordinates & velocities using most recent c.m.
      IF (NPAIRS.GT.0) THEN
          JJ = -1
          DO 88 JPAIR = 1,NPAIRS
          JJ = JJ + 2
          IF (LIST(1,JJ).GT.0) THEN
              ZZ = 1.0
              IZZ = -3
* Distinguish between low and high-order prediction of U & UDOT.
              IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
              CALL KSRES2(JPAIR,J1,J2,ZZ,IZZ)
          END IF
   88 CONTINUE
      END IF
*
       RETURN
*
      END
