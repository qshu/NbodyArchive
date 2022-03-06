      SUBROUTINE INTGRT
*
*
*       N-body integrator flow control.
*       -------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      INTEGER JHIST,JHISTR
      COMMON/BLKLVL/  JHIST(NMAX),JHISTR(NMAX)
      INTEGER  NXTLST(NMAX),IREG(NMAX),NBLIST(NMAX)
      logical llpr
*
      INTEGER IKSRL,JKSRL
      COMMON/LOCKSR/IKSRL(NMAX),JKSRL(NMAX)
c
#ifdef PARALLEL
      integer inum(maxpe),ista(maxpe)
      REAL*8 XMPI(19,NMAX),YMPI(41,NMAX)
#endif
      INTEGER IMPI(LMAX,NMAX),JMPI(10,maxpe)
*
*       Reset control & regularization indicators.
  999 IPHASE = 0
      IKS = 0
*       Initialize end-point of integration times.
      DO 1000 I = IFIRST,NTOT
         TIMENW(I) = T0(I) + STEP(I)
 1000 CONTINUE
      TMIN = TPREV
*
*       Find all particles due at next block time.
    1 CONTINUE
*
      NXTLEN = 0
*
*     determine next block particles without assuming
*     sorted time step list (R.Sp.)
*
      TMIN = TIMENW(IFIRST)
      IMIN = IFIRST
*
      DO 3 J = IFIRST+1, NTOT
         IF(TIMENW(J).LT.TMIN)THEN
         TMIN = TIMENW(J)
         IMIN = J
         END IF
  3   CONTINUE
*
      DO 5 J = IFIRST, NTOT
         IF(DABS(TIMENW(J)-TMIN).LT.DTK(64)) THEN
            NXTLEN = NXTLEN + 1
            NXTLST(NXTLEN) = J
         END IF
  5   CONTINUE
*
*       Set new time and save block time (for regularization terminations).
      TIME = T0(IMIN) + STEP(IMIN)
      TBLOCK = TIME
*
      IF (STEP(IMIN).LT.1.0E-19.OR.DMOD(TIME,STEP(IMIN)).NE.0.0D0) THEN
          WRITE (6,1005) IMIN, NAME(IMIN), NSTEPI, TIME, STEP(IMIN), 
     &                   TIME/STEP(IMIN)
 1005     FORMAT (' DANGER!   I NM # TIME STEP T/DT ',
     &                        2I5,I11,F12.5,1P,E9.1,0P,F16.4)
          CALL ABORT
      END IF
*
*       Check for new regularization at end of block.
      IF (IKS.GT.0) THEN
          TIME = TPREV
          IPHASE = 1
          GO TO 100
      END IF
*
*       Check next adjust time before beginning a new block.
      IF (TIME.GT.TADJ) THEN
          TIME = TPREV
          IPHASE = 3
          GO TO 100
      END IF
*
*       Also check output time in case DTADJ & DELTAT not commensurate.
      IF (TIME.GT.TNEXT) THEN
          TIME = TPREV
          CALL OUTPUT
          GO TO 1
      END IF
*
*       See whether to advance any KS solutions first.
      IF (NPAIRS.GT.0) THEN
      call cputim(tt1)
   10     TKSMIN = 1.0E+20
          DO 15 JPAIR = 1,NPAIRS
              J1 = 2*JPAIR - 1
              TKS = T0(J1) + STEP(J1)
              IF (TKS.LE.TKSMIN) THEN
                  IPAIR = JPAIR
                  TKSMIN = TKS
              END IF
   15     CONTINUE
*
         IIKS = 0
*       See whether the smallest KS time falls due before next block step.
          IF (TKSMIN.LE.TBLOCK) THEN
              TIME = TKSMIN
              I1 = 2*IPAIR - 1
              CALL KSINT(I1)
*       Check for termination, multiple regularization or merger.
              IF (IPHASE.NE.0) THEN
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
                  GO TO 100
              END IF
              GO TO 10
          END IF
          TIME = TBLOCK
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
      END IF
*
*       Check time for advancing any triple, quad or chain regularization.
      IF (NSUB.GT.0) THEN
   20     TSUB = 1.0D+10
          DO 25 L = 1,NSUB
              IF (TS(L).LT.TSUB) THEN
                  ISUB = L
                  TSUB = TS(L)
              END IF
   25     CONTINUE
*
*       See whether the smallest subsystem time is due before next block.
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
*       Specify maximum integration interval for chain (unless TCOMP > CPU).
                  IF (STEPS(ISUB).GT.0.0D0) THEN
                      STEPS(ISUB) = TBLOCK - TS(ISUB)
                  END IF
                  CALL CHAIN(ISUB)
              END IF
*       Resort the block after termination or membership change.
              IF (ISUB.LT.0.OR.IPHASE.LT.0) THEN
                  GO TO 999
              END IF
              GO TO 20
          END IF
          TIME = TBLOCK
      END IF
*
*       Check regularization criterion for single particles.
      IKS = 0
      ISMIN = 0
      DSMIN = DTMIN
      DO 50 L = 1,NXTLEN
          I = NXTLST(L)
*       Search for minimum timestep candidate for not ordered steps (R.Sp.)
          IF (STEP(I).LT.DTMIN.AND.STEP(I).LT.DSMIN.AND.I.LE.N) THEN
              DSMIN = STEP(I)
              ISMIN = I
          END IF
   50 CONTINUE
*
*       See whether dominant body can be regularized.
              IF(ISMIN.GT.0)CALL SEARCH(ISMIN,IKS)
*
*       Choose between predicting all neighbours or full N.
*       Warning do not distribute prediction on PE's for consistency
          call cputim(tt1)
      IF (NXTLEN.LT.1) THEN
*
*       Merge all neighbour lists (with absent members of NXTLST added).
          CALL NBSORT(NXTLEN,NXTLST,NNB,NBLIST)
*
*       Predict coordinates & velocities of neighbours and #I to order FDOT.
          NBPRED = NBPRED + NNB
          NBFLAG = 1
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
   35     CONTINUE
      ELSE
          NNPRED = NNPRED + 1
          NBFLAG = 1
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
 40       CONTINUE
      END IF
*
*       Resolve any KS coordinates & velocities using most recent c.m.
      IF (NPAIRS.GT.0) THEN
          DO 45 JPAIR = 1,NPAIRS
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
              CALL KSRES2(JPAIR,J1,J2,0.0)
   45     CONTINUE
      END IF
*
		  call cputim(tt2)
          ttnbp = ttnbp + (tt2-tt1)*60.
*
*       Initialize counters for irregular & regular integrations.
   55 CONTINUE
      NREG = 0
*
*       Advance the irregular step for all particles in the current block.
*       Block-Step Level Diagnostics (R.Sp. 29.Apr. 1993)
          IF(KZ(34).GT.0)JHIST(NXTLEN) = JHIST(NXTLEN) + 1
*
#ifdef PARALLEL
      IF(NXTLEN.LT.isernb)THEN
#endif
*
      DO 701 L = 1,NXTLEN
*
          I = NXTLST(L)
      CALL NBINT(I,NBFLAG)
*
 701    CONTINUE
*
#ifdef PARALLEL
*       start PARALLEL section for nbint
        ELSE
*
      nl = NXTLEN
*
      inl = nl/isize
      jsize = isize*inl
      idiff = nl - jsize
      irun = 0
*
      do 1003 ix = 1,isize
      inum(ix)=inl
      if(ix.le.idiff)inum(ix) = inum(ix) + 1
      ista(ix) = irun+1
      if(ista(ix).gt.nl)inum(ix) = 0
 1003 irun = irun + inum(ix)
*
      istart = ista(rank+1)
      iend = ista(rank+1) + inum(rank+1) - 1
*
*     if(time.lt.0.1d0)then
*     print*,' rank ',rank,' NXTLEN ',NXTLEN,TIME,
*    *    ' istart,iend=',istart,iend
*     end if
*
          call mpi_barrier(MPI_COMM_WORLD,ierr)
      DO 70 L = istart,iend
*
*     if(time.lt.0.1d0)then
*     print*,' rank ',rank,' NXT in DO L ',L,TIME,
*    *    ' istart,iend=',istart,iend
*     end if
          I = NXTLST(L)
      CALL NBINT(I,NBFLAG)
*
          DO 69 K = 1,3
          XMPI(K,L) = XN(K,I)
          XMPI(K+3,L) = XNDOT(K,I)
          XMPI(K+6,L) = D2(K,I)
          XMPI(K+9,L) = D3(K,I)
          XMPI(K+12,L) = FI(K,I)
          XMPI(K+15,L) = FIDOT(K,I)
 69       CONTINUE

          XMPI(19,L) = STEP(I)
*
 70    CONTINUE
*
          call cputim(tt2)
          ttirr = ttirr + (tt2-tt1)*60.
*
*        Distribute variables into private vectors again T3E (R.Sp.)
*
      isend = rank + 1
      if(isend.eq.isize)isend = 0
      irecv = rank - 1
      if(irecv.eq.-1)irecv = isize - 1
*
      do 1001 ir = 0,isize-2
*
      irank = rank - ir
      if(irank.lt.0)irank=irank+isize
*
      istsen=ista(irank+1)
      icnt = inum(irank+1)
*
      if(irank.eq.0)irank=isize
      istrec = ista(irank)
      icnt2 = inum(irank)
*
*     if(time.lt.0.1d0.and.icnt.gt.0)then
*     print*,' NXT: rank t',rank,time,' ir ',ir,' send ',istsen,
*    *    ' thru ',istsen+icnt-1,' to ',isend,' cnt ',icnt,
*    *    ' istart,iend=',istart,iend
*     end if
*     if(time.lt.0.1d0.and.icnt2.gt.0)then
*     print*,' NXT: rank t',rank,time,' ir ',ir,' recv ',istrec,
*    *    ' thru ',istrec+icnt2-1,' fr ',irecv,' cnt2 ',icnt2,
*    *    ' istart,iend=',istart,iend
*     end if
*
#ifdef PUREMPI
          call mpi_barrier(MPI_COMM_WORLD,ierr)
      CALL MPI_SENDRECV(XMPI(1,istsen),19*icnt,MPI_REAL,isend,rank,
     *                  XMPI(1,istrec),19*icnt2,MPI_REAL,irecv,irecv,
     *                  MPI_COMM_WORLD,status,ierr)
#endif
*
#ifdef SHMEM
       call barrier()
       call shmem_get(XMPI(1,istrec),XMPI(1,istrec),19*icnt2,irecv)
#endif
*
 1001  continue
*
          call mpi_barrier(MPI_COMM_WORLD,ierr)
       DO 72 L = 1,NXTLEN
          I = NXTLST(L)
              DO 68 K = 1,3
                  XN(K,I) = XMPI(K,L)
                  XNDOT(K,I) = XMPI(K+3,L)
                  D2(K,I) = XMPI(K+6,L)
                  D3(K,I) = XMPI(K+9,L)
                  FI(K,I) = XMPI(K+12,L)
                  FIDOT(K,I) = XMPI(K+15,L)
   68         CONTINUE
          STEP(I) = XMPI(19,L)
  72      CONTINUE
*
            call cputim(tt3)
            ttmov = ttmov + (tt3-tt2)*60.
*
       END IF
*          End PARALLEL section for nbint
#endif
*
       DO 71 L = 1,NXTLEN
          I = NXTLST(L)
*       Save new block step and update T0 & next time
          T0(I) = TIME
          TIMENW(I) = T0(I) + STEP(I)
*
*       Set non-zero indicator for new regular force.
          IF (T0R(I) + STEPR(I).LE.TIME) THEN
              NREG = NREG + 1
              IREG(NREG) = I
          ELSE
*       Extrapolate regular force & first derivatives to obtain F & FDOT.
              DTR = TIME - T0R(I)
              DO 65 K = 1,3
                  F(K,I) = 0.5*(FRDOT(K,I)*DTR + FR(K,I) + FI(K,I))
                  FDOT(K,I) = ONE6*(FRDOT(K,I) + FIDOT(K,I))
* Higher order extrapolation?
*                 F(K,I) = FI(K,I) + FR(K,I) + DTR*(FRDOT(K,I)
*    *                + DTR*(D2R(K,I)/2.D0 + DTR*D3R(K,I)/6.D0))
*                 FDOT(K,I) = FIDOT(K,I) + FRDOT(K,I)
*    *                + DTR*(D2R(K,I) + DTR*D3R(K,I)/2.D0)
*                 F(K,I) = F(K,I)/2.D0
*                 FDOT(K,I) = FDOT(K,I)/6.D0
   65         CONTINUE
          END IF
*
              DO 67 K = 1,3
                  X0(K,I) = XN(K,I)
                  X0DOT(K,I) = XNDOT(K,I)
                  D0(K,I) = FI(K,I)
                  D1(K,I) = FIDOT(K,I)
   67         CONTINUE
*
 71     CONTINUE
*
*           CALL nemo_savestate(n,3,time,body,x,xdot)
*
*       See whether any KS candidates are in the same block.
      IF (IKS.GT.0) THEN
*       Accept same time, otherwise reduce STEP(ICOMP) and/or delay.
          IF (T0(JCOMP).EQ.T0(ICOMP)) THEN
              I = ICOMP
              ICOMP = MIN(ICOMP,JCOMP)
              JCOMP = MAX(I,JCOMP)
          ELSE IF (T0(JCOMP) + STEP(JCOMP).LT.
     &             T0(ICOMP) + STEP(ICOMP)) THEN
              STEP(ICOMP) = 0.5D0*STEP(ICOMP)
              TIMENW(ICOMP) = T0(ICOMP) + STEP(ICOMP)
              IKS = 0
          ELSE
              IKS = 0
          END IF
      END IF
*
*       Update the previous time and increase counters (irregular & block).
      TPREV = TIME
      NCOUNT(1) = NCOUNT(1) + NXTLEN
      NSTEPI = NSTEPI + NXTLEN
      NBLOCK = NBLOCK + 1
*
*       Obtain total force for all particles due in the current block.
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
 811      CONTINUE
*       Resolve any KS coordinates & velocities using most recent c.m.
      IF (NPAIRS.GT.0) THEN
          DO 86 JPAIR = 1,NPAIRS
              CALL KSRES2(JPAIR,J1,J2,0.0)
   86     CONTINUE
      END IF
*
          call cputim(tt8)
          ttpre = ttpre + (tt8-tt7)*60.
*
*       Block-Step Level Diagnostics (R.Sp. 29.Apr. 1993)
          IF(KZ(34).GT.1)JHISTR(NREG) = JHISTR(NREG) + 1
*
#ifdef PARALLEL
      IF(NREG.LT.iserreg)THEN
#endif
          call cputim(tt1)
*
      DO 801 L = 1,NREG
          I = IREG(L)
*
          DO 655 K = 1,LMAX
 655      IMPI(K,L) = LIST(K,I)
*
       CALL REGINT(I,IMPI(1,L))
*
         DO 615 K = 1,3
              F(K,I) = 0.5D0*(FI(K,I) + FR(K,I))
              FDOT(K,I) = ONE6*(FIDOT(K,I) + FRDOT(K,I))
  615    CONTINUE
  801 CONTINUE
*
      call cputim(tt2)
      ttreg = ttreg + (tt2-tt1)*60.
#ifdef PARALLEL
*       Start PARALLEL section for regint
      ELSE
*
	  call cputim(tt1)
      nl = NREG
*
      inl = nl/isize
      jsize = isize*inl
      idiff = nl - jsize
      irun = 0
*
      do 1103 ix = 1,isize
      inum(ix)=inl
      if(ix.le.idiff)inum(ix) = inum(ix) + 1
      ista(ix) = irun+1
      if(ista(ix).gt.nl)inum(ix) = 0
 1103 irun = irun + inum(ix)
*
      istart = ista(rank+1)
      iend = ista(rank+1) + inum(rank+1) - 1
*     PRINT*,' Start reg block rank ',rank,
*    * ' TIME=',TIME,' NREG=',NREG
*     PRINT*,' Start reg block rank ',rank,
*    * ' istart,iend=',istart,iend
*
          NC5 = NBFULL
          NC6 = NBVOID
          NC11 = NBSMIN
          NC15 = NBDIS2
          NC18 = NRCONV
          NC19 = NBDIS
          NC20 = NLSMIN
          NC30 = NFAST
          NCX = NICONV
          NCY = NBFAST
*
      DO 80 L = istart,iend
          I = IREG(L)
*
          DO 656 K = 1,LMAX
 656      IMPI(K,L) = LIST(K,I)
*
       CALL REGINT(I,IMPI(1,L))
*
          YMPI(37,L) = STEPR(I)
          YMPI(38,L) = STEP(I)
          YMPI(39,L) = TIMENW(I)
          YMPI(40,L) = RS(I)
          YMPI(41,L) = T0R(I)
*
         DO 814 K = 1,3
              YMPI(K,L) = X0(K,I)
              YMPI(K+3,L) = X0DOT(K,I)
              YMPI(K+6,L) = FI(K,I)
              YMPI(K+9,L) = FR(K,I)
              YMPI(K+12,L) = FIDOT(K,I)
              YMPI(K+15,L) = FRDOT(K,I)
              YMPI(K+18,L) = D1R(K,I)
              YMPI(K+21,L) = D2R(K,I)
              YMPI(K+24,L) = D3R(K,I)
*       Save cXrrections of force polynomials from routine fpcorr
              IF (KZ(38).GT.0.OR.I.GT.N) THEN
              YMPI(K+27,L) = D1(K,I)
              YMPI(K+30,L) = D2(K,I)
              YMPI(K+33,L) = D3(K,I)
              END IF
  814    CONTINUE

   80 CONTINUE
*
          j = rank + 1
          JMPI(1,j) = NBFULL - NC5
          JMPI(2,j) = NBVOID - NC6
          JMPI(3,j) = NBSMIN - NC11
          JMPI(4,j) = NBDIS2 - NC15
          JMPI(5,j) = NRCONV - NC18
          JMPI(6,j) = NBDIS - NC19
          JMPI(7,j) = NLSMIN - NC20
          JMPI(8,j) = NFAST - NC30
          JMPI(9,j) = NICONV - NCX
          JMPI(10,j) = NBFAST - NCY
*
          call cputim(tt2)
          ttreg = ttreg + (tt2-tt1)*60.
*
*        Distribute variables into private vectors again T3E (R.Sp.)
      isend = rank + 1
      if(isend.eq.isize)isend = 0
      irecv = rank - 1
      if(irecv.eq.-1)irecv = isize - 1
*
      do 1002 ir = 0,isize-2
*
      irank = rank - ir
      if(irank.lt.0)irank=irank+isize
*
      istsen=ista(irank+1)
      icnt = inum(irank+1)
*
      if(irank.eq.0)irank=isize
      istrec = ista(irank)
      icnt2 = inum(irank)
*
*     print*,' INTGRT-R: rank t',rank,time,' ir ',ir,' send ',istsen,
*    *    ' thru ',istsen+icnt-1,' to ',isend,' cnt ',icnt
*     print*,' INTGRT-R: rank t',rank,time,' ir ',ir,' recv ',istrec,
*    *    ' thru ',istrec+icnt2-1,' fr ',irecv,' cnt2 ',icnt2
*
#ifdef PUREMPI
          call mpi_barrier(MPI_COMM_WORLD,ierr)
      CALL MPI_SENDRECV(YMPI(1,istsen),41*icnt,MPI_REAL,isend,
     *         rank,YMPI(1,istrec),41*icnt2,MPI_REAL,irecv,irecv,
     *         MPI_COMM_WORLD,status,ierr)
*
      CALL MPI_SENDRECV(IMPI(1,istsen),lmax*icnt,MPI_INTEGER,isend,
     *         rank,IMPI(1,istrec),lmax*icnt2,MPI_INTEGER,irecv,irecv,
     *         MPI_COMM_WORLD,status,ierr)
#endif
*
#ifdef SHMEM
*
       call barrier()
       call shmem_get(YMPI(1,istrec),YMPI(1,istrec),41*icnt2,irecv)
       call shmem_get(IMPI(1,istrec),IMPI(1,istrec),lmax*icnt2,irecv)
#endif
*
 1002  continue
*
      CALL MPI_BCAST(JMPI(1,rank+1),10,MPI_INTEGER,0,
     *                                  MPI_COMM_WORLD,ierr)
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      DO 812 L = 1,NREG
          I = IREG(L)
         DO 815 K = 1,3
              X0(K,I) = YMPI(K,L)
              X0DOT(K,I) = YMPI(K+3,L)
              FI(K,I) = YMPI(K+6,L)
              FR(K,I) = YMPI(K+9,L)
              F(K,I) = 0.5D0*(FI(K,I) + FR(K,I))
              FIDOT(K,I) = YMPI(K+12,L)
              FRDOT(K,I) = YMPI(K+15,L)
              FDOT(K,I) = ONE6*(FIDOT(K,I) + FRDOT(K,I))
              D0(K,I) = FI(K,I)
              D0R(K,I) = FR(K,I)
              D1R(K,I) = YMPI(K+18,L)
              D2R(K,I) = YMPI(K+21,L)
              D3R(K,I) = YMPI(K+24,L)
*       Save corrections of force polynomials from routine fpcorr
              IF (KZ(38).GT.0.OR.I.GT.N) THEN
              D1(K,I) = YMPI(K+27,L)
              D2(K,I) = YMPI(K+30,L)
              D3(K,I) = YMPI(K+33,L)
              END IF
  815    CONTINUE
              STEPR(I) = YMPI(37,L)
              STEP(I) = YMPI(38,L)
              TIMENW(I) = YMPI(39,L)
              RS(I) = YMPI(40,L)
              T0R(I) = YMPI(41,L)
 812   CONTINUE
*
	  call cputim(tt3)
      ttmov = ttmov + (tt3-tt2)*60.
*
	  call cputim(tt1)
       DO 83 J = 1,isize
              NCOUNT(5) = NC5 + JMPI(1,J)
              NBFULL = NC5 + JMPI(1,J)
              NCOUNT(6) = NC6 + JMPI(2,J)
              NBVOID = NC6 + JMPI(2,J)
              NCOUNT(11) = NC11 + JMPI(3,J)
              NBSMIN = NC11 + JMPI(3,J)
              NCOUNT(15) = NC15 + JMPI(4,J)
              NBDIS2 = NC15 + JMPI(4,J)
              NCOUNT(18) = NC18 + JMPI(5,J)
              NRCONV = NC18 + JMPI(5,J)
              NCOUNT(19) = NC19 + JMPI(6,J)
              NBDIS = NC19 + JMPI(6,J)
              NCOUNT(20) = NC20 + JMPI(7,J)
              NLSMIN = NC20 + JMPI(7,J)
              NCOUNT(30) = NC30 + JMPI(8,J)
              NFAST = NC30 + JMPI(8,J)
              NICONV = NCX + JMPI(9,J)
              NBFAST = NCY + JMPI(10,J)
*
   83 CONTINUE
      call cputim(tt2)
      ttsub = ttsub + (tt2-tt1)*60.
*
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      END IF
*         End PARALLEL section for regint
#endif
*
      NCOUNT(2) = NCOUNT(2) + NREG
      NSTEPR = NSTEPR + NREG
      NBLCKR = NBLCKR + 1
*
      DO 81 L = 1,NREG
          I = IREG(L)
*
              DO 816 K = 1,LMAX
 816          LIST(K,I) = IMPI(K,L)
*
*       Check minimum neighbor sphere since last output
              RSMIN = MIN(RSMIN,RS(I))
 81    CONTINUE
*
      END IF
*
*       Copy all corrected coordinates & velocities (NB! only at the end).
      DO 85 L = 1,NXTLEN
          I = NXTLST(L)
          DO 82 K = 1,3
              X(K,I) = X0(K,I)
              XDOT(K,I) = X0DOT(K,I)
   82     CONTINUE
   85 CONTINUE
*
*       Advance counters and check timer & optional COMMON save (NSUB = 0).
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
*       Check option for general binary search.
      IF (KZ(4).NE.0.AND.TIME - TLASTS.GT.DELTAS) THEN
          CALL EVOLVE(0,0)
      END IF
*
*       Include facility for termination of run (create dummy file STOP).
      IF(rank.EQ.0)THEN
      OPEN (99,FILE='STOP',STATUS='OLD',FORM='FORMATTED',IOSTAT=IO)
      IF (IO.EQ.0) THEN
          CLOSE (99)
          IF (NSUM.EQ.0) WRITE (6,90)
   90     FORMAT  (/,9X,'TERMINATION BY MANUAL INTERVENTION')
          CPU = 0.0
      END IF
      END IF
*
*       Repeat cycle until elapsed computing time exceeds the limit.
      CALL CPUTIM(TCOMP)
      IF (TCOMP.LT.CPU) GO TO 1
*
*       Do not terminate during triple, quad or chain regularization.
      IF (NSUB.GT.0) THEN
*       Specify zero step to enforce termination.
          DO 95 L = 1,NSUB
              STEPS(ISUB) = 0.0D0
   95     CONTINUE
          NTIMER = NMAX
          GO TO 1
      END IF
*
*       Terminate run with optional COMMON save.
      IF (KZ(1).NE.0) THEN
          CPUTOT = CPUTOT + TCOMP - CPU0
          CALL MYDUMP(1,1)
          WRITE (6,98)  TIME, TCOMP, CPUTOT/60.0, ERRTOT, DETOT
   98     FORMAT (//,9X,'COMMON SAVED AT TIME =',F8.2,'  TCOMP =',F7.1,
     &                  '  CPUTOT =',F6.1,'  ERRTOT =',F10.6,
     &                  '  DETOT =',F10.6)
      END IF
*
          STOP
*
 100      CONTINUE
*          Full prediction at end of intgrt to preserve consistency after
*          parallel execution
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
 400      CONTINUE
*       Resolve any KS coordinates & velocities using most recent c.m.
      IF (NPAIRS.GT.0) THEN
          DO 88 JPAIR = 1,NPAIRS
              CALL KSRES2(JPAIR,J1,J2,0.0)
   88     CONTINUE
      END IF
*
       RETURN
*
      END
