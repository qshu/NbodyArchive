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
      LOGICAL LSHRNK
      EXTERNAL SHORT
      COMMON/BLKLVL/JHIST(0:NMAX),JHISTR(0:NMAX)
      INTEGER  NXTLST(NMAX),IREG(NMAX),NBLIST(NMAX),IBL(LMAX)
      real*8 r2(nmax)
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
      COMMON/CHAINC/  DF1(3*NCMAX),DF2(3*NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      LOGICAL UPDATEF(0:MAXPE-1)
      COMMON/XCUPD/  UPDATEF
      COMMON/STSTAT/  TINIT,NIR,NIB,NRGL,NKS
*
#ifdef PARALLEL
      integer inum(maxpe),ista(maxpe)
      REAL*8 XMPI(19,NMAX),YMPI(41,NMAX)
#endif
      INTEGER IMPI(LMAX,NMAX),JMPI(11,maxpe)
*
      SAVE IQ,ICALL,STEPM
      DATA IQ,ICALL,STEPM /0,2,0.03125/
*
*       Search for high velocities after escape or KS/chain termination.
  999 IF (KZ(18).GT.0.AND.(IPHASE.EQ.-1.OR.IPHASE.GE.2)) THEN
          CALL HIVEL(0)
      END IF
*
*       Reset control & regularization indicators.
      IPHASE = 0
      IKS = 0
*     CALL ENERGY
*     ETOT = ZKIN-POT+EBIN
*     IF(TIME.EQ.0.D0)ZKIN0=ZKIN
*       Initialize end-point of integration times and set TMIN.
      TMIN = 1.0D+10
*       Initialize end-point of integration times and set TMIN and DTM.
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
*       Find all particles due at next block time.
    1 CONTINUE
*
      NXTLEN = 0
*
*       Reset TSMALL second time after main change to catch new small steps.
      ICALL = ICALL + 1
      IF (ICALL.EQ.2) GO TO 999
*
*     determine next block particles without assuming
*     sorted time step list (R.Sp.)
*
      DO 5 J = IFIRST, NTOT
         IF(DABS(TIMENW(J)-TMIN).LT.DTK(40)) THEN
            NXTLEN = NXTLEN + 1
            NXTLST(NXTLEN) = J
         END IF
  5   CONTINUE
*
*       Update short timestep list for regularization search.
      CALL SHORT(NXTLEN,NXTLST)
*
*       Set new time and save block time (for regularization terminations).
      TIME = TMIN
      TBLOCK = TIME
*
*       Check option for advancing interstellar clouds.
      IF (KZ(13).GT.0) THEN
          CALL CLINT
      END IF
*
*       Include commensurability test (may be suppressed if no problems).
*     IF (STEP(IMIN).LT.1.0E-15.OR.DMOD(TIME,STEP(IMIN)).NE.0.0D0) THEN
*         WRITE (6,1005) IMIN, NAME(IMIN), NSTEPI, TIME, STEP(IMIN), 
*    &                   TIME/STEP(IMIN)
*1005     FORMAT (' DANGER!   I NM # TIME STEP T/DT ',
*    &                        2I5,I11,F12.5,1P,E9.1,0P,F16.4)
*         CALL ABORT
*     END IF
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
*       See whether to advance any close encounters at first new time.
      IF (TIME.GT.TPREV) THEN
            call cputim(tt5)
          CALL SUBINT(IQ,I10)
            call cputim(tt6)
            ttks = ttks + (tt6-tt5)*60.
*
          IF (IQ.LT.0) GO TO 999
      END IF
*
*       Check regularization criterion for single particles.
      IKS = 0
      ISMIN = 0
      DSMIN = DTMIN
*       Search only in prepared list of short-step particles. (R.Sp.)
      ISHORT = LSHORT(1)
      DO 50 L = 2,ISHORT+1
          I = LSHORT(L)
*       Search for minimum timestep candidate for not ordered steps (R.Sp.)
*       Beware that members of LSHORT may be members of KS pair (R.Sp.)
          IF (STEP(I).LT.DTMIN.AND.STEP(I).LT.DSMIN.AND.I.LE.N.AND.
     &        I.GE.IFIRST) THEN
              DSMIN = STEP(I)
              ISMIN = I
          END IF
   50 CONTINUE
*
*       See whether dominant body can be regularized.
      IF(ISMIN.GT.0) THEN
          CALL SEARCH(ISMIN,IKS)
*
*       Include close encounter search for low-eccentric massive binaries.
      IF (IKS.EQ.0.AND.STEP(ISMIN).LT.4.0*DTMIN) THEN
*       Consider massive single bodies in absence of subsystems.
          IF (ISMIN.LE.N.AND.BODY(I).GT.2.0*BODYM.AND.NSUB.EQ.0) THEN
*
*       Obtain two-body elements and relative perturbation.
              JMIN = 0
              CALL ORBIT(ISMIN,JMIN,SEMI,ECC,GI)
*
              EB = -0.5*BODY(ISMIN)*BODY(JMIN)/SEMI
              IF (EB.LT.EBH.AND.GI.LT.0.25.AND.JMIN.GE.IFIRST) THEN
                  APO = SEMI*(1.0 + ECC)
*       Check eccentricity (cf. max perturbation) and neighbour radius.
                  IF (ECC.LT.0.5.AND.APO.LT.0.02*RS(ISMIN)) THEN
*                     PRINT*, ' KS TRY: NAM E A EB ',
*    *       NAME(ISMIN), NAME(JMIN), ECC, SEMI, EB
*            CALL FLUSH(6)
                      IKS = IKS + 1
                      ICOMP = ISMIN
                      JCOMP = JMIN
                  END IF
              END IF
          END IF
      END IF
      END IF
*
*         OPEN(97,STATUS='OLD',ERR=124)
*         print*,' last irr block t=',time,' l=',nxtlen,
*    *    ' n=',n,' npairs=',npairs
*         ix=min(nxtlen,10)
*         print*,' first parts. in block (max 10) =',
*    *      (name(nxtlst(l)),l=1,ix)
*         if(ix.le.3)then
*         do 873 ii=1,ix
*         ik = nxtlst(ii)
*         print*,' i,name,stepi,r=',ik,name(ik),step(ik),stepr(ik)
*         if(ik.gt.n)then
*         ipair=ik-n
*         i1=2*ipair-1
*         i2=2*ipair
*         print*,' c.m.',ik,' w.members ',i1,i2,' names ',
*    *     name(i1),name(i2),' masses ',body(i1),body(i2)
*         end if
*873      continue
*         do 876 iy=1,ix
*         ipart = nxtlst(iy)
*         nny=list(1,ipart)
*         if(nny.gt.0)then
*         do 875 iz=1,nny
*         innb = list(1+iz,ipart)
*         r2(iz) = (x(1,ipart) - x(1,innb))**2 +
*    *     (x(2,ipart) - x(2,innb))**2 + (x(3,ipart) - x(3,innb))**2
*         jlist(iz) = innb
*875      continue
*         call sort1(nny,r2,jlist)
*         print*,' part ',ipart,' name=',name(ipart),
*    *     ' nnb=',nny,'-------------'
*         do 874 k=1,min(nny,3)
*         print*,' part, name, dist, step, mass=',jlist(k),
*    *    name(jlist(k)),r2(k),step(jlist(k)),body(jlist(k))
*         if(jlist(k).gt.n)then
*         ipair=jlist(k)-n
*         i1=2*ipair-1
*         i2=2*ipair
*         print*,'       w.members ',name(i1),name(i2),
*    *     ' masses ',body(i1),body(i2)
*         end if
*874      continue
*         end if
*         print*,'----------------------------------------'
*876      continue
*         end if
*         call flush(6)
*         CLOSE(97)
*124      CONTINUE
*
*       Choose between predicting all neighbours or full N.
*       Warning do not distribute prediction on PE's for consistency
          call cputim(tt1)
      IF (NXTLEN.LT.10) THEN
*
*       Initialize pointers for neighbour lists.
          DO 30 L = 1,NXTLEN
              IBL(L) = NXTLST(L)
   30     CONTINUE
*
*       Merge all neighbour lists (with absent members of NXTLST added).
          CALL NBSORT(NXTLEN,IBL,NNB,NBLIST)
*
*       Predict coordinates & velocities of neighbours and #I to order FDOT.
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
   35     CONTINUE
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
 40       CONTINUE
      END IF
*
*       Resolve any KS coordinates & velocities using most recent c.m.
      IF (NPAIRS.GT.0) THEN
*       Resolve perturbed KS pairs with c.m. prediction after NBSORT.
      JJ = -1
      DO 45 JPAIR = 1,NPAIRS
      JJ = JJ + 2
      IF (LIST(1,JJ).GT.0) THEN
*       Ignore c.m. prediction after full N loop (all active KS needed).
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
          ZZ = 1.0
          IZZ = -1
*       Distinguish between low and high-order prediction of U & UDOT.
          IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
          CALL KSRES2(JPAIR,J1,J2,ZZ,IZZ)
      END IF
   45 CONTINUE
*
      END IF
*
*       Save new time (output time at TIME> TADJ) and increase # of blocks.
      TPREV = TIME
      NBLOCK = NBLOCK + 1
      TMIN = 1.0D+10
*
		  call cputim(tt2)
          ttnbp = ttnbp + (tt2-tt1)*60.
*
*       Initialize counters for irregular & regular integrations.
      NREG = 0
*
*       Advance the irregular step for all particles in the current block.
*       Block-Step Level Diagnostics (R.Sp. 29.Apr. 1993)
          IF(KZ(33).GT.0)JHIST(NXTLEN) = JHIST(NXTLEN) + 1
*
#ifdef PARALLEL
      IF(NXTLEN.LE.isernb)THEN
#endif
*
      DO 701 L = 1,NXTLEN
*
          I = NXTLST(L)
*
      CALL NBINT(I,NBFLAG)
*             PRINT*,' Irr Step I,STEP,T0,T0+STEP=',
*    *       I,STEP(I),T0(I),T0(I)+STEP(I)
*            IF(I.GT.N)PRINT*,' I1,STEPU=',2*(I-N)-1,STEP(2*(I-N)-1)
*
 701    CONTINUE
*
#ifdef PARALLEL
*       start PARALLEL section for nbint
        ELSE
*
          call mpi_barrier(MPI_COMM_WORLD,ierr)
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
*
*        Reset updatef to check Chain updates in parallel part.
      updatef(rank)=.false.                         
*
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
*
*       Distribute updatef among all PE's
      CALL MPI_ALLGATHER(updatef(rank),1,MPI_LOGICAL,
     *  updatef(0),1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
*       Has any PE updated Chain-Variables ?
      do iy=0,isize-1
*
*         If so, does any PE have not ?
      if(updatef(iy)) then
          call mpi_barrier(MPI_COMM_WORLD,ierr)
*         If so, distribute Chain variables
       call MPI_BCAST(DF2(1),3*NCMAX,MPI_REAL,iy,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DF1(1),3*NCMAX,MPI_REAL,iy,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(LISTC(1),LMAX,MPI_INTEGER,iy,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(GPERT,1,MPI_REAL,iy,MPI_COMM_WORLD,ierr) 
*
       end if
       end do
*
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
*             PRINT*,' Reg Step I,STEPR,T0R,T0R+STEPR,STEP=',
*    *       I,STEPR(I),T0R(I),T0R(I)+STEPR(I),STEP(I)
*            IF(I.GT.N)PRINT*,' I1,STEPU=',2*(I-N)-1,STEP(2*(I-N)-1)
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
      NSTEPI = NSTEPI + NXTLEN
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
          JJ = -1
          DO 86 JPAIR = 1,NPAIRS
          JJ = JJ + 2
          IF (LIST(1,JJ).GT.0) THEN
              ZZ = 1.0
              IZZ = -2
*       Distinguish between low and high-order prediction of U & UDOT.
              IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
              CALL KSRES2(JPAIR,J1,J2,ZZ,IZZ)
          END IF
   86     CONTINUE
      END IF
*
          call cputim(tt8)
          ttpre = ttpre + (tt8-tt7)*60.
*
*       Block-Step Level Diagnostics (R.Sp. 29.Apr. 1993)
          IF(KZ(33).GT.1)JHISTR(NREG) = JHISTR(NREG) + 1
*
#ifdef PARALLEL
      IF(NREG.LE.iserreg)THEN
#endif
          call cputim(tt1)
*
      DO 801 L = 1,NREG
          I = IREG(L)
*
          DO 655 K = 1,LMAX
 655      IMPI(K,L) = LIST(K,I)
*
          NBSUM = 0
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
          call mpi_barrier(MPI_COMM_WORLD,ierr)
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
          NBSUM = NBFLUX
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
              YMPI(K,L) = XN(K,I)
              YMPI(K+3,L) = XNDOT(K,I)
              YMPI(K+6,L) = FI(K,I)
              YMPI(K+9,L) = FR(K,I)
              YMPI(K+12,L) = FIDOT(K,I)
              YMPI(K+15,L) = FRDOT(K,I)
              YMPI(K+18,L) = D1R(K,I)
              YMPI(K+21,L) = D2R(K,I)
              YMPI(K+24,L) = D3R(K,I)
*       Save corrections of force polynomials from routine fpcorr
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
          JMPI(11,j) = NBFLUX - NBSUM
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
      CALL MPI_BCAST(JMPI(1,rank+1),11,MPI_INTEGER,0,
     *                                  MPI_COMM_WORLD,ierr)
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      DO 812 L = 1,NREG
          I = IREG(L)
         DO 815 K = 1,3
              XN(K,I) = YMPI(K,L)
              XNDOT(K,I) = YMPI(K+3,L)
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
              NBFULL = NC5 + JMPI(1,J)
              NBVOID = NC6 + JMPI(2,J)
              NBSMIN = NC11 + JMPI(3,J)
              NBDIS2 = NC15 + JMPI(4,J)
              NRCONV = NC18 + JMPI(5,J)
              NBDIS = NC19 + JMPI(6,J)
              NLSMIN = NC20 + JMPI(7,J)
              NFAST = NC30 + JMPI(8,J)
              NICONV = NCX + JMPI(9,J)
              NBFAST = NCY + JMPI(10,J)
              NBFLUX = NBSUM + JMPI(11,J)
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
      NSTEPR = NSTEPR + NREG
      NBLCKR = NBLCKR + 1
      LSHRNK = .FALSE.
*
      DO 81 L = 1,NREG
          I = IREG(L)
*
              DO 816 K = 1,LMAX
 816          LIST(K,I) = IMPI(K,L)
*
*       Check minimum neighbor sphere since last output
              IF(LIST(1,I).GT.0)RSMIN = MIN(RSMIN,RS(I))
*
 81   CONTINUE
*
*          OPEN(98,STATUS='OLD',ERR=123)
*         print*,' last reg block t=',time,' length=',nreg
*         print*,' first 10 =',(name(ireg(l)),l=1,min(nreg,10))
*         call flush(6)
*         CLOSE(98)
*123      CONTINUE
*
          DO 54 L = 1,NREG
          I = IREG(L)
          NNB1 = LIST(1,I) + 1
          DO 56 K = 2,NNB1
          J = LIST(K,I)
          IPAIR = KVEC(J)
          IF(J.LT.IFIRST)THEN
          PRINT*,K,' 1 Warning I,J=',I,J,' IPAIR=',IPAIR,
     *       ' names,t=',NAME(I),NAME(J),TIME
          RIJ2=(X(1,I)-X(1,J))**2+(X(2,I)-X(2,J))**2+(X(3,I)-X(3,J))**2
          RSP2 = RS(N+IPAIR)**2
          RCRIT2 = 1.59*RS(I)**2
          PRINT*,' RIJ2/RS(N+IPAIR)=',RIJ2/RS(N+IPAIR),' RIJ2/RS2=',
     *    RIJ2/RCRIT2
          END IF
 56       CONTINUE
 54       CONTINUE

*     ETOTOL = ETOT
*     CALL ENERGY
*     ETOT = ZKIN-POT+EBIN
*     DETOTT = ETOT - ETOTOL
*     DETTTT = DETTTT + DETOTT
*     DETREL = DETOTT/MAX(ZKIN0,ETOT)
*     DETRTT = DETTTT/MAX(ZKIN0,ETOT)
*     WRITE(55,*)'---------------------------------------'
*     WRITE(55,*)' This was block NREG=',NREG
*     WRITE(55,*)' members ',(NAME(IREG(L)),L=1,NREG)
*     WRITE(55,*)' NREG=',NREG
*     WRITE(55,*)' t,ETOT,kin,pot,bin=',TIME,ZKIN-POT+EBIN,ZKIN,POT,EBIN
*     PRINT*,' NREG=',NREG,' members=',(NAME(IREG(L)),L=1,NREG)
*     PRINT*,' t,ETOT,kin,pot,bin=',TIME,ZKIN-POT+EBIN,ZKIN,POT
*     PRINT*,' t,DE,DECUM=',TIME,DETOTT,DETTTT
*     PRINT*,' t,DErel,DErelCUM=',TIME,DETREL,DETRTT
*     PRINT*,' ---------------------------------------'

      END IF
*
*       Copy all corrected coordinates & velocities (NB! only at the end).
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
   82     CONTINUE
   85 CONTINUE
*
*       Exit on KS termination, new multiple regularization or merger.
      IF (IQ.GT.0) THEN
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
      IF (KZ(18).GT.0.AND.LISTV(1).GT.0) THEN
          IF (DMOD(TIME,STEPM).EQ.0.0D0) THEN
              CALL SHRINK
              IF (LISTV(1).GT.0) THEN
                  CALL HIVEL(-1)
              END IF
          END IF
      END IF
*
*       Check optional mass loss time.
      IF (KZ(19).GT.0) THEN
*       Delay until time commensurate with 1000-year step (new polynomials).
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
          IF (NSUM.EQ.0.and.rank.eq.0) WRITE (6,90)
   90     FORMAT  (/,9X,'TERMINATION BY MANUAL INTERVENTION')
          CPU = 0.0
      END IF
      END IF
*
*       Repeat cycle until elapsed computing time exceeds the limit.
      CALL CPUTIM(TCOMP)
          TCOMP = (TCOMP-TTOTA)*60.
*
      IF (TCOMP.LT.CPU) GO TO 1
*
*       Do not terminate during triple, quad or chain regularization.
      IF (NSUB.GT.0) THEN
*       Specify zero step to enforce termination.
          DO 95 L = 1,NSUB
              STEPS(L) = 0.0D0
   95     CONTINUE
          NTIMER = NMAX
          GO TO 1
      END IF
*
*       Terminate run with optional COMMON save.
      IF (KZ(1).NE.0) THEN
          CPUTOT = CPUTOT + TCOMP - CPU0
          CALL MYDUMP(1,1)
          if(rank.eq.0)
     &    WRITE (6,98)  TOFF, TIME, TIME+TOFF, TCOMP, CPUTOT/60.0, 
     &                  ERRTOT, DETOT
   98     FORMAT (//,9X,'COMMON SAVED AT TOFF/TIME/TTOT =',3F16.8,
     &                  '  TCOMP =',F7.1,'  CPUTOT =',F6.1,
     &                  '  ERRTOT =',F10.6,'  DETOT =',F10.6)
      END IF
*
*        Determine time interval and step numbers per time unit
      TIMINT = TIME + TOFF - TINIT
*
#ifdef PARALLEL
          IF(rank.EQ.0)THEN
#endif
      WRITE (6,195)  rank,TIMINT,NSTEPI,NSTEPB,NSTEPR,NSTEPU
     &                          NSTEPU
      WRITE (6,195)  rank,TIMINT,NSTEPI-NIR,NSTEPB-NIB,NSTEPR-NRGL,
     &                          NSTEPU-NKS
  195 FORMAT (//,I9,' INTEGRATION INTERVAL =',F8.2,3X,' NIRR=',I11,
     &     ' NIRRB=',I11,' NREG=',I11,' NKS=',I11)
      WRITE (6,196)  (NSTEPI-NIR)/TIMINT,(NSTEPB-NIB)/TIMINT,
     &               (NSTEPR-NRGL)/TIMINT,(NSTEPU-NKS)/TIMINT
  196 FORMAT (//,9X,' PER TIME UNIT: NIRR=',1P,D12.5,' NIRRB=',
     &   D12.5,' NREG=',D12.5,' NKS=',D12.5)
#ifdef PARALLEL
          END IF
#endif
          STOP
*
 100      CONTINUE
*
*       Set current global time.
          TTOT = TIME + TOFF
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
          JJ = -1
          DO 88 JPAIR = 1,NPAIRS
          JJ = JJ + 2
          IF (LIST(1,JJ).GT.0) THEN
              ZZ = 1.0
              IZZ = -3
*       Distinguish between low and high-order prediction of U & UDOT.
              IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
              CALL KSRES2(JPAIR,J1,J2,ZZ,IZZ)
          END IF
   88     CONTINUE
      END IF
*
       RETURN
*
      END
