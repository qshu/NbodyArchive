      SUBROUTINE FPOLY0(RS0)
*
*
*       GPU initialization of neighbour lists and forces.
*       -------------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BARR/ ibarcount
      PARAMETER (NIMAX=1024)
      REAL*8   H2I(NIMAX),GPUACC(3,NIMAX),GPUJRK(3,NIMAX),GPUPHI(NIMAX),
     &         DTR(NIMAX)
      REAL*8   DV(3)
      INTEGER  LISTGP(LMAX,NIMAX)
      SAVE  LISTGP
#ifdef PARALLEL
      integer inum(NIMAX),ista(NIMAX)
#endif
*
*
*       Open the GPU libraries on new run (note nnbmax = NN is printed).
      NN = N + 10
      CALL GPUNB_OPEN(NN,rank)
*       Set larger value for GPUIRR (note further possible increase of NTOT).
*      NNN = N + 10
*      CALL GPUIRR_OPEN(NNN,LMAX)
*
*       Set provisional neighbour radius and regular time-step for GPUNB.
**!$omp parallel do private(I)
      DO 10 I = 1,N
*          IREG(I) = I
          T0(I) = 0.0
*       Modify neighbour radius according to NBLIST procedure.
          RI2 = X(1,I)**2 + X(2,I)**2 + X(3,I)**2
          RS(I) = RS0*SQRT(1.0 + RI2)
*       Set an estimated r-dependent regular time-step.
          STEPR(I) = SMAX/8.0D0*SQRT(1.0 + RI2)
          STEPR(I) = MIN(STEPR(I),SMAX)
          DO 5 K = 1,3
              F(K,I) = 0.0
              FDOT(K,I) = 0.0
    5     CONTINUE
   10 CONTINUE
**!$omp end parallel do
*
*       Send all particles (X0 & X0DOT) to the GPU for prediction.
**!$omp parallel do private(I)
*      DO 15 I = 1,N
*          CALL GPUIRR_SET_JP(I,X(1,I),XDOT(1,I),F(1,I),FDOT(1,I),
*     &                                          BODY(I),T0(I))
*   15 CONTINUE
**!$omp end parallel do
*
*
*       Send all single particles to the GPU.
      CALL GPUNB_SEND(N,BODY(IFIRST),X(1,IFIRST),XDOT(1,IFIRST))
*
*       Define maximum GPU neighbour number and initialize counters.
      NBMAX = MIN(NNBMAX + 150,LMAX-5)
      NOFL2 = 0
*
*       Loop over all particles split into NIMAX blocks.
#ifdef PARALLEL
      nl = n
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
      JNEXT = istart - 1      
#else
      JNEXT = 0
      istart = 1
      iend = n
#endif      
      DO 100 II = istart,iend,NIMAX
         NI = MIN(iend-JNEXT,NIMAX)
*       Copy neighbour radius, STEPR and state vector for each block.
**!$omp parallel do private(LL, I, K)
          DO 30 LL = 1,NI
              I = JNEXT + LL
              H2I(LL) = RS(I)**2
              DTR(LL) = STEPR(I)
   30     CONTINUE
**!$omp end parallel do
*
*       Evaluate forces, first derivatives and neighbour lists for new block.
 31       I = JNEXT + 1
          CALL GPUNB_REGF(NI,H2I,DTR,X(1,I),XDOT(1,I),GPUACC,GPUJRK,
     &                                     GPUPHI,LMAX,NNBMAX,LISTGP)
*       Check neighbour lists for overflow or zero membership (NNB = 1).
          DO 50 LL = 1,NI
              NNB = LISTGP(1,LL)
*       Repeat last block with reduced RS(I) on NNB < 2 (at end of loop).
              IF (NNB.LT.2) THEN
                  I = JNEXT + LL
                  RI2 = (X(1,I)-RDENS(1))**2 + (X(2,I)-RDENS(2))**2 +
     &                                         (X(3,I)-RDENS(3))**2
                  WRITE (41,40)  NSTEPR, NAME(I), NNB,
     &                           RS(I), SQRT(RI2)
   40             FORMAT (' OVERFLOW!   #R NAME NB RS ri ',
     &                                  I11,I7,I5,2F8.3)
                  CALL FLUSH(41)
                  IF (NNB.LT.0) THEN
                      RS(I) = 0.9*RS(I)
                  ELSE
                      RS(I) = 1.5*RS(I)
                  END IF
                  H2I(LL) = RS(I)**2
                  NOFL2 = NOFL2 + 1
              END IF
   50     CONTINUE
*
*       Repeat the last block for rare exceptions (NNB < 0 and = 1).
          IF (NOFL2.GT.0) THEN
              NOFL2 = 0
              GO TO 31
           END IF
*
*      Copy regular force and neighbour list from GPU.
**!$omp parallel do private(LL, I, ITEMP, NNB, L1, L)
          DO 70 LL = 1,NI
              I = JNEXT + LL
              DO 55 K = 1,3
                  FR(K,I) = GPUACC(K,LL)
                  D1R(K,I) = GPUJRK(K,LL)
   55         CONTINUE
              NNB = LISTGP(1,LL)
              L1 = 1
              DO 60 L = 2,NNB+1
*       Note GPU address starts from 0 (hence add IFIRST to neighbour list).
                  ITEMP = LISTGP(L,LL) + IFIRST
                  IF (ITEMP.NE.I) THEN
                      L1 = L1 + 1
                      LISTGP(L1,LL) = ITEMP
                  END IF
   60         CONTINUE
              LISTGP(1,LL) = L1 - 1
*              CALL GPUIRR_SET_LIST(I,LISTGP(1,LL))
              DO 65 L = 1,L1
                  LIST(L,I) = LISTGP(L,LL)
   65         CONTINUE
   70     CONTINUE
**!$omp end parallel do
*
*       Evaluate current irregular forces by vector procedure.
          DO 334 LL = 1,NI
             I = JNEXT + LL
*             PHII(I) = 0
             Do 3331 K = 1,3
                FI(K,I) = 0
                D1(K,I) = 0
 3331        CONTINUE
             NNB = LIST(1,I) + 1
             DO 335 L = 2,NNB
                K = LIST(L,I)
                A1 = X(1,K) - X(1,I)
                A2 = X(2,K) - X(2,I)
                A3 = X(3,K) - X(3,I)
                DV(1) = XDOT(1,K) - XDOT(1,I)
                DV(2) = XDOT(2,K) - XDOT(2,I)
                DV(3) = XDOT(3,K) - XDOT(3,I)
                RIJ2 = A1*A1 + A2*A2 + A3*A3
*     
                DR2I = 1.0/RIJ2
                DR3I = BODY(K)*DR2I*SQRT(DR2I)
                DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
                FI(1,I) = FI(1,I) + A1*DR3I
                FI(2,I) = FI(2,I) + A2*DR3I
                FI(3,I) = FI(3,I) + A3*DR3I
*                PHII(I) = PHII(I) - DR3I*RIJ2
                D1(1,I) = D1(1,I) + (DV(1) - A1*DRDV)*DR3I
                D1(2,I) = D1(2,I) + (DV(2) - A2*DRDV)*DR3I
                D1(3,I) = D1(3,I) + (DV(3) - A3*DRDV)*DR3I
 335         CONTINUE
 334      CONTINUE
*          CALL GPUIRR_FIRR_VEC(NI,IREG(II),GF(1,1),GFD(1,1))
*
*       Copy irregular force and first derivative.
c$$$*!$omp parallel do private(I)
c$$$          DO 80 LL = 1,NI
c$$$              I = JNEXT + LL
c$$$              DO 75 K = 1,3
c$$$                  FI(K,I) = GF(K,LL)
c$$$                  D1(K,I) = GFD(K,LL)
c$$$   75         CONTINUE
c$$$   80     CONTINUE
c$$$*!$omp end parallel do
*
          JNEXT = JNEXT + NI
  100 CONTINUE

*       Check option for external force.
      IF (KZ(14).GT.0) THEN
          CALL XTRNLD(istart,iend,1)
      END IF
*
*
#ifdef PARALLEL
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

      if(irank.eq.0)irank=isize
      istrec = ista(irank)
      icnt2 = inum(irank)

c$$$      print*, 'rank ',rank,' istsen ',istsen,' istrec-1 ',istrec-1,
c$$$     &     ' before fi(1) ',fi(1,istsen), ' fi(e) ',fi(1,istrec-1),
c$$$     &     ' size',icnt,' isend',isend,' irecv',irecv
*

#ifdef PUREMPI
      call cputim(tta)
      CALL MPI_SENDRECV(FI(1,istsen),3*icnt,MPI_REAL8,isend,
     *         rank,FI(1,istrec),3*icnt2,MPI_REAL8,irecv,irecv,
     *         MPI_COMM_WORLD,status,ierr)
      CALL MPI_SENDRECV(FR(1,istsen),3*icnt,MPI_REAL8,isend,
     *         rank,FR(1,istrec),3*icnt2,MPI_REAL8,irecv,irecv,
     *         MPI_COMM_WORLD,status,ierr)
      CALL MPI_SENDRECV(D1(1,istsen),3*icnt,MPI_REAL8,isend,
     *         rank,D1(1,istrec),3*icnt2,MPI_REAL8,irecv,irecv,
     *         MPI_COMM_WORLD,status,ierr)
      CALL MPI_SENDRECV(D1R(1,istsen),3*icnt,MPI_REAL8,isend,
     *         rank,D1R(1,istrec),3*icnt2,MPI_REAL8,irecv,irecv,
     *         MPI_COMM_WORLD,status,ierr)
      CALL MPI_SENDRECV(RS(istsen),icnt,MPI_REAL8,isend,
     *         rank,RS(istrec),icnt2,MPI_REAL8,irecv,irecv,
     *         MPI_COMM_WORLD,status,ierr)
      CALL MPI_SENDRECV(LIST(1,istsen),LMAX*icnt,MPI_INTEGER,isend,
     *         rank,LIST(1,istrec),LMAX*icnt2,MPI_REAL8,irecv,irecv,
     *         MPI_COMM_WORLD,status,ierr)
      call cputim(ttb)
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call cputim(tt999)
      ibarcount=ibarcount+1
      ttbar = ttbar + (tt999-ttb)*60
      xtsub2 = xtsub2 + dble((13*8+lmax*8)*(icnt+icnt2))
      ttsub2 = ttsub2 + (ttb-tta)*60.
c$$$      print*, 'rank ',rank,' istsen ',istsen,' istrec ',istrec,
c$$$     &     ' after fi(1) ',fi(1,istsen), ' fi(e) ',fi(1,istrec)
#endif
*
#ifdef SHMEM
       call barrier()
       call shmem_get(FI(1,istren),FI(1,istrec),3*icnt2,irecv)
       call shmem_get(FR(1,istren),FR(1,istrec),3*icnt2,irecv)
       call shmem_get(D1(1,istren),D1(1,istrec),3*icnt2,irecv)
       call shmem_get(D1R(1,istren),D1R(1,istrec),3*icnt2,irecv)
       call shmem_get(RS(1,istren),RS(1,istrec),icnt2,irecv)
       call shmem_get(LIST(1,istren),LIST(1,istrec),lmax*icnt2,irecv)
#endif
*
 1002 continue
#endif
      
*     Form total force & force derivative and extra variables for XVPRED.
**!$omp parallel do private(I)
      DO 110 I = 1,N
          DO 105 K = 1,3
              F(K,I) = FI(K,I) + FR(K,I)
              FDOT(K,I) = D1(K,I) + D1R(K,I)
              D0(K,I) = FI(K,I)
              D0R(K,I) = FR(K,I)
              FIDOT(K,I) = D1(K,I)
              FRDOT(K,I) = D1R(K,I)
  105     CONTINUE
  110 CONTINUE
**!$omp end parallel do
*

*       Close the GPU libraries (limits change in INTGRT).
      CALL GPUNB_CLOSE
*      CALL GPUIRR_CLOSE
*
      RETURN
*
      END
