      SUBROUTINE INTGRT
*
*
*       N-body integrator flow control.
*       -------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'timing.h'
*      include 'omp_lib.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/BLKLVL/JHIST(0:NMAX),JHISTR(0:NMAX),JHISTU(0:NMAX)
      COMMON/STSTAT/  TINIT,NIR,NIB,NRGL,NKS
      INTEGER JHIST,JHISTR,LSHORT(NMAX),TBLST,TNXLST(NMAX),NTCOUNT
      INTEGER NTLSTC
      INTEGER NXTLST(NMAX),IREG(NMAX)
      REAL*8 TREGSIT
      LOGICAL LSTEPM,NTLOOP,IPSTART
#ifdef PARALLEL
      integer inum(maxpe),ista(maxpe),NXTLST_MPI(NMAX+500),NPIECE
      integer istart_lst,iend_lst
      REAL*8 XMPI(19,NMAX),YMPI(41,NMAX)
      SAVE NPIECE
#endif
      INTEGER IMPI(LMAX,NMAX),JMPI(12,maxpe), JMPILOCAL(12,maxpe)
      DATA LSHORT(1),NDTK,NTLOOP,NTLSTC,NTCOUNT/0,11,.false.,0,0/
      DATA IQ,ICALL,LSTEPM,STEPM,IPSTART /0,2,.FALSE.,0.03125,.true./
      SAVE IQ,ICALL,LSTEPM,STEPM,IPSTART,TBLST,TREGSIT,TNXLST
      SAVE NDTK,NTLOOP,NTLSTC,NTCOUNT
*     --09/26/13 22:52-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      INTEGER NGPUC
c$$$      SAVE NGPUC
c$$$      DATA NGPUC/0/
c$$$      REAL*8 TCCT
c$$$      DATA TCCT/0.0/
c$$$      SAVE TCCT
*     --09/26/13 22:52-lwang-end----------------------------------------*

*
*
*       Update quantized value of STEPM for large N (first time only).
      call cputim(tttbegin)
      IF (.NOT.LSTEPM.AND.NZERO.GT.1024) THEN
         K = (FLOAT(NZERO)/1024.0)**0.333333
         STEPM = 0.03125D0/2**(K-1)
         LSTEPM = .TRUE.
      END IF
*

*       Open GPU for all single particles
      IF (IPSTART) THEN
#ifdef GPU
         NN = NTOT + 10
         CALL GPUNB_OPEN(NN,rank)
         CALL GPUPOT_INIT_FLOAT(rank)
         CALL GPUPOT_INIT(rank,NN)
#endif
#ifdef PARALLEL         
*       Get the offset number for NXTLST_MPI
         IF (MOD(NMAX,ISIZE).EQ.0) THEN
            NPIECE = NMAX / ISIZE
         ELSE
            NPIECE = NMAX / ISIZE + 1
         END IF
#endif         
         IPSTART = .false.
      END IF
      call cputim(tttiend)
      ttintb = ttintb +(tttiend-tttbegin)*60
*
*       Search for high velocities after escape or KS/chain termination.
  999 IF (KZ(37).GT.0.AND.(IPHASE.EQ.-1.OR.IPHASE.GE.2)) THEN
          CALL HIVEL(0)
      END IF
*
      call cputim(tttbegin)
*       Reset control & regularization indicators.
      IPHASE = 0
      IKS = 0
*       Initialize end-point of integration times and set TMIN.
      TMIN = 1.0D+10
*       Initialize end-point of integration times and set TMIN & IMIN.
      DO I = IFIRST,NTOT
         TIMENW(I) = T0(I) + STEP(I)
         IF(TIMENW(I).LT.TMIN)THEN
            TMIN = TIMENW(I)
            IMIN = I
         END IF
      END DO
      TREGSIT = TMIN - DTK(NDTK-1)
*      if(rank.eq.0) print*,'Tmin ',tmin
*
      IF (IQ.LT.0) ICALL = 0
      IQ = 0
      call cputim(tttiend)
      ttintb = ttintb +(tttiend-tttbegin)*60
*
*       Find all particles due at next block time.
    1 CONTINUE
*
      NXTLEN = 0
*
*       Redetermine TMIN after main change to catch new small steps after chain.
      ICALL = ICALL + 1
      IF (ICALL.EQ.2) GO TO 999
*
*     determine next block particles without assuming
*     sorted time step list (R.Sp.)
*
      call cputim(tttbegin)

*     Check regular force condition for small block memberships.
*     Use TNXLST array to store the next small time step list
      NREG = 0
      If(TMIN-TREGSIT.GE.DTK(NDTK)) THEN
         IF (NTLOOP) THEN
            IF (NTCOUNT.LE.1) THEN
               NDTK = NDTK - 3
            ELSEIF (NTCOUNT.LE.7) THEN
               NDTK = NDTK - 1
            END IF
            IF (NTLSTC.GE.NTOT/8) THEN
               NDTK = NDTK + 3
            ELSEIF (NTLSTC.GE.NTOT/128) THEN
               NDTK = NDTK + 1
            END IF
            NTCOUNT = 0
            NTLSTC = 0
         ELSE
            NTLOOP=.true.
         END IF
*     --08/12/13 13:04-lwang-debug--------------------------------------*
***** Note: print whether it is full prediction------------------------**
*         if(rank.eq.0) then
*            print*,'T-RENEW TMIN TBLST NDTK',
*     &           TMIN,TBLST,NDTK
*         END if
*     --08/12/13 13:04-lwang-end-debug----------------------------------*
         TBLST=0
         DO J = IFIRST, NTOT
            IF(TIMENW(J)-TMIN.LT.DTK(NDTK)) THEN
               TBLST = TBLST + 1
               TNXLST(TBLST)=J
               IF(DABS(TIMENW(J)-TMIN).LT.DTK(40)) THEN
#ifdef PARALLEL                  
                  NXTLST_MPI(MOD(NXTLEN,isize)*NPIECE + NXTLEN/ISIZE +1)
     &                 = J
#endif                  
                  NXTLEN = NXTLEN + 1
                  NXTLST(NXTLEN) = J
                  IF (TIMENW(J).GE.T0R(J) + STEPR(J)) THEN
                     NREG = NREG + 1
                     IREG(NREG) = J
                  END IF
               END IF
            END IF
         END DO
         TREGSIT = TMIN
         NTLSTC = NTLSTC + TBLST
      ELSE
*     If not reach the maximum stored time step, search array
*     --08/12/13 13:03-lwang-debug--------------------------------------*
***** Note: print whether it work--------------------------------------**
*         if(rank.eq.0) then
*            print*,'T-Sublist TMIN TREGSIT TBLST',TMIN,TREGSIT,TBLST
*         end if
*     --08/12/13 13:03-lwang-end-debug----------------------------------*
         NTCOUNT = NTCOUNT + 1
         DO L = 1, TBLST
            J = TNXLST(L)
            IF(DABS(TIMENW(J)-TMIN).LT.DTK(40)) THEN
#ifdef PARALLEL                  
               NXTLST_MPI(MOD(NXTLEN,isize)*NPIECE + NXTLEN/isize + 1)
     *              = J
#endif               
               NXTLEN = NXTLEN + 1
               NXTLST(NXTLEN) = J
               IF (TIMENW(J).GE.T0R(J) + STEPR(J)) THEN
                  NREG = NREG + 1
                  IREG(NREG) = J
               END IF
            END IF
         END DO
      END IF

*
      LSKIP = NXTLEN 
      IF (LSKIP.LE.50) THEN
*       Update short timestep list for regularization search.
          CALL SHORT(NXTLEN,NXTLST,LSHORT)
      END IF
*
*     --11/08/13 23:15-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(rank.eq.0.and.time.ge.TCCT) then
c$$$         write(120,*) time,ifirst,n,ntot
c$$$         do L=1,NTOT
c$$$            write(120,*) name(L),tev(L)-tev0(L),step(L),stepr(L)
c$$$         END DO
c$$$         TCCT = TCCT + 0.1
c$$$      end if
*     --11/08/13 23:15-lwang-end----------------------------------------*
*       Set new time and save block time (for regularization terminations).
      TIME = TMIN
      TBLOCK = TIME
*     --11/09/13 19:12-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(rank.eq.0) print*,'NXT ',time,nxtlen,nxtlst(1:nxtlen)
c$$$      call flush(6)
*     --11/09/13 19:12-lwang-end----------------------------------------*
*     --11/03/13 10:22-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
*      if(time+toff.ge.968.46875) then
*         if(rank.eq.0) call mydump(1,1)
*         call mpi_barrier(MPI_COMM_WORLD,ierr)                     
*         stop
*      end if
*     --11/03/13 10:22-lwang-end----------------------------------------*
*
*     Check option for advancing interstellar clouds.
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
      call cputim(tttiend)
      ttintb = ttintb +(tttiend-tttbegin)*60

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
          TIME = TADJ
          IPHASE = 3
          GO TO 100
      END IF
*
*       Also check output time in case DTADJ & DELTAT not commensurate.
      IF (TIME.GT.TNEXT) THEN
          TIME = TNEXT
          CALL OUTPUT
          GO TO 1
      END IF
*
*       See whether to advance any close encounters at first new time.
      IF (TIME.GT.TPREV) THEN
         call cputim(tt5)
*     --08/27/13 16:25-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$            if (time.ge.3.90625000000000000E-003) then
c$$$               do iii=1,ntot
*         if(name(iii).eq.5287) then
*         write(113+rank,*)'5287a',iii,
*     *           t0(iii),step(iii),ibarcount,time,tprev
*           call flush(6)
*         end if
c$$$                  if(name(iii).eq.23871) then
*     print*,rank,'suba',iii,x0(1,iii),x0dot(1,iii),time,tprev
c$$$                    print*,rank,'suba',iii,t0(iii),step(iii),time,tprev,
c$$$     *                    x(1,iii),x0(1,iii)
c$$$                    call flush(6)
c$$$                  end if
c$$$               end do
c$$$            end if
*          if(time.ge.6E-002) print*,rank,time
*          if(time.ge.6.6E-02) stop
*     --09/04/13 15:10-lwang-end----------------------------------------*
c$$$      call flush(6)
c$$$      call mpi_barrier(MPI_COMM_WORLD,ierr)            
c$$$      if(time.ge.TCCT) print*,rank,'asub ',iq,nsub,time,tprev
*     --08/27/13 16:25-lwang-end-debug----------------------------------*
         CALL SUBINT(IQ,I10)
*     --08/27/13 16:31-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(time.ge.TCCT) print*,rank,'bsub ',iq,nsub,time,tprev
c$$$            if (time.ge.3.90625000000000000E-003) then
c$$$               do iii=1,ntot
c$$$                  if(name(iii).eq.23871) then
c$$$c$$$*         write(113+rank,*)'5287b',iii,
c$$$c$$$*     *           t0(iii),step(iii),ibarcount,time,tprev
c$$$c$$$*         end if
c$$$c$$$         if(time.ge.1.41601562500000000E-002.and.name(iii).eq.5287) then
c$$$                    print*,rank,'subb',iii,t0(iii),step(iii),time,tprev,
c$$$     *                    x(1,iii),x0(1,iii)
c$$$                    call flush(6)
c$$$c$$$           stop
c$$$                 end if
c$$$              end do
c$$$              end if
*     --08/27/13 16:31-lwang-end-debug----------------------------------*
         call cputim(tt6)
         ttks = ttks + (tt6-tt5)*60.
*
         IF (IQ.LT.0) GO TO 999
      END IF

*
*       Check regularization criterion for single particles.
      call cputim(tttbegin)
      IKS = 0
      ISMIN = 0
      DSMIN = DTMIN
      IF (LSKIP.LE.50) THEN
*       Search only in prepared list of short-step particles. (R.Sp.)
         ISHORT = LSHORT(1)
         DO L = 2,ISHORT+1
            I = LSHORT(L)
*       Search for minimum timestep candidate for not ordered steps (R.Sp.)
*       Beware that members of LSHORT may be members of KS pair (R.Sp.)
c$$$          If(step(i).LT.1.0D-8) THEN
c$$$             WRITE (6,24)  TIME, I, NAME(I), NXTLEN, NSTEPI,
c$$$     &            STEP(I), STEPR(I), NTOT, NPAIRS
c$$$ 24          FORMAT (' SMALL STEP!!  T I NAME LEN #I SI SR NT NP',F10.5,
c$$$     &            3I6,I11,1P,2E10.2,2I11)
c$$$             CALL FLUSH(6)
c$$$             STOP
c$$$  END IF
            
            IF (STEP(I).LT.DTMIN.AND.STEP(I).LT.DSMIN.AND.I.LE.N.AND.
     &           I.GE.IFIRST) THEN
               DSMIN = STEP(I)
               ISMIN = I
            END IF
         END DO
      END IF
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
*     Check eccentricity (cf. max perturbation) and neighbour radius.
                  IF (ECC.LT.0.25.AND.APO.LT.0.02*RS(ISMIN)) THEN
                     PRINT*, ' KS TRY: NAM E A EB ',
     *                    NAME(ISMIN), NAME(JMIN), ECC, SEMI, EB
                     CALL FLUSH(6)
                     IKS = IKS + 1
                     ICOMP = ISMIN
                     JCOMP = JMIN
                  END IF
               END IF
            END IF
         END IF
      END IF

      call cputim(tt1)
      ttintb=ttintb +(tt1-tttbegin)*60.
*
*       Initialize counters for irregular & regular integrations.
      TPREV = TIME
      NBLOCK = NBLOCK + 1
      TMIN = 1.0D+10

*     Initialize counters for irregular & regular integrations.
*      NREG = 0
*
*     write(6,*)' irr ',time,nxtlen,(nxtlst(k),step(k),k=1,5)
*     call flush(6)
*       Advance the irregular step for all particles in the current block.
*       Block-Step Level Diagnostics (R.Sp. 29.Apr. 1993)
      IF(KZ(33).GT.0)JHIST(NXTLEN) = JHIST(NXTLEN) + 1
*
*     --08/27/13 16:38-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
*****      GET NEIGHBOR LIST      
c$$$      if(rank.eq.0.and.nxtlen.ge.5000) then
c$$$         do j=1,nxtlen
c$$$            jj = nxtlst(j)
c$$$            write(102+NGPUC,*) jj, list(1,jj)
c$$$         end do
c$$$*         tcct = tcct + 0.1
c$$$         NGPUC = NGPUC + 1
c$$$      end if
c$$$      if(time.ge.3.5373535156250000) then
c$$$         do j=ifirst,ntot
c$$$      j=105112
c$$$            if(name(j).eq.3173) then
c$$$             write(102+rank,*),j,name(j),'x0',(x0(kk,j),kk=1,3),
c$$$     *           'x0dot',(x0dot(kk,j),kk=1,3),'t0',t0(j),'step',step(j),
c$$$     *           'stepr',stepr(j),
c$$$     *           'f',(f(kk,j),kk=1,3),'fdot',(fdot(kk,j),kk=1,3),
c$$$     *           'd0',(d0(kk,j),kk=1,3),'d1',(d1(kk,j),kk=1,3),
c$$$     *           'd2',(d2(kk,j),kk=1,3),'d3',(d3(kk,j),kk=1,3),
c$$$     *           'd0r',(d0r(kk,j),kk=1,3),'d1r',(d1r(kk,j),kk=1,3),
c$$$     *           'd2r',(d2r(kk,j),kk=1,3),'d3r',(d3r(kk,j),kk=1,3),
c$$$     *           'body',body(j),'time',time,'list',list(1,j)
c$$$            end if
c$$$         end do
c$$$      end if
c$$$      call mpi_barrier(MPI_COMM_WORLD,ierr)
c$$$      if(time.ge.TCCT) then
c$$$         do i=1,ifirst-1
c$$$            if(name(i).eq.16107) then
c$$$            write(100+rank,*) i,name(i),step(i),t0(i),time,body(i),
c$$$     *           gamma((i+1)/2),list(1,i),name((i+1)/2+N)
c$$$c$$$*               write(100+rank,*) i,name(i),x0(1,i),x0dot(1,i),
c$$$c$$$*     *              f(1,i),fdot(1,i),time
c$$$c$$$            end if
c$$$c$$$         end do
c$$$         do j=ifirst,ntot
c$$$            if(name(j).eq.16107) then
c$$$            write(100+rank,*),j,name(j),'x0',(x0(kk,j),kk=1,3),
c$$$     *        'x0dot',(x0dot(kk,j),kk=1,3),'t0',t0(j),
c$$$     &           'step',step(j),'stepr',stepr(j),
c$$$     *        'f',(f(kk,j),kk=1,3),'fdot',(fdot(kk,j),kk=1,3),
c$$$     *       'fi',(fi(kk,j),kk=1,3),'fidot',(fidot(kk,j),kk=1,3),
c$$$     *        'd0',(d0(kk,j),kk=1,3),'d1',(d1(kk,j),kk=1,3),
c$$$     *        'd2',(d2(kk,j),kk=1,3),'d3',(d3(kk,j),kk=1,3),
c$$$     *        'd0r',(d0r(kk,j),kk=1,3),'d1r',(d1r(kk,j),kk=1,3),
c$$$     *        'd2r',(d2r(kk,j),kk=1,3),'d3r',(d3r(kk,j),kk=1,3),
c$$$     *        'body',body(j),'time',time,'list',list(1,j)
c$$$            call flush(100+rank)
c$$$            end if
c$$$         end do
c$$$         call mpi_barrier(MPI_COMM_WORLD,ierr)
c$$$*        if (TCCT.gt.0.111) stop
c$$$*        TCCT = TCCT + 0.001
c$$$      end if
c$$$      if(time.ge.0.10729980468750000) stop
c$$$      if(time.ge.12.920412063598633) stop
c$$$      if (time.ge.TCCT) then
c$$$         do i=1,nreg
c$$$            j=ireg(i)
c$$$            write(105+rank,*),j,name(j),
c$$$            'x0',(x0(kk,j),kk=1,3),
c$$$     *        'x0dot',(x0dot(kk,j),kk=1,3),'t0',t0(j),'step',step(j),
c$$$     *        'f',(f(kk,j),kk=1,3),'fdot',(fdot(kk,j),kk=1,3),
c$$$     *        'd0',(d0(kk,j),kk=1,3),'d1',(d1(kk,j),kk=1,3),
c$$$     *        'd2',(d2(kk,j),kk=1,3),'d3',(d3(kk,j),kk=1,3),
c$$$     *        'd0r',(d0r(kk,j),kk=1,3),'d1r',(d1r(kk,j),kk=1,3),
c$$$     *        'd2r',(d2r(kk,j),kk=1,3),'d3r',(d3r(kk,j),kk=1,3),
c$$$     *        'body',body(j),'time',time,'list',list(1,j)
c$$$         end do
c$$$      end if
c$$$c$$$      if (time.ge.0.02) stop
c$$$c$$$      if(time.ge.0.31) then
c$$$      if(time.ge.60.37) then
c$$$        print*, 'r NXTLEN :',rank,nxtlen,time,ibarcount,nreg,ireg(1)
c$$$        call flush(6)
c$$$        call mpi_barrier(MPI_COMM_WORLD,ierr)        
c$$$      end if
c$$$         if (time.ge.20.852127075195313) then
c$$$            print*,rank,time,'nxtlst ',(nxtlst(k),k=1,nxtlen)
c$$$         end if
c$$$      if(time.ge.54.192613631486893) stop
*     --08/27/13 16:38-lwang-end-debug----------------------------------*
#ifdef PARALLEL
      IF(NXTLEN.LE.isernb)THEN
#endif
*
*     Predict x and xdot
         call cputim(ttnb1)
         call xbpred(NREG,1,NXTLEN,NXTLST)
*       Resolve Chain if it is in block or in neighbour lists (RS Nov. 03).
         IF (NCH.GT.0) THEN
            ICHPR = 0
*     First check whether chain c.m. is in current block (ICHPR = 2).
            IF (TIMENW(ICH).EQ.TIME) ICHPR = 2
*     Second check whether neighbour lists contain chain c.m. (ICHPR = 1).
            IF (ICHPR.EQ.0) THEN
               DO 48 L = 1,NXTLEN
                  I = NXTLST(L)
                  NNB1 = LIST(1,I) + 1
*     Second check whether neighbour lists contain chain c.m. (ICHPR = 1).
                  DO K = 2,NNB1
                     J = LIST(K,I)
                     IF (J.GT.ICH) GO TO 48
                     IF (J.EQ.ICH) ICHPR = 1
                  END DO
 48            CONTINUE
            END IF
            IF (ICHPR.GT.1) CALL CHLIST(ICH)
            IF (ICHPR.GT.0) CALL XCPRED(0)
         END IF
*
         call cputim(ttnb2)
         ttpre = ttpre + (ttnb2-ttnb1)*60.

*     Neighbor Integration
*!$omp parallel do if(NXTLEN.GE.ITHREAD) private(L,I)
         DO L = 1,NXTLEN
            I = NXTLST(L)
            CALL NBINT(I)
            DO K = 1,3
               X0(K,I) = XN(K,I)
               X0DOT(K,I) = XNDOT(K,I)
            END DO
         END DO
*!$omp end parallel do

         call cputim(tt3)
         ttnbi = ttnbi +(tt3-ttnb2)*60
         ttirr = ttirr + (tt3-ttnb1)*60.
*
#ifdef PARALLEL
*       start PARALLEL section for nbint
        ELSE
*
           call cputim(tt998)
           call mpi_barrier(MPI_COMM_WORLD,ierr)
           call cputim(tt999)
           ttbar = ttbar + (tt999-tt998)*60
           ibarcount=ibarcount+1
           ttbarnb = ttbarnb + (tt999-tt998)*60
*           print*, 'Nbint Barrier 1: rank iphase ttbarnb ttbar',
*     &          ' dt nxtlen bcount',rank,ixxx,ttbarnb,ttbar,tttdet,
*     &          nxtlen,ibarcount
*
           istart_lst = rank * NPIECE + 1
           ishift = MOD(NXTLEN,isize)
           idivide = NXTLEN / isize
           if (rank.LT.ishift) then
              iend_lst = idivide + istart_lst
           else
              iend_lst = idivide + istart_lst - 1
           end if
           
           nl = NXTLEN
*     
           inl = nl/isize
           jsize = isize*inl
           idiff = nl - jsize
           irun = 0
*     
           do ix = 1,isize
              inum(ix)=inl
              if(ix.le.idiff)inum(ix) = inum(ix) + 1
              ista(ix) = irun+1
              if(ista(ix).gt.nl)inum(ix) = 0
              irun = irun + inum(ix)
           end do
*     
           istart = ista(rank+1)
           iend = ista(rank+1) + inum(rank+1) - 1
*
*     if(time.lt.0.1d0)then
*     print*,' rank ',rank,' NXTLEN ',NXTLEN,TIME,
*    *    ' istart,iend=',istart,iend
*     end if
*
c$$$      call cputim(tt998)
c$$$      call mpi_barrier(MPI_COMM_WORLD,ierr)
c$$$      call cputim(tt999)
c$$$      tttdet=(tt999-tt998)*60
c$$$      ttbar = ttbar + (tt999-tt998)*60
c$$$      ibarcount=ibarcount+1
c$$$      ttbarnb = ttbarnb + (tt999-tt998)*60
*      print*, 'Nbint Barrier 2: rank iphase ttbarnb ttbar',
*     &     ' dt nxtlen bcount',rank,ixxx,ttbarnb,ttbar,tttdet,
*     &     nxtlen,ibarcount
*
*     if(time.lt.0.1d0)then
*     print*,' rank ',rank,' NXT in DO L ',L,TIME,
*    *    ' istart,iend=',istart,iend
*     end if

*     Predict x and xdot
           call cputim(ttnb1)
*     --11/09/13 19:09-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$           print*,'TEST ',time,rank,istart_lst,iend_lst,
c$$$     &          nxtlst_mpi(istart_lst:iend_lst)
c$$$           call flush(6)
c$$$           call mpi_barrier(MPI_COMM_WORLD,ierr)
*     --11/09/13 19:09-lwang-end----------------------------------------*
           call xbpred(NREG,istart_lst,iend_lst,NXTLST_MPI)
*       Resolve Chain if it is in block or in neighbour lists (RS Nov. 03).
           IF (NCH.GT.0) THEN
              ICHPR = 0
*       First check whether chain c.m. is in current block (ICHPR = 2).
              IF (TIMENW(ICH).EQ.TIME) ICHPR = 2
*       Second check whether neighbour lists contain chain c.m. (ICHPR = 1).
              IF (ICHPR.EQ.0) THEN
                 DO 47 L = 1,NXTLEN
                    I = NXTLST(L)
                    NNB1 = LIST(1,I) + 1
*       Second check whether neighbour lists contain chain c.m. (ICHPR = 1).
                    DO K = 2,NNB1
                       J = LIST(K,I)
                       IF (J.GT.ICH) GO TO 47
                       IF (J.EQ.ICH) ICHPR = 1
                    END DO
 47              CONTINUE
              END IF
              IF (ICHPR.GT.1) CALL CHLIST(ICH)
              IF (ICHPR.GT.0) CALL XCPRED(0)
           END IF
*
           call cputim(ttnb2)
           ttpre = ttpre + (ttnb2-ttnb1)*60.

*     Neighbor Integration
*!$omp parallel do if(iend-istart.GE.ITHREAD) private(L,I)
           LL = istart
           DO L = istart_lst,iend_lst
              I = NXTLST_MPI(L)
              CALL NBINT(I)
              XMPI(1:3,LL) = XN(1:3,I)
              XMPI(4:6,LL) = XNDOT(1:3,I)
              XMPI(7:9,LL) = D2(1:3,I)
              XMPI(10:12,LL) = D3(1:3,I)
              XMPI(13:15,LL) = FI(1:3,I)
              XMPI(16:18,LL) = FIDOT(1:3,I)
              XMPI(19,LL) = STEP(I)
              LL = LL + 1
           END DO
*!$omp end parallel do       
*
           call cputim(tt3)
           ttnbi = ttnbi + (tt3-ttnb2)*60.
           ttirr = ttirr + (tt3-ttnb1)*60.
*
*        Distribute variables into private vectors again T3E (R.Sp.)
*
           isend = rank + 1
           if(isend.eq.isize)isend = 0
           irecv = rank - 1
           if(irecv.eq.-1)irecv = isize - 1
*
           do ir = 0,isize-2
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
              call cputim(tta)
              CALL MPI_SENDRECV(XMPI(1,istsen),19*icnt,MPI_REAL8,isend,
     *             rank,XMPI(1,istrec),19*icnt2,MPI_REAL8,irecv,irecv,
     *             MPI_COMM_WORLD,status,ierr)
              call cputim(ttb)
              call mpi_barrier(MPI_COMM_WORLD,ierr)
              call cputim(tt999)
              ttbar = ttbar + (tt999-ttb)*60
              ttbarnb = ttbarnb + (tt999-ttb)*60
              ibarcount=ibarcount+1
*     print*, 'Nbint Barrier 3: rank iphase ttbarnb ttbar',
*     &     ' dt nxtlen bcount',rank,ixxx,ttbarnb,ttbar,tttdet,
*     &     nxtlen,ibarcount
              xtsub1 = xtsub1 + dble(19*8*(icnt+icnt2))
              ttsub = ttsub + (ttb-tta)*60.
#endif
*
#ifdef SHMEM
              call barrier()
              call shmem_get(XMPI(1,istrec),XMPI(1,istrec),19*icnt2,
     &             irecv)
#endif
*
           end do
*
*          call mpi_barrier(MPI_COMM_WORLD,ierr)
*
           call cputim(ttmov3)
**!$omp parallel do if(NXTLEN.GE.ithread)
**!$omp& default(shared) private(L,I,K,LL)
           LL = 1
           LLrank = 1
           IF (ishift.ne.0) then
              LLoffset = idivide +1
           else
              LLoffset = idivide
           end if
           DO L = 1,NXTLEN
              I = NXTLST_MPI(LL)
*     --11/09/13 19:29-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
*              if(rank.eq.0) print*,'XMPI',L,LL,I,XMPI(1,L)
*     --11/09/13 19:29-lwang-end----------------------------------------*
              X0(1:3,I) = XMPI(1:3,L)
              X0DOT(1:3,I) = XMPI(4:6,L)
              D2(1:3,I) = XMPI(7:9,L)
              D3(1:3,I) = XMPI(10:12,L)
              FI(1:3,I) = XMPI(13:15,L)
              FIDOT(1:3,I) = XMPI(16:18,L)
              STEP(I) = XMPI(19,L)
              LL = LL + 1
              if (LL.GT.LLoffset) then
                 LL = LLrank*NPIECE + 1
                 if (LLrank.ge.ishift.and.(LLrank-1).lt.ishift) then
                    LLoffset = LLoffset + NPIECE - 1
                 else
                    LLoffset = LLoffset + NPIECE
                 end if
                 LLrank = LLrank + 1
              end if
           END DO
**!$omp end parallel do
*     --11/09/13 19:30-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
*           stop
*     --11/09/13 19:30-lwang-end----------------------------------------*
*
           call cputim(tt4)
           ttmov = ttmov + (tt4-ttmov3)*60.
*     
        END IF
*          End PARALLEL section for nbint
#endif
*

        call cputim(ttt32)
**!$omp parallel do if(NXTLEN.GE.ithread)
**!$omp& default(shared) private(L,I,DTR,K)
        DO L = 1,NXTLEN
           I = NXTLST(L)
*       Save new block step and update T0 & next time
           T0(I) = TIME
           TIMENW(I) = T0(I) + STEP(I)
           
*     Set non-zero indicator for new regular force.
           IF (T0R(I) + STEPR(I).GT.TIME) THEN
*     *!$omp critical 
*              NREG = NREG + 1
*              IREG(NREG) = I
**!$omp end critical
*          ELSE
*       Extrapolate regular force & first derivatives to obtain F & FDOT.
              DTR = TIME - T0R(I)
              F(1,I) = 0.5*(FRDOT(1,I)*DTR + FR(1,I) + FI(1,I))
              F(2,I) = 0.5*(FRDOT(2,I)*DTR + FR(2,I) + FI(2,I))
              F(3,I) = 0.5*(FRDOT(3,I)*DTR + FR(3,I) + FI(3,I))
              FDOT(1,I) = ONE6*(FRDOT(1,I) + FIDOT(1,I))
              FDOT(2,I) = ONE6*(FRDOT(2,I) + FIDOT(2,I))
              FDOT(3,I) = ONE6*(FRDOT(3,I) + FIDOT(3,I))
* Higher order extrapolation?
*                 F(K,I) = FI(K,I) + FR(K,I) + DTR*(FRDOT(K,I)
*    *                + DTR*(D2R(K,I)/2.D0 + DTR*D3R(K,I)/6.D0))
*                 FDOT(K,I) = FIDOT(K,I) + FRDOT(K,I)
*    *                + DTR*(D2R(K,I) + DTR*D3R(K,I)/2.D0)
*                 F(K,I) = F(K,I)/2.D0
*                 FDOT(K,I) = FDOT(K,I)/6.D0
           END IF
          
*     X0(K,I) = XN(K,I)
*     X0DOT(K,I) = XNDOT(K,I)
           D0(1:3,I) = FI(1:3,I)
           D1(1:3,I) = FIDOT(1:3,I)
        END DO
**!$omp end parallel do
*
*     --09/26/13 12:46-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$       if(time.ge.12.920379638671875) then
c$$$         do i=1,ifirst
c$$$c$$$            if(name(i).eq.3233) then
c$$$            write(100+rank,*) i,name(i),step(i),t0(i),time
c$$$c$$$*               write(100+rank,*) i,name(i),x0(1,i),x0dot(1,i),
c$$$c$$$*     *              f(1,i),fdot(1,i),time
c$$$c$$$            end if
c$$$         end do
c$$$         do i=ifirst,ntot
c$$$c$$$            if(name(i).eq.3233) then
c$$$c$$$c$$$            if(name(i).eq.467.or.name(i).eq.6599.
c$$$c$$$c$$$c$$$     &           or.name(i).eq.6662.or.name(i).eq.35667) then
c$$$           write(100+rank,*),i,name(i),x0(1,i),x0dot(1,i),
c$$$     *      step(i),t0(i),f(1,i),fdot(1,i),time,'nb'
c$$$c$$$            
c$$$c$$$c$$$     &              ,time,ibarcount
c$$$c$$$*               call flush(6)
c$$$c$$$            end if
c$$$         end do
c$$$      call mpi_barrier(MPI_COMM_WORLD,ierr)
c$$$         stop
c$$$      end if
*     --09/26/13 12:46-lwang-end----------------------------------------*
*           CALL nemo_savestate(n,3,time,body,x,xdot)
C$$$       if (time.ge.0.41) then
C$$$          print*,rank,'4821a',x(1,4821),xdot(1,4821),f(1,4821)
C$$$     &        ,fdot(1,4821),time,timenw(4821),t0(4821),step(4821)
C$$$          print*,rank,'2411+Na',x(1,2411+N),xdot(1,2411+N),
C$$$     &   f(1,2411+N),fdot(1,2411+N),time,timenw(2411+N),t0(2411+N),
C$$$     &     step(2411+N)
C$$$      end if

*
*       See whether any KS candidates are in the same block.
        IF (IKS.GT.0) THEN
*       Accept same time, otherwise reduce STEP(ICOMP) and/or delay.
           IF (T0(JCOMP).EQ.T0(ICOMP)) THEN
*             call cputim(ttnb1)
*             call xbpredall
*             call cputim(ttnb2)
*             ttpre = ttpre + (ttnb2-ttnb1)*60.
              I = ICOMP
              ICOMP = MIN(ICOMP,JCOMP)
              JCOMP = MAX(I,JCOMP)
*     --09/25/13 21:44-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$              if(time.ge.0.3) then
c$$$                 print*,rank,'icomp',icomp,name(icomp),
c$$$     &            x(1,icomp),xdot(1,icomp),x0(1,icomp),x0dot(i,icomp)
c$$$              end if
*     --09/25/13 21:44-lwang-end----------------------------------------*
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
*       Obtain total force for all particles due in the current block.
*
      call cputim(ttt33)
      ttnewt = ttnewt + (ttt33 - ttt32)*60.

***** --Regular force calculation--------------------------------------**
      IF(NREG.GT.0)THEN
*
#ifdef GPU
*       Send all single particles to GPU memory
         call cputim(tt51)
         NN = NTOT - IFIRST + 1
         CALL GPUNB_SEND(NN,BODY(IFIRST),X(1,IFIRST),XDOT(1,IFIRST))
         call cputim(tt52)
*     --09/26/13 16:58-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(rank.eq.0) print*,'NREG',NXTLEN,NREG,TIME,NGPUC,tt52-tt51
c$$$      NGPUC = NGPUC + 1
c$$$      if(rank.eq.0) then
c$$$         do i=1,nreg
c$$$            L=ireg(i)
c$$$            write(100,*) NREG, Name(L), STEPR(L),TIME
c$$$         end do
c$$$      end if
*     --09/26/13 16:58-lwang-end----------------------------------------*
         ttgrcomm = ttgrcomm + (tt52-tt51)*60.0
#endif

*
*       Block-Step Level Diagnostics (R.Sp. 29.Apr. 1993)
         IF(KZ(33).GT.1)JHISTR(NREG) = JHISTR(NREG) + 1
*

#ifdef PARALLEL
         IF(NREG.LE.iserreg)THEN
#endif
*
            call cputim(tt1)
#ifdef GPU
            CALL CALC_REG_ON_GPU(IREG,1,NREG)
            call cputim(tt54)
            ttgrcalc = ttgrcalc + (tt54-tt1)*60.0
#endif
*
**!$omp parallel do if(NREG.GE.ithread) private(L,I,K)
            DO L = 1,NREG
               I = IREG(L)
*     
               IMPI(1:LMAX,L) = LIST(1:LMAX,I)
*     
               CALL REGINT(I,IMPI(1,L))
*     
               F(1,I) = 0.5D0*(FI(1,I) + FR(1,I))
               F(2,I) = 0.5D0*(FI(2,I) + FR(2,I))
               F(3,I) = 0.5D0*(FI(3,I) + FR(3,I))
               FDOT(1,I) = ONE6*(FIDOT(1,I) + FRDOT(1,I))
               FDOT(2,I) = ONE6*(FIDOT(2,I) + FRDOT(2,I))
               FDOT(3,I) = ONE6*(FIDOT(3,I) + FRDOT(3,I))               
               X0(1:3,I) = XN(1:3,I)
               X0DOT(1:3,I) = XNDOT(1:3,I)
            END DO
**!$omp end parallel do
*     
            NBSUM = 0
            call cputim(tt2)
            ttreg = ttreg + (tt2-tt1)*60.

#ifdef PARALLEL
*       Start PARALLEL section for regint
         ELSE
*
            call cputim(tt998)
            call mpi_barrier(MPI_COMM_WORLD,ierr)
            call cputim(tt999)
            ibarcount=ibarcount+1
            ttbar = ttbar + (tt999-tt998)*60
            ttbarreg = ttbarreg +(tt999-tt998)*60

*
            nl = NREG
*     
            inl = nl/isize
            jsize = isize*inl
            idiff = nl - jsize
            irun = 0
*     
            do ix = 1,isize
               inum(ix)=inl
               if(ix.le.idiff)inum(ix) = inum(ix) + 1
               ista(ix) = irun+1
               if(ista(ix).gt.nl)inum(ix) = 0
               irun = irun + inum(ix)
            end do
*     
            istart = ista(rank+1)
            iend = ista(rank+1) + inum(rank+1) - 1
c$$$      if (ixxxx.le.10000000.and.time.gt.1.33.and.rank.eq.0) then
c$$$         PRINT*,' Start reg block rank ',rank,
c$$$     *        ' TIME=',TIME,' NREG=',NREG
c$$$         PRINT*,' Start reg block rank ',rank,
c$$$     *        ' istart,iend=',istart,iend
c$$$         call flush(6)
c$$$      end if
*
            NC5 = NBFULL
            NC6 = NBVOID
            NC11 = NBSMIN
            NC15 = NBDIS2
*     NC18 = NRCONV
            NC19 = NBDIS
            NC20 = NLSMIN
            NC30 = NFAST
            NCX = NICONV
            NCY = NBFAST
            NBSUM = NBFLUX
            NC12 = NBPRED
*

#ifdef GPU
            call cputim(tt1)
            CALL CALC_REG_ON_GPU(IREG,istart,iend)
            call cputim(tt56)
            ttgrcalc2 = ttgrcalc2 + (tt56-tt1)*60.0
            ttreg = ttreg + (tt56-tt1)*60.0
#endif
*
**!$omp parallel do private(L,I,K) lastprivate(IMPI)
            DO L = istart,iend
               call cputim(tt1)
               I = IREG(L)
*
               IMPI(1:LMAX,L) = LIST(1:LMAX,I)
*
               CALL REGINT(I,IMPI(1,L))
               call cputim(tt2)
               ttreg = ttreg + (tt2-tt1)*60.
*
               YMPI(37,L) = STEPR(I)
               YMPI(38,L) = STEP(I)
               YMPI(39,L) = TIMENW(I)
               YMPI(40,L) = RS(I)
               YMPI(41,L) = T0R(I)
*
               YMPI(1:3,L) = XN(1:3,I)
               YMPI(4:6,L) = XNDOT(1:3,I)
               YMPI(7:9,L) = FI(1:3,I)
               YMPI(10:12,L) = FR(1:3,I)
               YMPI(13:15,L) = FIDOT(1:3,I)
               YMPI(16:18,L) = FRDOT(1:3,I)
               YMPI(19:21,L) = D1R(1:3,I)
               YMPI(22:24,L) = D2R(1:3,I)
               YMPI(25:27,L) = D3R(1:3,I)
*     Save corrections of force polynomials from routine fpcorr
               IF (KZ(38).GT.0.OR.I.GT.N) THEN
                  YMPI(28:30,L) = D1(1:3,I)
                  YMPI(31:33,L) = D2(1:3,I)
                  YMPI(34:36,L) = D3(1:3,I)
               END IF
               CONTINUE
               call cputim(tt3)
               ttmov = ttmov +(tt3-tt2)*60.

            END DO
**!$omp end parallel do      
            call cputim(tt3)
            do j=1,isize
               if(j.eq.rank+1)then
                  JMPI(1,j) = NBFULL - NC5
                  JMPI(2,j) = NBVOID - NC6
                  JMPI(3,j) = NBSMIN - NC11
                  JMPI(4,j) = NBDIS2 - NC15
*     JMPI(5,j) = NRCONV - NC18
                  JMPI(6,j) = NBDIS - NC19
                  JMPI(7,j) = NLSMIN - NC20
                  JMPI(8,j) = NFAST - NC30
                  JMPI(9,j) = NICONV - NCX
                  JMPI(10,j) = NBFAST - NCY
                  JMPI(11,j) = NBFLUX - NBSUM
                  JMPI(12,j) = NBPRED - NC12
               else
                  JMPI(1:12,j) = 0
               end if
            END DO
            call cputim(tt4)
            ttmov = ttmov +(tt4-tt3)*60.
*
*
*        Distribute variables into private vectors again T3E (R.Sp.)
            isend = rank + 1
            if(isend.eq.isize)isend = 0
            irecv = rank - 1
            if(irecv.eq.-1)irecv = isize - 1
*
            do ir = 0,isize-2
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
c$$$      if (ixxxx.le.10000000.and.time.gt.1.33) then
c$$$         print*,' INTGRT-R: rank t',rank,time,' ir ',ir,' send ',istsen,
c$$$     *        ' thru ',istsen+icnt-1,' to ',isend,' cnt ',icnt
c$$$         print*,' INTGRT-R: rank t',rank,time,' ir ',ir,' recv ',istrec,
c$$$     *        ' thru ',istrec+icnt2-1,' fr ',irecv,' cnt2 ',icnt2
c$$$         call flush(6)
c$$$      end if
*
#ifdef PUREMPI
               call cputim(tta)
               CALL MPI_SENDRECV(YMPI(1,istsen),41*icnt,MPI_REAL8,isend,
     *              rank,YMPI(1,istrec),41*icnt2,MPI_REAL8,irecv,irecv,
     *              MPI_COMM_WORLD,status,ierr)
*     
c$$$  call mpi_barrier(MPI_COMM_WORLD,ierr)
               CALL MPI_SENDRECV(IMPI(1,istsen),lmax*icnt,MPI_INTEGER,
     *              isend,rank,IMPI(1,istrec),lmax*icnt2,MPI_INTEGER,
     *              irecv,irecv,MPI_COMM_WORLD,status,ierr)
               call cputim(ttb)
               call mpi_barrier(MPI_COMM_WORLD,ierr)
               call cputim(tt999)
               ibarcount=ibarcount+1
               ttbar = ttbar + (tt999-ttb)*60
               ttbarreg = ttbarreg +(tt999-ttb)*60
               xtsub2 = xtsub2 + dble((41*8+lmax*4)*(icnt+icnt2))
               ttsub2 = ttsub2 + (ttb-tta)*60.
#endif
*
#ifdef SHMEM
*
               call barrier()
               call shmem_get(YMPI(1,istrec),YMPI(1,istrec),41*icnt2,
     *              irecv)
               call shmem_get(IMPI(1,istrec),IMPI(1,istrec),lmax*icnt2,
     *              irecv)
#endif
*
            END DO
*
            CALL MPI_ALLGATHER(JMPI(1,rank+1),11,MPI_INTEGER,
     *           JMPILOCAL(1,1),11,MPI_INTEGER,MPI_COMM_WORLD,ierr)
*     
            call cputim(tt998)
            call mpi_barrier(MPI_COMM_WORLD,ierr)
            call cputim(tt999)
            ibarcount=ibarcount+1
            ttbar = ttbar + (tt999-tt998)*60
*!$omp parallel do if(NREG.GE.ITHREAD) private(L,I,K)     
            DO L = 1,NREG
               I = IREG(L)
               X0(1:3,I) = YMPI(1:3,L)
               X0DOT(1:3,I) = YMPI(4:6,L)
               FI(1:3,I) = YMPI(7:9,L)
               FR(1:3,I) = YMPI(10:12,L)
               FIDOT(1:3,I) = YMPI(13:15,L)
               FRDOT(1:3,I) = YMPI(16:18,L)
               D0(1:3,I) = FI(1:3,I)
               D0R(1:3,I) = FR(1:3,I)
               D1R(1:3,I) = YMPI(19:21,L)
               D2R(1:3,I) = YMPI(22:24,L)
               D3R(1:3,I) = YMPI(25:27,L)
*     Save corrections of force polynomials from routine fpcorr
               IF (KZ(38).GT.0.OR.I.GT.N) THEN
                  D1(1:3,I) = YMPI(28:30,L)
                  D2(1:3,I) = YMPI(31:33,L)
                  D3(1:3,I) = YMPI(34:36,L)
               END IF
               STEPR(I) = YMPI(37,L)
               STEP(I) = YMPI(38,L)
               TIMENW(I) = YMPI(39,L)
               RS(I) = YMPI(40,L)
               T0R(I) = YMPI(41,L)
               F(1,I) = 0.5D0*(FI(1,I) + FR(1,I))
               F(2,I) = 0.5D0*(FI(2,I) + FR(2,I))
               F(3,I) = 0.5D0*(FI(3,I) + FR(3,I))               
               FDOT(1,I) = ONE6*(FIDOT(1,I) + FRDOT(1,I))
               FDOT(2,I) = ONE6*(FIDOT(2,I) + FRDOT(2,I))
               FDOT(3,I) = ONE6*(FIDOT(3,I) + FRDOT(3,I))
            END DO
*!$omp end parallel do       
*
            NBFULL = NC5
            NBVOID = NC6
            NBSMIN = NC11
            NBDIS2 = NC15
*     NRCONV = NC18
            NBDIS = NC19
            NLSMIN = NC20
            NFAST = NC30
            NICONV = NCX
            NBFAST = NCY
            NBFLUX = NBSUM
            NBPRED = NC12
            DO J=1,isize
               NBFULL = NBFULL + JMPILOCAL(1,J)
               NBVOID = NBVOID + JMPILOCAL(2,J)
               NBSMIN = NBSMIN + JMPILOCAL(3,J)
               NBDIS2 = NBDIS2 + JMPILOCAL(4,J)
*     NRCONV = NRCONV + JMPILOCAL(5,J)
               NBDIS = NBDIS + JMPILOCAL(6,J)
               NLSMIN = NLSMIN + JMPILOCAL(7,J)
               NFAST = NFAST + JMPILOCAL(8,J)
               NICONV = NICONV + JMPILOCAL(9,J)
               NBFAST = NBFAST + JMPILOCAL(10,J)
               NBFLUX = NBFLUX + JMPILOCAL(11,J)
               NBPRED = NBPRED + JMPILOCAL(12,J)
            END DO
*     
            call cputim(tt3)
            ttmov = ttmov + (tt3-tt999)*60.
*
*      call mpi_barrier(MPI_COMM_WORLD,ierr)
*      call cputim(tt999)
*      ibarcount=ibarcount+1
*      ttbar = ttbar + (tt999-tt3)*60
         END IF
*         End PARALLEL section for regint
#endif
*
         call cputim(tt333) 
         NSTEPR = NSTEPR + NREG
         NBLCKR = NBLCKR + 1
*
**!$omp parallel do if(NREG.GE.ithread) private(L,I,K)      
         DO L = 1,NREG
            I = IREG(L)
*     
            LIST(1:LMAX,I) = IMPI(1:LMAX,L)
*
*       Check minimum neighbor sphere since last output
**!$omp critical              
            IF(LIST(1,I).GT.0)RSMIN = MIN(RSMIN,RS(I))
**!$omp end critical              
*
         END DO
**!$omp end parallel do      
         call cputim(tt334)
         ttmov = ttmov +(tt334-tt333)*60
*
*          OPEN(98,STATUS='OLD',ERR=123)
*         print*,' last reg block t=',time,' length=',nreg
*         print*,' first 10 =',(name(ireg(l)),l=1,min(nreg,10))
*         call flush(6)
*         CLOSE(98)
*123      CONTINUE
*
      END IF

*       Copy all corrected coordinates & velocities (NB! only at the end).
      call cputim(tt335)
      DO L = 1,NXTLEN
         I = NXTLST(L)
         IF (I.GT.N) THEN
            IPAIR = I - N
            IF (LIST(1,2*IPAIR-1).GT.0) NSTEPB = NSTEPB + 1
         END IF
*     
         IF(TIMENW(I).LT.TMIN)THEN
            TMIN = TIMENW(I)
            IMIN = I
         END IF
*         X(1:3,I) = X0(1:3,I)
*         XDOT(1:3,I) = X0DOT(1:3,I)
      END DO

*     --10/27/13 10:33-lwang-add----------------------------------------*
***** if ks need terminate, do prediction for icomp and jcomp-----------*
      if (iks.gt.0) then
         call xvpred(ICOMP,0)
         call xvpred(JCOMP,0)
      end if
*     --10/27/13 10:33-lwang-end----------------------------------------*
*
**!$omp parallel do if (nreg.ge.ithread) private(L,I)      
*       DO L =1,NREG
*          I = IREG(L)
*          DO 82 K = 1,3
*             X0(K,I) = XN(K,I)
*             X0DOT(K,I) = XNDOT(K,I)
* 82       CONTINUE
*       END DO
**!$omp end parallel do      
*     Needed to update only active particles, avoid sending all at reg block
*     Take care is this really efficient?
*         CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),F(1,I),FDOT(1,I),
*    &                                            BODY(I),T0(I))
      
      call cputim(tt336)
      ttcopyx =ttcopyx + (tt336-tt335)*60
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
      IF (KZ(37).GT.0.AND.LISTV(1).GT.0) THEN
          IF (DMOD(TIME,STEPM).EQ.0.0D0) THEN
              CALL SHRINK(TMIN)
              IF (LISTV(1).GT.0) THEN
                  CALL HIVEL(-1)
              END IF
          END IF
      END IF
*
*       Check optional mass loss time.
*     --09/27/13 10:39-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      IF (NREG.GT.0.or.(TIME.GT.TMDOT.AND.DMOD(TIME,STEPX).EQ.0.0D0))
c$$$     &     THEN
c$$$         if(TIME.GT.TMDOT.AND.DMOD(TIME,STEPX).EQ.0.0D0) then
c$$$            NGPUC = NGPUC + 1
c$$$         end if
c$$$         if(rank.eq.0) then
c$$$            write(140,*) NREG,NGPUC,TIME
c$$$            call flush(140)
c$$$         end if
c$$$      end if
c$$$      if(time.ge.3.5373535156250000) then
c$$$      ii=3173
c$$$*      print*,rank,'bm',ii,name(ii),f(1,ii),fdot(1,ii),body(ii),time
c$$$      print*,rank,'bm',ii,name(ii),d3(1,ii),d2r(1,ii),d3r(1,ii),time
c$$$      call flush(6)
c$$$      end if
*     --09/27/13 10:39-lwang-end----------------------------------------*
      IF (KZ(19).GT.0) THEN
*       Delay until time commensurate with 1000-year step (new polynomials).
          IF (NREG.NE.0.AND.TIME.GT.TMDOT
     *        .AND.DMOD(TIME,STEPX).EQ.0.0D0) THEN
c$$$          IF (TIME.GT.TMDOT.AND.DMOD(TIME,STEPX).EQ.0.0D0) THEN         
             call cputim(tttmdota)
             IF (KZ(19).GE.3) THEN
                CALL MDOT
             ELSE
                CALL MLOSS
             END IF
             call cputim(tttmdotb)
             ttmdot = ttmdot +(tttmdotb-tttmdota)*60
             IF (IPHASE.LT.0) GO TO 999
          END IF
      END IF
*
*     --09/27/13 10:39-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if (time.ge.3.5373535156250000) then
c$$$      ii=3173
c$$$*      print*,rank,'bm',ii,name(ii),f(1,ii),fdot(1,ii),body(ii),time
c$$$      print*,rank,'bm',ii,name(ii),d3(1,ii),d2r(1,ii),d3r(1,ii),time
c$$$      call flush(6)
c$$$      end if
c$$$      call mpi_barrier(MPI_COMM_WORLD,ierr)      
*     --09/27/13 10:39-lwang-end----------------------------------------*
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
*     Check option for general binary search.
      IF (KZ(4).NE.0.AND.TIME - TLASTS.GT.DELTAS) THEN
         CALL EVOLVE(0,0)
      END IF
*     
*     Include facility for termination of run (create dummy file STOP).
      IF(rank.EQ.0)THEN
         OPEN (99,FILE='STOP',STATUS='OLD',FORM='FORMATTED',IOSTAT=IO)
         IF (IO.EQ.0) THEN
            CLOSE (99)
*     --03/18/13 16:38-lwang-bug-fix------------------------------------*
***** Note:change nsum to nsub-----------------------------------------**
            IF (NSUB.EQ.0.and.rank.eq.0) WRITE (6,90)
*     --03/18/13 16:39-lwang-end-bug-fix--------------------------------*
 90         FORMAT  (/,9X,'TERMINATION BY MANUAL INTERVENTION')
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
*     Do not terminate during triple, quad or chain regularization.
      IF (NSUB.GT.0) THEN
*     Specify zero step to enforce termination.
         STEPS(1:NSUB) = 0.0D0
         NTIMER = NMAX
         GO TO 1
      END IF
*     
*     Terminate run with optional COMMON save.
      IF (KZ(1).NE.0) THEN
         CPUTOT = CPUTOT + TCOMP - CPU0
         CALL MYDUMP(1,1)
         if(rank.eq.0)
     &        WRITE (6,98)  TOFF, TIME, TIME+TOFF, TCOMP, CPUTOT/60.0, 
     &        ERRTOT, DETOT
 98      FORMAT (//,9X,'COMMON SAVED AT TOFF/TIME/TTOT =',1P,E16.8,
     &        '  TCOMP =',E16.8,'  CPUTOT =',E16.8,
     &                  '  ERRTOT =',F10.6,'  DETOT =',F10.6)
      END IF
*
*     Determine time interval and step numbers per time unit
      TIMINT = TIME + TOFF - TINIT
*
#ifdef PARALLEL
      IF(rank.EQ.0)THEN
#endif
         WRITE (6,195)  rank,TIMINT,NSTEPI-NIR,NSTEPB-NIB,NSTEPR-NRGL,
     &        NSTEPU-NKS
 195     FORMAT (//,I9,' INTEGRATION INTERVAL =',F8.2,3X,' NIRR=',I11,
     &        ' NIRRB=',I11,' NREG=',I11,' NKS=',I11)
         WRITE (6,196)  (NSTEPI-NIR)/TIMINT,(NSTEPB-NIB)/TIMINT,
     &        (NSTEPR-NRGL)/TIMINT,(NSTEPU-NKS)/TIMINT
 196     FORMAT (//,9X,' PER TIME UNIT: NIRR=',1P,D12.5,' NIRRB=',
     &        D12.5,' NREG=',D12.5,' NKS=',D12.5)
#ifdef PARALLEL
      END IF
#endif
#ifdef GPU
      CALL GPUNB_CLOSE
      CALL GPUPOT_CLOSE
      CALL GPUPOT_CLOSE_FLOAT
#endif
#ifdef PARALLEL
      call cputim(tt998)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call cputim(tt999)
      ibarcount=ibarcount+1
      ttbar = ttbar + (tt999-tt998)*60
      CALL MPI_ABORT(MPI_COMM_WORLD,ierr)
#endif
      STOP
*
 100  CONTINUE
*
*     Set current global time.
      TTOT = TIME + TOFF

      RETURN
*
      END
