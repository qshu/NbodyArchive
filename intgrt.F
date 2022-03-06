c#define HARP3DEBUG 1
c#define HARP2TEST 1      
c
c   NOTE : precent code DOES NOT work on HARP2 at all...
c

/*#define PEC2*/
#ifdef HARP2
#    undef HARP2
#    define HARP2 1
#else
#    define HARP2 0
#endif
#ifdef HARP3
#    undef HARP3
#    define HARP3 1
#else
#    define HARP3 0
#endif
#ifdef GRAPE6
#    undef GRAPE6
#    define GRAPE6 1
#else
#    define GRAPE6 0
#endif


#if (HARP2 || HARP3 || GRAPE6)
#define USE_GRAPE 1      
#   ifdef HARP2TEST
#      define HOST_FORCE 1
#   else
#      define HOST_FORCE 0
#   endif      
#else
#   define USE_GRAPE 0      
#   define HOST_FORCE 1
#endif      
      SUBROUTINE accel(i, acc, adot, phii)
*
*
*       ------------------
*
      INCLUDE 'common1.h'
      integer i
      integer j,  k
c     automatic j,  k
      real*8 acc(3), adot(3), phii
      real*8  dr2, dx(3), dv(3)
      real*8 drdv, dr2i, dr3i,dri
c     automatic  dr2, dx, dv
c     automatic drdv, dr2i, dr3i
      
            DO 10 K = 1,3
               acc(k) = 0.0
               adot(k) = 0.
 10         CONTINUE
            phii = 0.0
            DO 15 J = 1,N
               IF (J.NE.I) THEN
*              
	       if(j.le.nbh .and. i.le.nbh) then
                   dr2 = 0
               else if(j.le.nbh .or. i.le.nbh) then
                   dr2 = meps2
               else
                   dr2 = EPS2
               endif
               drdv = 0.
               DO 11 K = 1,3
                  dx(K) = X(K,J) - X(K,I)
                  dv(K) = XDOT(K,J) - XDOT(K,I)
                  dr2 = dr2 + dx(k)**2
                  drdv = drdv + dx(k)*dv(k)
 11            CONTINUE
*              
               dr2i = 1.0/dr2
               dri =  body(j)*sqrt(dr2i)
               phii = phii - dri
               dr3i = dr2i*dri
               drdv = 3.*drdv*dr2i
*              
               DO 12 K = 1,3
                  acc(k) = acc(k) + dx(K)*dr3i
                  adot(k) =adot(k) + (dv(K) - dx(K)*drdv)*dr3i
 12            CONTINUE
               endif
 15         CONTINUE
      end

      SUBROUTINE accel0(i, acc, adot, phii, jend)
*
*
*       ------------------
*
      INCLUDE 'common1.h'
      integer i, jend
      integer j,  k
c     automatic j,  k
      real*8 acc(3), adot(3), phii
      real*8  dr2, dx(3), dv(3)
      real*8 drdv, dr2i, dr3i,dri
c     automatic  dr2, dx, dv
c     automatic drdv, dr2i, dr3i
            DO 15 J = 1, jend
               IF (J.NE.I) THEN
	       if(j.le.nbh .and. i.le.nbh) then
                   dr2 = 0
               else if(j.le.nbh .or. i.le.nbh) then
                   dr2 = meps2
               else
                   dr2 = EPS2
               endif
               drdv = 0.
               DO 11 K = 1,3
                  dx(K) = X(K,J) - X(K,I)
                  dv(K) = XDOT(K,J) - XDOT(K,I)
                  dr2 = dr2 + dx(k)**2
                  drdv = drdv + dx(k)*dv(k)
 11            CONTINUE
*              
               dr2i = 1.0/dr2
               dri =  BODY(J)*SQRT(dr2i)
               phii = phii - dri
               dr3i = dr2i*dri
               drdv = 3.*drdv*dr2i
*              
               DO 12 K = 1,3
                  acc(k) = acc(k) + dx(K)*dr3i
                  adot(k) =adot(k) + (dv(K) - dx(K)*drdv)*dr3i
 12            CONTINUE
               endif
 15         CONTINUE
      end

      SUBROUTINE accel2(i, acc, adot, phii, jfirst,jend)
*
*
*       ------------------
*
      INCLUDE 'common1.h'
      integer i, jend,jfirst
      integer j,  k
c     automatic j,  k
      real*8 acc(3), adot(3), phii
      real*8  dr2, dx(3), dv(3)
      real*8 drdv, dr2i, dr3i,dri
c     automatic  dr2, dx, dv
c     automatic drdv, dr2i, dr3i
            DO 15 J = jfirst, jend
               IF (J.NE.I) THEN
	       if(j.le.nbh .and. i.le.nbh) then
                   dr2 = 0
               else if(j.le.nbh .or. i.le.nbh) then
                   dr2 = meps2
               else
                   dr2 = EPS2
               endif
               drdv = 0.
               DO 11 K = 1,3
                  dx(K) = X(K,J) - X(K,I)
                  dv(K) = XDOT(K,J) - XDOT(K,I)
                  dr2 = dr2 + dx(k)**2
                  drdv = drdv + dx(k)*dv(k)
 11            CONTINUE
*              
               dr2i = 1.0/dr2
               dri =  BODY(J)*SQRT(dr2i)
               phii = phii - dri
               dr3i = dr2i*dri
               drdv = 3.*drdv*dr2i
*              
               DO 12 K = 1,3
                  acc(k) = acc(k) + dx(K)*dr3i
                  adot(k) =adot(k) + (dv(K) - dx(K)*drdv)*dr3i
 12            CONTINUE
               endif
 15         CONTINUE
      end

      SUBROUTINE PREDICT(npred)
*
*
*       N-body integrator.
*       ------------------
*
      INCLUDE 'common1.h'
       real*8  s, s1, s2
      integer j, npred
c$par doall
         DO 5 J = 1,Npred
            S = TIME - T0(J)
            S1 = 1.5*S
            S2 = 2.0*S
            X(1,J) = ((FDOT(1,J)*S + F(1,J))*S + X0DOT(1,J))*S + X0(1,J)
            X(2,J) = ((FDOT(2,J)*S + F(2,J))*S + X0DOT(2,J))*S + X0(2,J)
            X(3,J) = ((FDOT(3,J)*S + F(3,J))*S + X0DOT(3,J))*S + X0(3,J)
            XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
            XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
            XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
#if HARP2TEST
            write(6,*)'i,x,v', j, X(1,J),X(2,J),X(3,J)
            write(6,*)'     ', XDOT(1,J),XDOT(2,J),XDOT(3,J)
#endif
    5    CONTINUE

      end
*
#if USE_GRAPE
      subroutine predict_for_harp(nxtlen,nxtlst, xharp, vharp,
     $                   aharp,jharp,pharp, tharp, mharp, jpmax)
      INCLUDE 'common1.h'
      INTEGER  NXTLEN,NXTLST(nmax), ii, j, ip, jpfirst
       real*8  s, s1, s2, ssq, scu
      integer  npipemax, jpmax, k
      parameter (npipemax=nmax)
      real * 8 xharp(3,npipemax),vharp(3,npipemax)
      real * 8 aharp(3,npipemax),jharp(3,npipemax)

#if USE_GRAPE
      real * 8 mharp(npipemax),tharp(npipemax)
      integer  pharp(npipemax)
#endif
         ip = 0
         jpfirst = 1
         DO 1100 ii = 1, nxtlen
            ip = ip + 1
            j = nxtlst(ii)
            S = TIME - T0(J)
            ssq = s*s
            scu = ssq*s
            S1 = 3.0*Ssq
            S2 = 2.0*S
            do 1120 k = 1, 3
            X(k,j) = FDOT(k,J)*scu+F(k,J)*ssq+X0DOT(k,J)*S
     $              +X0(k,J)
            xdot(k,j) = FDOT(k,J)*S1 + F(k,J)*S2 + X0DOT(k,J)
 1120       continue
            j = nxtlst(ii)
            xharp(1,ii)=x(1,j)
            xharp(2,ii)=x(2,j)
            xharp(3,ii)=x(3,j)
            vharp(1,ii)=xdot(1,j)
            vharp(2,ii)=xdot(2,j)
            vharp(3,ii)=xdot(3,j)
c            write(6,*)ii,vharp(1,ii),vharp(2,ii),vharp(3,ii)
#if HARP2            
            call h2setxvt(j-1,x(1,j),xdot(1,j),time)
#endif
#if HARP3
            tharp(j)=time
            if(ip .eq. jpmax .or. ii .eq. nxtlen) then
                call h3mjpdma_indirect(ip,nxtlst(jpfirst),
     $                     x,xdot,f,fdot,mharp,tharp,1,0)
                call h3mjpdma_start(0)
                call h3wait
                ip = 0
                jpfirst = ii+1
            endif
#endif            
 1100    CONTINUE
      end
#endif


      subroutine corrector(nxtlen,nxtlst, xpred, vpred,firr, fd,ph,
     $                    mharp,jpmax,stepmin,
     $                    ierror)
      INCLUDE 'common1.h'
      INTEGER  NXTLEN,NXTLST(nmax), ii,  k, i,jpmax,ip, jpfirst
      real*8   dt, dtold, dtsq, dt2, dt6,stepmin
      real*8   at3xc, bt2xc, at3vc, bt2vc
      real*8   fd6, sum, at3, bt2, df
      real*8 step2, step4, ttmp, dt0, tstep
      REAL*8  FIRR(3,nmax),fd(3,nmax),ph(nmax), d2t(3),d3t(3)
      real * 8 xpred(3,nmax), vpred(3,nmax)
      integer ierror
c     automatic   dt, dtold, dtsq, dt2, dt6
c     automatic  at3xc, bt2xc, at3vc, bt2vc
c     automatic   fd6, sum, at3, bt2, df
c     automatic step2, step4, ttmp, dt0
      real * 8 j2(3)
      
#if USE_GRAPE
      integer  npipemax, id
      parameter (npipemax=nmax)
#if (HARP3 || GRAPE6)
      real * 8 mharp(npipemax)
#endif
#else
      real *8 mharp
#endif
      integer nsmall
      data nsmall /0/
      data j2 /0.0,0.0,0.0/
      dtold = 0.0
      ip = 0
      ierror = 0
      jpfirst = 1
      do 1010 ii = 1, nxtlen
         ip = ip + 1
         i = nxtlst(ii)
         dt = time - t0(i)
         if(dtold .ne. dt) then
            dtsq = dt**2
            dt6 = 6.0/(dt*dtsq)
            dt2 = 2.0/dtsq
*           
            at3xc = 0.05*dtsq
            bt2xc = one12*dtsq
            at3vc = 0.25*dt
            bt2vc = one3*dt
         endif
         t0(i) = time
         do 20 k = 1,3
            df = 2.0*f(k,i) - firr(k,ii)
            
            fd6 = 6.0*fdot(k,i)
            sum = fd6 + fd(k,ii)
            
            at3 = 2.0*df + dt*sum
            bt2 = -3.0*df - dt*(sum + fd6)
*           
            x(k,i) = xpred(k,i) + at3xc*at3 + bt2xc*bt2
            xdot(k,i) = vpred(k,i) + at3vc*at3 + bt2vc*bt2
*           
            x0(k,i) = x(k,i)
            x0dot(k,i) = xdot(k,i)
*           
            f(k,i) = 0.5*firr(k,ii)
            fdot(k,i) = one6*fd(k,ii)
            phi(i)=ph(ii)
*           
c            d3(k,i) = at3*dt6
c            d2(k,i) = bt2*dt2 + dt*d3(k,i)
            d3t(k) = at3*dt6
            d2t(k) = bt2*dt2 + dt*d3t(k)
*           NOTE: These are real derivatives!
 20      continue
#define USE_INDIRECT_JPDMA               

*              
*        Specify new time-step (use fast function instead of routine STEPI).
c         TTMP = TSTEP(FIRR(1,ii),FD(1,ii),D2(1,I),D3(1,I),ETA)
         TTMP = TSTEP(FIRR(1,ii),FD(1,ii),D2t,D3t,ETA)
         if(i.le.nbh) ttmp = ttmp*0.33
         DT0 = TTMP
         if(ttmp .lt. stepmin*stepmin*256) then
c           very small timestep - see what is going on
            write(6,*)'small step, i, t, dt = ', i, t0(i),
     $           step(i),ttmp
            write(6,*)'F :',FIRR(1,ii),FIRR(2,ii),FIRR(3,ii)
            write(6,*)'F1:',Fd(1,ii),Fd(2,ii),Fd(3,ii)
c           write(6,*)'d2:',d2(1,i),d2(2,i),d2(3,i)
c           write(6,*)'d3:',d3(1,i),d3(2,i),d3(3,i)
            nsmall = nsmall + 1
            if(nsmall .gt. 200) then
#if HARP2      
               call h2free
#endif         
#if HARP3      
               call h3close
#endif         
#if GRAPE6     
               call g6_close(0)
#endif         
               write(6,*)'Small steps : restart at ',time
               call mydump(0,2)
               write(6,*)'Small steps : new time ',time
               ierror = 1
               return
            endif
            
         endif
         step2 = step(i)**2
         step4 = step2**2
*        
*        Select discrete value .
         if (ttmp .gt. 4.0*step2) then
            if (dmod(time,2.0*step(i)) .eq. 0d0) then 
               ttmp = MIN(2.0*step(i),1.0d0)
            else
               ttmp = step(i) 
            end if
         else if (ttmp .lt. step2 
     $           .and.ttmp.gt.stepmin*stepmin) then
            dt0 = sqrt(dt0)
            call stepk(dt0,ttmp)
         else
            ttmp = step(i)
         end if
         
*        
*        Set new block step, update next time and save current time.
         if(ttmp .gt. deltat) ttmp = ttmp * 0.5
         step(i) = ttmp
         timenw(i) = step(i) + t0(i)
         NSTEPN = NSTEPN + 1
         if(i .le. nbh) nstepbh = nstepbh + 1
#if GRAPE6
         id = i  - 1
         call g6_set_j_particle(0,id, i,
     $        t0(i), step(i), mharp(i), j2, fdot(1,i),
     $        f(1,i), x0dot(1,i), x0(1,i ))
#endif   
#if HARP3
         if(ip .eq. jpmax .or. ii .eq. nxtlen) then
            call h3mjpdma_indirect(ip,nxtlst(jpfirst),
     $           x0,x0dot,f,fdot,mharp,t0,1,0)
            call h3mjpdma_start(0)
            call h3wait
            ip = 0
            jpfirst = ii+1
         endif
#endif   
 1010 continue
      end      
*
      SUBROUTINE INTGRT
*
*
*       N-body integrator.
*       ------------------
*
      INCLUDE 'common1.h'
      REAL*8  FIRR(3,nmax),fd(3,nmax)
      INTEGER  NXTLEN,ILOC,NXTLST(nmax)
       integer li,i,j,k,ii,iii
       real*8 tmin
       real * 8 xpred(3,nmax), vpred(3,nmax),ph(nmax), epsinv
       real*8 stepmin
      integer ierror
      real *8 j2(3)
      integer nfield, jpmax, npipemax
      parameter (npipemax=nmax)
      real * 8 mharp(npipemax),tharp(npipemax)
#if USE_GRAPE
      integer npipe, np, ip,itmp, id, nf, ipmove
      real * 8 xharp(3,npipemax),vharp(3,npipemax)
      integer h2getnpipe
      REAL*8 fh(3,npipemax),ptmp(npipemax),
     $     fdh(3,npipemax), demax, fe, fdemax, fde
      REAL*8 zeroeps2(npipemax)
      real * 8 aharp(3,npipemax),jharp(3,npipemax)
      real * 8 eps2harp(npipemax),h2harp(npipemax)
      integer  pharp(npipemax)
      integer h3npipe, h3jpmax,  nj, ifirst,ncalc
      integer calldum
#if GRAPE6
      REAL*8 fg6old(3,n), jg6old(3,n)
#endif
#endif
*
      do k = 1, 3
         j2(k) = 0
      enddo
      nfield = n - nbh
 2000 call printcpu
      epsinv = 1.0d0/sqrt(eps2)
      stepmin = 2.0**(-50)/sqrt(body(1))
#if USE_GRAPE
#ifdef HARP3DEBUG
      write(6,*)'call h3open'
      call flush(6)
#endif
#if HARP3      
      call h3open
      call h3setmode(255,15)
      npipe = h3npipe()
#endif
#if GRAPE6
      call g6_open(0)
      npipe = 48
#endif      
      ifirst = 1
#ifdef HARP3DEBUG
      write(6,*)'call h3setmode'
      call flush(6)
#endif
#ifdef HARP3DEBUG
      write(6,*)' return h3setmode'
      call flush(6)
#endif
      write(6,*)'Number of pipes = ', npipe
#ifdef HARP2TEST
      write(6,*)'Number of pipes = ', npipe
#endif
      do 1042 i = 1,n
         h2harp(i) =0.0
         eps2harp(i) =eps2*0.5
         pharp(i) = i
         aharp(1,i)=f(1,i)
         aharp(2,i)=f(2,i)
         aharp(3,i)=f(3,i)
         jharp(1,i)=fdot(1,i)
         jharp(2,i)=fdot(2,i)
         jharp(3,i)=fdot(3,i)
         mharp(i) = body(i)
         if (i .le. nbh)          mharp(i) = 0.0
 1042 continue
#if HARP3

      do 1041 i = n,n+npipe
         h2harp(i) =0.0
         eps2harp(i) =eps2*0.5
         zeroeps2(i) =meps2*0.5
 1041 continue
      
      jpmax = h3jpmax()
      do 1044 iii = 1, 2
      do 1040 i = 1, n, jpmax
         nj = min(jpmax,n-i+1)
         call h3mjpdma_indirect(nj,pharp(i),x0,x0dot,f,fdot,
     $        mharp,t0, 1,0)
         call h3mjpdma_start(0)
         call h3wait
 1040 continue
 1044 continue
#endif      
#if GRAPE6
      do i = 1, n
         call g6_set_j_particle(0,i-1, i,
     $        t0(i), step(i), mharp(i), j2, fdot(1,i),
     $        f(1,i), x0dot(1,i), x0(1,i))
      enddo
#endif
#endif
*     Initialize block times.
      IF (TIME.EQ.0D0) THEN
          DO 1000 I = 1,N
              TIMENW(I) = T0(I) + STEP(I)
 1000     CONTINUE
      END IF
*
      NXTLEN = 0
      TMIN = TPREV
*
*       Determine all particles due at next block time.
 1001 LI = 0
    1	continue
	if(nxtlen .gt. 0) nxtlen = - nxtlen
c
c       i gives the first particle to move, nxtlen gives
c       the number of particles to be moved.
c      write(6,*) 'before: i,nxtlen = ', i, nxtlen, ipos
 	call inext(i,nxtlen,iloc,nxtlst,timenw,tmin)
c      write(6,*) 'after:i,nxtlen = ', i, nxtlen, ipos
      TIME = T0(I) + STEP(I)
c      write(6,*) i, time, t0(i), step(i)
*
*       Check next output time before beginning a new block.
      IF (TIME.GT.TNEXT) THEN
         write(6,*) 'intgrt, ', time, tnext
         TIME = TPREV
         GO TO 30
      END IF
*
*       Predict all coordinates and velocities to first order.
         NBLOCK = NBLOCK + 1
#if HOST_FORCE         
         call predict(n)
#endif         
#if USE_GRAPE
         call predict_for_harp(nxtlen,nxtlst, xharp, vharp,
     $         aharp,jharp, pharp, tharp, mharp, jpmax)
#if GRAPE6
c        set fh, fdh and ph ...
         do ii = 1, nxtlen
            j = nxtlst(ii)
            do k = 1, 3
               fh(k,ii) = f(k,j)
               fdh(k,ii) = fdot(k,j)
            enddo
            ph(ii) = phi(j)
            eps2harp(ii) = eps2
            if (j .le. nbh) eps2harp(ii) = meps2
         enddo
#endif   
c        
c        prdict the point-mass particles anyway...
c
             call predict(nbh)
#endif
         DO 1110 ii = 1, nxtlen
            j = nxtlst(ii)
            do 1120 k = 1, 3
               xpred(k,j)=x(k,j)
               vpred(k,j)=xdot(k,j)
 1120       continue
 1110    CONTINUE
         
 1500    continue
C        in the case of hardware error, redo the force calculation
*        

*       Obtain the current force and derivative.
c         
#if HOST_FORCE      
c$par doall
         do 1053 ii = 1, nxtlen
            i = nxtlst(ii)
            call accel(i, firr(1,ii), fd(1,ii),ph(ii))
 1053    continue
#endif      
#if USE_GRAPE
#if HARP3         
            call h3wait
#ifdef HARP3DEBUG
      write(6,*)'call h3setti ', time
      call flush(6)
#endif
            call h3setti(time)
            call h3wait
c           write(6,*)'before h3calc', ii, xharp(1,ii), npipe 
#endif         
#if GRAPE6
            call g6_set_ti(0,time)
c           write(6,*)'before h3calc', ii, xharp(1,ii), npipe 
#endif         
         do 1050 ii = 1, nxtlen, npipe
            np = nxtlen - ii + 1
            if(np .gt. npipe) np = npipe
#if HARP3
#ifdef HARP3DEBUG
      write(6,*)'call h3calc ', npipe
      call flush(6)
#endif
c      ncalc = ((np+5)/6)*6
      ncalc = npipe
            call h3calc(n,ncalc, xharp(1,ii),vharp(1,ii),eps2harp(ii),
     $                  h2harp(ii),
     $           fh(1,ii),fdh(1,ii),ph(ii))
            if (ifirst .eq. 1) then
               ifirst = 0
               call h3calc(n,npipe, xharp(1,ii),vharp(1,ii),
     $              eps2harp(ii),
     $              h2harp(ii),
     $              fh(1,ii),fdh(1,ii),ph(ii))
            endif
#endif
#if GRAPE6
#ifdef HARP3DEBUG
            write(6,*)'call h3calc ', npipe
            call flush(6)
#endif
            ncalc = np
            if (ncalc .lt. 48) ncalc = 48
C
C     The above is to see if the hardware error reduces...
C     actually, it does not look effective, though
C            
            call g6calc_firsthalf0(0,n, ncalc, nxtlst(ii),
     $           xharp(1,ii),
     $           vharp(1,ii), fh(1,ii), fdh(1,ii), ph(ii),
     $           eps2harp(ii), h2harp(ii),0)
            call g6calc_lasthalf(0, n, ncalc,nxtlst(ii),
     $           xharp(1,ii),
     $           vharp(1,ii), eps2harp(ii), h2harp(ii), fh(1,ii),
     $           fdh(1,ii), ph(ii))
*              write(6,*) 'GRAPE-6 hardware error'
*              call g6_reset(0)
*              call g6_reset_fofpga(0)
*              call g6_close(0)
*              call g6_open(0)
*              call g6_set_ti(0,time)
*              do i = 1, n
*                 call g6_set_j_particle(0,i-1, i,
*    $                 t0(i), step(i), mharp(i), j2,
*    $                 fdot(1,i), f(1,i), x0dot(1,i), x0(1,i))
*              enddo
*              
*              goto 1500
*              write(6,*) 'error recovery not implemented yet'
*              stop
*           endif
#endif

            do 1060 ip=1,np
               i=ii+ip-1
#ifdef HARP2TEST            
               write(6,600) i,(firr(k,i),k=1,3),
     $            (fh(k,i),k=1,3),
     $    (fd(k,i),k=1,3),
     $    (fdh(k,i),k=1,3)
 600           format('Particle ', i5, /,
     $              'host a:', 3g20.10,/,'harp a:',3g20.10,/,
     $              'host j:', 3g20.10,/,'harp j:',3g20.10)
               fe = ((fh(1,i)-firr(1,i))**2
     $              +(fh(2,i)-firr(2,i))**2
     $              +(fh(3,i)-firr(3,i))**2)
     $              /(fh(1,i)**2+fh(2,i)**2+fh(3,i)**2)
               fde = ((fdh(1,i)-fd(1,i))**2
     $              +(fdh(2,i)-fd(2,i))**2
     $              +(fdh(3,i)-fd(3,i))**2)
     $              /(fdh(1,i)**2+fdh(2,i)**2+fdh(3,i)**2)
               write(6,601) sqrt(fe), sqrt(fde)
 601           format('Rel errors = ', 2g15.6)
#endif            
               do 1070 k = 1, 3
                  firr(k,i) = fh(k,i)
                  fd(k,i) = fdh(k,i)
 1070          continue
 1060       continue
 1050    continue
#if HARP3
         do ii = 1, nxtlen
            i = nxtlst(ii)
            if(i .le. nbh) then
               if(step(i) .gt. 1e-11) then
#if HARP3         
                  call h3wait
                  call h3calc(n,npipe, xharp(1,ii),vharp(1,ii),
     $                 zeroeps2(1),
     $                 h2harp(ii),
     $                 fh(1,1),fdh(1,1),ptmp)

                  do 1058 k = 1, 3
                     firr(k,ii) = fh(k,1)
                     fd(k,ii) = fdh(k,1)
 1058             continue
                  ph(ii) = ptmp(1)
#endif
#if GRAPE6                  
                  do 1058 k = 1, 3
                     firr(k,ii) = fh(k,ii)
                     fd(k,ii) = fdh(k,ii)
 1058             continue
#endif
               else
                  write(6,*)'BH force on host', i, time

                  call predict(n)
                  call accel2(i, firr(1,ii), fd(1,ii),ph(ii),
     $                 nbh+1,n)
               endif
            else
#if HARP3
               ph(ii) = ph(ii)+ body(i)*epsinv
#endif
            endif
         enddo
#endif   
c        
c        include the contribution of point-mass particles
c        $par doall
         do 1052 ii = 1, nxtlen
            i = nxtlst(ii)
            nf = nbh
c            if (i .eq. 3) then
c               write(6,666)'after h3ca', i,ii,
c     $              (firr(k,ii),k=1,3),(fd(k,ii),k=1,3)
 666           format(A, 2i5, 6e15.6)
c            endif
            call accel0(i, firr(1,ii), fd(1,ii),ph(ii),nf)
c            if (i.eq.3) then
c               write(6,666)'after acc0', i,ii,
c     $              (firr(k,ii),k=1,3),(fd(k,ii),k=1,3)
c              call predict(n)
c              call accel(i, firr(1,ii), fd(1,ii),ph(ii))
c              write(6,666)'after acc ', i,ii,
c              $              (firr(k,ii),k=1,3),(fd(k,ii),k=1,3)
c            endif
 1052    continue
         
#endif   
*        
*        Include the corrector and set new T0, F, FDOT, D2 & D3.
         call corrector(nxtlen,nxtlst, xpred, vpred,firr, fd,ph,
     $          mharp,jpmax,stepmin,  ierror)
         if(ierror .ne. 0) goto 2000
         tprev = time
         GO TO 1
*        
 30      continue
         call printcpu
#if HARP2
         call h2free
#endif   
#if GRAPE6
         call g6_close(0)
#endif   
#if HARP3
         call h3close
#endif   
         RETURN
*        
         end
#ifdef AHO
      subroutine touch_dummy
      INCLUDE 'common1.h'
      integer i
      do i = 1, nmax
        zzzzdummy(i) = i
      end do
      end
      
#endif

