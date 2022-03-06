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
C      automatic j,  k
      real*8 acc(3), adot(3), phii
      real*8  dr2, dx(3), dv(3)
      real*8 drdv, dr2i, dr3i,dri
C      automatic  dr2, dx, dv
C      automatic drdv, dr2i, dr3i
      
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
C      automatic j,  k
      real*8 acc(3), adot(3), phii
      real*8  dr2, dx(3), dv(3)
      real*8 drdv, dr2i, dr3i,dri
C      automatic  dr2, dx, dv
C      automatic drdv, dr2i, dr3i
c*$* noconcurrentize
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
C      automatic j,  k
      real*8 acc(3), adot(3), phii
      real*8  dr2, dx(3), dv(3)
      real*8 drdv, dr2i, dr3i,dri
C      automatic  dr2, dx, dv
C      automatic drdv, dr2i, dr3i
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

      real * 8 mharp(npipemax),tharp(npipemax)
      integer  pharp(npipemax)
         ip = 0
         jpfirst = 1
c*$* concurrentize
         DO 1100 ii = 1, nxtlen
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
            xharp(1,ii)=x(1,j)
            xharp(2,ii)=x(2,j)
            xharp(3,ii)=x(3,j)
            vharp(1,ii)=xdot(1,j)
            vharp(2,ii)=xdot(2,j)
            vharp(3,ii)=xdot(3,j)
c            write(6,*)ii,vharp(1,ii),vharp(2,ii),vharp(3,ii)
 1100    CONTINUE
      end
#endif


      subroutine corrector(nxtlen,nxtlst, xpred, vpred,firr, fd,ph,
     $                    mharp,jpmax,stepmin,
     $                    ic, ierror)
      INCLUDE 'common1.h'
#if USE_GRAPE
      INCLUDE 'harp3.h'
#endif
      INTEGER  NXTLEN,NXTLST(nmax), ii,  k, i,jpmax,ip, jpfirst
      real*8   dt, dtold, dtsq, dt2, dt6,stepmin
      real*8   at3xc, bt2xc, at3vc, bt2vc
      real*8   fd6, sum, at3, bt2, df
      real*8 step2, step4, ttmp, dt0, tstep
      REAL*8  FIRR(3,nmax),fd(3,nmax),ph(nmax), d2t(3),d3t(3)
      real * 8 xpred(3,nmax), vpred(3,nmax)
      integer ierror, ic
C      automatic   dt, dtold, dtsq, dt2, dt6
C      automatic  at3xc, bt2xc, at3vc, bt2vc
C      automatic   fd6, sum, at3, bt2, df
C      automatic step2, step4, ttmp, dt0
      
#if USE_GRAPE
      integer  npipemax, id
      parameter (npipemax=nmax)
#if (HARP3 || GRAPE6)
      real * 8 mharp(npipemax)
      real * 8 j2(3)
#endif
#else
      real mharp
#endif
      integer nsmall
      data nsmall /0/
      data j2 /0.0,0.0,0.0/
      dtold = 0.0
      ip = 0
      ierror = 0
      jpfirst = 1
c*$* noconcurrentize
      do 1010 ii = 1, nxtlen
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
c              d3(k,i) = at3*dt6
c              d2(k,i) = bt2*dt2 + dt*d3(k,i)
               d3t(k) = at3*dt6
               d2t(k) = bt2*dt2 + dt*d3t(k)
*              NOTE: These are real derivatives!
 20         continue
#define USE_INDIRECT_JPDMA
            
*           
*           Specify new time-step (use fast function instead of routine STEPI).
c           TTMP = TSTEP(FIRR(1,ii),FD(1,ii),D2(1,I),D3(1,I),ETA)
            TTMP = TSTEP(FIRR(1,ii),FD(1,ii),D2t,D3t,ETA)
            if(i.le.nbh) ttmp = ttmp*0.5
            if (ttmp .gt. 0.03125) ttmp = 0.03125
            DT0 = TTMP
            step2 = step(i)**2
            step4 = step2**2
*           
*           Select discrete value .
            if (ttmp .gt. 4.0*step2) then
               if (dmod(time,2.0*step(i)) .eq. 0d0) then 
                  ttmp = MIN(2.0*step(i),1.0d0)
               else
                  ttmp = step(i) 
               end if
            else if (ttmp .lt. step2 
     $              .and.ttmp.gt.stepmin*stepmin) then
               dt0 = sqrt(dt0)
               call stepk(dt0,ttmp)
            else
               ttmp = step(i)
            end if
            
*           
*           Set new block step, update next time and save current time.
            if(ttmp .gt. deltat) ttmp = ttmp * 0.5
            step(i) = ttmp
            timenw(i) = step(i) + t0(i)
#if GRAPE6 
            if (i .gt. nbh) then
               id = i  - 1 - nbh
               call g6_set_j_particle(g6jcid(i),g6jcloc(i), i,
     $              t0(i), step(i), mharp(i), j2, fdot(1,i),
     $              f(1,i), x0dot(1,i), x0(1,i ))
            endif
#endif      
#ifdef SMP_HOST
         endif
#endif
 1010 continue
#if GRAPE6 
      call g6_flush_jp_buffer(0)

#endif      
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
#if USE_GRAPE
      INCLUDE 'harp3.h'
      integer npipe, npipemax, np, ip,itmp, id, nf, ipmove,ic,jj
      parameter (npipemax=nmax)
      real * 8 xharp(3,npipemax),vharp(3,npipemax)
      integer h2getnpipe
      REAL*8 fh(3,npipemax),ptmp(npipemax),
     $     fdh(3,npipemax), demax, fe, fdemax, fde
      REAL*8 zeroeps2(npipemax)
      real * 8 aharp(3,npipemax),jharp(3,npipemax)
      real * 8 mharp(npipemax),tharp(npipemax)
      real * 8 eps2harp(npipemax),h2harp(npipemax)
      integer  pharp(npipemax)
      integer jpmax, nj, ifirst,ncalc
      integer calldum
      real *8 j2(3)
      integer nfield
#if GRAPE6
      REAL*8 fg6old(3,n), jg6old(3,n)
      integer g6_npipes
      integer iinext, npnext, ncalcnext
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
#if GRAPE6
      ic=1
      call g6_open(ic-1)
*     call g6_setfoerrchkmode(ic-1, 1)
      call g6_initialize_jp_buffer(ic-1,n+128)

      npipe = g6_npipes()
#endif      
      ifirst = 1
      write(6,*)'Number of pipes = ', npipe
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
#if GRAPE6
      do i = nbh+1, n
         call g6_set_j_particle(g6jcid(i),g6jcloc(i), i,
     $        t0(i), step(i), mharp(i), j2, fdot(1,i),
     $        f(1,i), x0dot(1,i), x0(1,i))
      enddo
      write(6,*)'end sending j'
#endif
#endif
*     Initialize block times.
c      IF (TIME.EQ.0D0) THEN
      DO 1000 I = 1,N
         TIMENW(I) = T0(I) + STEP(I)
 1000 CONTINUE
c      END IF
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

#if GRAPE6
      /* test code for non-overlap call ... */
      call g6_set_ti(0,time)

      do 1050 ii = 1, nxtlen, npipe
         np = nxtlen - ii + 1
         if(np .gt. npipe) np = npipe
         iinext = ii + np
         npnext = nxtlen - iinext + 1
         if(npnext .gt. npipe) npnext = npipe
         
#ifdef HARP3DEBUG
         write(6,*)'call h3calc ', npipe
         call flush(6)
#endif   
         ncalc = np
         ncalcnext = npnext
         call g6calc_firsthalf0(0,g6nj(ic), ncalc, 
     $           nxtlst(ii),xharp(1,ii),
     $           vharp(1,ii), fh(1,ii), fdh(1,ii), ph(ii),
     $           eps2harp(ii), h2harp(ii),0)
         ic = 1
         call g6calc_lasthalf(0, g6nj(ic), ncalc,
     $        nxtlst(ii),xharp(1,ii),
     $        vharp(1,ii), eps2harp(ii), h2harp(ii),
     $        g6fwork(1,1,ic),g6jwork(1,1,ic),
     $        g6pwork(1,ic))
*           g6error(ic) = 1
*        else
*           g6error(ic) = 0
*        endif
*        if (g6error(ic) .eq. 1) then
*           write(6,*) 'GRAPE-6 hard error on cluster', ic
*           jj=1
*           call g6_reset(jj-1)
*           The following is added on May 23.... JM
*           call g6_reset_fofpga(jj-1)
*           call g6_close(jj-1)
*           call g6_open(jj-1)
*           call g6_set_ti(jj-1,time)
*           do i = nbh+1, n
*              call g6_set_j_particle(g6jcid(i),g6jcloc(i), i,
*    $              t0(i), step(i), mharp(i), j2,
*    $              fdot(1,i), f(1,i), x0dot(1,i), x0(1,i))
*           enddo
*           call g6_flush_jp_buffer(0)
*           
*           goto 1500
*        endif
         do ip = 1, np
            do k = 1, 3
               fh(k,ii+ip-1) = g6fwork(k,ip,1) 
               fdh(k,ii+ip-1) =g6jwork(k,ip,1) 
            enddo
            ph(ii+ip-1) = g6pwork(ip,1)
         enddo
         do 1060 ip=1,np
            i=ii+ip-1
            do 1070 k = 1, 3
               firr(k,i) = fh(k,i)
               fd(k,i) = fdh(k,i)
 1070       continue
 1060    continue
 1050 continue
#endif
c        
c        include the contribution of point-mass particles

         NSTEPN = NSTEPN + nxtlen
         do ii = 1, nxtlen
            i = nxtlst(ii)
            if(i .le. nbh) nstepbh = nstepbh + 1
         enddo
c$par doall
c*$* assert concurrent call
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
#ifndef SMP_HOST         
         call corrector(nxtlen,nxtlst, xpred, vpred,firr, fd,ph,
     $          mharp,jpmax,stepmin, 0, ierror)
#else

         call corrector(nxtlen,nxtlst, xpred, vpred,firr, fd,ph,
     $        mharp,jpmax,stepmin, 0, ierror)
#endif         
         
         if(ierror .ne. 0) then
#ifdef GRAPE6
               call g6_close(0)
#endif            
            goto 2000
         endif
         tprev = time
         GO TO 1
*        
 30      continue
         call printcpu
#if GRAPE6
         call g6_close(0)
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

