#ifdef HARP3
#define USE_GRAPE 1
#endif      
#ifdef GRAPE6
#define USE_GRAPE 1
#endif      
#ifndef USE_GRAPE
#define USE_GRAPE 0
#endif      

      SUBROUTINE START
*
*
*       Polynomial initialization.
*       --------------------------
*
      INCLUDE 'common1.h'
*
*
*       COMMON variables
*       ****************
*
*       ------------------------------------------------------------------
*       BE      Total energy (BE(1-3) = initial, previous & current value).
*       BODY    Mass of a particle.
*       CPU     Maximum computing time in minutes.
*       DELTAT  Output time interval in units of the crossing time.
*       DTLIST  Interval for updating the time-step list.
*       D1      First divided force difference.
*       D2      Second divided force difference.
*       D3      Third divided force difference.
*       EPS2    Soft potential parameter (square is saved).
*       ETA     STEP = (ETA*(F*F2DOT+FDOT**2)/(FDOT*F3DOT+F2DOT**2))**0.5.
*       F       One-half the total force at last step. 
*       FDOT    One-sixth the force derivative at last step.
*       KZ      Options for alternative paths at run time (see table).
*       N       Total particle number.
*       NSTEPN  Total number of individual time-steps.
*       ONE3    The constant 1.0/3.0 (ONE6, ONE9 & ONE12 are similar).
*       QE      Relative energy tolerance for error control (option 2).
*       RSCALE  Half-mass radius determined from potential energy.
*       STEP    Time-step of a particle.
*       TCR     Standard crossing time.
*       TCRIT   Termination time in units of the standard crossing time.
*       TIME    Physical integration time in scaled units.
*       TLIST   Time for next updating of the time-step list.
*       TNEXT   Time for next output.
*       T0      Time of the last force calculation.
*       T1      Time of first backwards force calculation.
*       T2      Time of second backwards force calculation.
*       T3      Time of third backwards force calculation.
*       X       Coordinates at current time and at output.
*       XDOT    Velocities at beginning of step and at output time.
*       X0      Coordinates at the beginning of a step.
*       X0DOT   Velocities at the beginning of a step.
*       ZMASS   Total mass.
*       ------------------------------------------------------------------
*
*
*       Initialize times, integration step counter & constants.
      TIME = 0.0
      TNEXT = 0.0
      TLIST = 0.0
      NSTEPN = 0
      NSTEPBH = 0
      ONE3 = 1.0/3.0D0
      ONE6 = 1.0/6.0D0
      ONE9 = 1.0/9.0D0
      ONE12 = 1.0/12.0D0
*
*       Read input parameters.
      CALL INPUT
*
*       Set initial conditions: BODY(I), X(K,I), XDOT(K,I); I=1,N & K=1,3.
      CALL DATA
      write(6,*)'exit data'
*
*       Scale initial conditions to new units.
      if (kz(4) .eq. 0) CALL mySCALE
*
*       Initialize fixed block steps (32 levels).
      CALL IBLOCK
*

#ifndef USE_GRAPE      
*       Obtain force & first derivative.
      CALL FPOLY1(1,N)
*
*       Obtain second & third force derivatives.
*     CALL FPOLY2(1,N)
*
*       Determine initial steps.
      CALL STEPS(1,N)
#else
      call fpoly0
#endif      
*
      RETURN
*
      END

#ifdef HARP3      
            SUBROUTINE FPOLY0
*
*
*       Force & first derivative on HARP-2.
*       -----------------------------------
*
      INCLUDE 'common1.h'
#include "harp3.h"      
      integer  jpmax, jsend
      integer i, j, ii, idum, npipe, ncpert, k, nn, np, ip, jj, j1,
     $     ifirst, nfield
      real * 8 eps
      eps = sqrt(eps2)
#ifndef MONOPOLIZEHARP
      write(6,*) 'call h3open'
      call h3open()
      write(6,*) 'return  h3open'
      call h3wait
#endif
      call h3setti(time)
      call h3setmode(255,15)
      npipe = h3npipe()
      write(6,*)'fpoly0, eps2 = ', eps2, ' npipe = ', npipe
      do  i = 1,npipe
         h3eps2(i) = eps2*0.5
      enddo
      
*
      ifirst = nbh + 1
      nfield = n - nbh
      do i = 1,nfield
         h3index(i) = i
      enddo
      do i = 1,n
         do k = 1, 3
            f(k,i) = 0.0
            fdot(k,i) = 0.0 
         enddo
         t0(i) = time
      enddo
      jpmax = h3jpmax()
      do i = 1, nfield, jpmax
         jsend = min(nfield - i + 1, jpmax)
#if 0         
         call h3jpdma_indirect(jsend,h3index(i),x(1,ifirst),
     $        xdot(1,ifirst), f(1,ifirst),fdot(1,ifirst),
     $        body(ifirst),t0(ifirst),1)
#endif
         call h3mjpdma_indirect(jsend,h3index(i),x(1,ifirst),
     $        xdot(1,ifirst), f(1,ifirst),fdot(1,ifirst),
     $        body(ifirst),t0(ifirst),1,0)
         call h3mjpdma_start(0)
         call h3wait
      enddo
         
*
*       Calculate F & FDOT on all single particles using HARP-2.
      I = 1
       DO 50 II = 1,N,NPIPE
          NP = N - II + 1
          IF (NP.GT.NPIPE) NP = NPIPE
          CALL H3CALC(Nfield,NPIPE,X(1,II),XDOT(1,II),h3eps2, h3h2,
     &         h3acc,h3jerk,h3pot)
          DO 40 IP = 1,NP
              JJ = II + IP - 1
              if(i .ne. jj) then
                 write(6,*)' fpoly0 , i <> jj with ', i, jj, ii, ip
              endif
*       Copy F & FDOT into COMMON variables.
              DO 35 K = 1,3
                  F(K,I) = h3acc(K,IP)
                  FDOT(K,I) = h3jerk(K,IP)
   35         CONTINUE
 
              phi(i) = h3pot(ip) + body(i)/eps
              
              I = I + 1
   40     CONTINUE
   50 CONTINUE
#ifndef MONOPOLIZEHARP
      call h3close
#endif
*
      do i = 1, n
         call fpoly1a(i)
c         write(6,666) i, phi(i), (f(k,i), k = 1, 3)
c 666     format(i5, 4e16.7)
      enddo
*       Obtain new time-step for all single particles.
      CALL mySTEPS(1,N,1)
*
      RETURN
*
      END
#endif
#ifdef GRAPE6
            SUBROUTINE FPOLY0
*
*
*       Force & first derivative on HARP-2.
*       -----------------------------------
*
      INCLUDE 'common1.h'
#include "harp3.h"      
      integer  jpmax, jsend
      integer i, j, ii, idum, npipe, ncpert, k, nn, np, ip, jj, j1,
     $     ifirst, nfield, ic, js, je, kk
      real * 8 eps
      real * 8 j2(3)
      integer iloc, g6_npipes
      eps = sqrt(eps2)
      write(6,*) 'call g6open'
      do i = 1, g6ncluster
         call g6_open(i-1)
      enddo
      write(6,*) 'return  g6open'
 1500 continue
      do i = 1, g6ncluster 
         call g6_set_ti(i-1,time)
      enddo
      npipe = g6_npipes()
      write(6,*)'fpoly0, eps2 = ', eps2, ' npipe = ', npipe
*
      ifirst = nbh + 1
      nfield = n - nbh
      do i = 1,n
         h3index(i) = i
      enddo
      do k = 1, 3
         j2(k) = 0.0
      enddo
      do i = 1,n
         do k = 1, 3
            f(k,i) = 0.0
            fdot(k,i) = 0.0 
         enddo
         t0(i) = time
         step(i)=1.0/1024
      enddo
      do  i = 1,npipe
         h3eps2(i) = eps2
         h3h2(i) = 0.0d0
      enddo
      j1 = (n-ifirst+g6ncluster)/g6ncluster
      if (g6ncluster .eq. 4) then
         j1 = (n-ifirst+g6ncluster)/g6ncluster * 64.0/62.0
      endif
      do ic = 1, g6ncluster
         js = ifirst+(ic-1)*j1
         je = js + j1 -1
         if (je .gt. n) je = n
         g6nj(ic) = je-js+1
         write(6,*)' ic, nj  = ', ic, g6nj(ic)
         do i = js, je
            iloc = i-js
            g6jcid(i)=ic-1
            g6jcloc(i)=iloc
            call g6_set_j_particle(g6jcid(i),g6jcloc(i), i,
     $           t0(i), step(i), body(i), j2, fdot(1,i),
     $           f(1,i), xdot(1,i), x(1,i))
         enddo
      enddo
      
*     
*       Calculate F & FDOT on all single particles using HARP-2.
*
*     initialize f, fdot and phi with some reasonable guess...
*      
      
      do i = 1,n
         do k = 1, 3
            f(k,i) = 0.01
            fdot(k,i) = 100.0
         enddo
         phi(i) = 100.0
      enddo
      do kk = 1, 3
         DO  II = 1,N,NPIPE
            NP = N - II + 1
            IF (NP.GT.NPIPE) NP = NPIPE
            do ic = 1, g6ncluster
               do ip = 1, np
                  do k = 1, 3
                     g6fwork(k,ip,ic) = f(k,ii+ip-1)
                     g6jwork(k,ip,ic) = fdot(k,ii+ip-1)
                  enddo
                  g6pwork(ip,ic) = phi(ii+ip-1)
               enddo
            enddo
            do ic = 1, g6ncluster
               call g6calc_firsthalf0(ic-1,g6nj(ic), np, h3index(ii),
     $              x(1,ii), xdot(1,ii), g6fwork(1,1,ic),
     $              g6jwork(1,1, ic), g6pwork(1,ic), h3eps2,h3h2,0)
            enddo
#if 0
            do ic = 1, g6ncluster
               call g6_set_njp(ic-1, g6nj(ic))
            enddo
#endif            
            do ic = 1, g6ncluster
               call g6calc_lasthalf(ic-1, g6nj(ic), np,h3index(ii),
     $              x(1,ii),xdot(1,ii), h3eps2, h3h2, g6fwork(1,1,ic),
     $              g6jwork(1,1, ic), g6pwork(1,ic),0)
*                 write(6,*) 'GRAPE-6 hardware error cl=', ic
*                 do jj = 1, g6ncluster
*                    call g6_reset(jj-1)
*                 enddo
*                 goto 1500
*              endif
            enddo
            do ip = 1, np
               do k = 1, 3
                  f(k,ii+ip-1) = g6fwork(k,ip,1) 
                  fdot(k,ii+ip-1) =g6jwork(k,ip,1) 
               enddo
               phi(ii+ip-1) = g6pwork(ip,1)
            enddo
            do ic = 2, g6ncluster 
               do ip = 1, np
                  do k = 1, 3
                     f(k,ii+ip-1) = f(k,ii+ip-1) + g6fwork(k,ip,ic) 
                     fdot(k,ii+ip-1) = fdot(k,ii+ip-1) +g6jwork(k,ip,ic) 
                  enddo
                  phi(ii+ip-1) = phi(ii+ip-1) +g6pwork(ip,ic)
               enddo
            enddo
         enddo
C
C check possible  underflow...
C            
            do i = 1, n
               if ((fdot(1,i) .eq. 0.0d0) .and.
     $              (fdot(2,i) .eq. 0.0d0) .and.
     $              (fdot(3,i) .eq. 0.0d0)) then
                  fdot(1,i) = 1e4*(1e-4 ** kk)
                  if (kk .eq. 3)write(6,*) 'underflow ', i, kk
               endif
            enddo
      enddo
      do i = 1, g6ncluster
         call g6_close(i-1)
      enddo
*      
       do i = 1, n
          call fpoly1a(i)
c         write(6,666) i, phi(i), (f(k,i), k = 1, 3)
c         666     format(i5, 4e16.7)
       enddo
*      Obtain new time-step for all single particles.
       CALL mySTEPS(1,N,1)
*      
       RETURN
*      
       END
#endif
      
      
      
#ifdef USE_GRAPE
      SUBROUTINE mySTEPS(I1,I2,KCASE)
*
*
*       Initialization of time-steps & prediction variables.
*       ----------------------------------------------------
*
      INCLUDE 'common1.h'
      integer i, i1, i2, iter, k, kcase
      real*8 fi2, fd2, dt, dtn
*
*
*       Set new steps and initialize prediction variables.
      DO 40 I = I1,I2
*
*       Obtain time-step using DT = 0.5*ETA*F/FDOT (D2 & D3 not available).
          FI2 = F(1,I)**2 + F(2,I)**2 + F(3,I)**2
          FD2 = FDOT(1,I)**2 + FDOT(2,I)**2 + FDOT(3,I)**2
          if(fd2 .ne. 0.0d0) then
             DT = 0.125*ETA*SQRT(FI2/FD2)
          else
             DT = 0.125*ETA*SQRT(FI2)
             write(6,*) '(steps) ZERO FDOT for particle ', i
          endif
c          write(6,*)i, fi2, fd2, dt
*       Convert predicted step to nearest block time-step (truncated down).
          dtn = 0.0078125
          do while (dtn .gt. dt)
             dtn = dtn * 0.5
          enddo
          IF (TIME.LE.0.0D0) THEN
              STEP(I) = DTN
          ELSE 
*       Reduce step by factor 2 until commensurate with current time.
              STEP(I) = DTN
              ITER = 0
   10         IF (DMOD(TIME,STEP(I)).NE.0.0D0) THEN
                  STEP(I) = 0.5D0*STEP(I)
                  ITER = ITER + 1
                  IF (ITER.LT.12) GO TO 10
                  WRITE (6,15)  I, ITER, TIME/STEP(I), DT, STEP(I)
   15             FORMAT (' WARNING!   I ITER T/STEP DT STEP ',
     &                                 I5,I4,F16.4,1P,2E9.1)
                  IF (STEP(I).GT.1e-13) GO TO 10
                  STOP
              END IF
          END IF
*
*
*       Set prediction variables
          DO 30 K = 1,3
              X0(K,I) = X(K,I)
              X0dot(K,I) = Xdot(K,I)
              F(K,I) = 0.5D0*F(K,I)
              FDOT(K,I) = ONE6*FDOT(K,I)
#if ( ! USE_GRAPE)
              D2(K,I) = 0.0
              D3(K,I) = 0.0
#endif
   30     CONTINUE
   40 CONTINUE
*
      RETURN
*
      END
#endif
      
