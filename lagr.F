      SUBROUTINE LAGR(C)
c     C serves as the center of the system...
*
*
*       Lagrangian radii.
*       -----------------
*
      INCLUDE 'common1.h'
      REAL*4  R2
      integer * 4 jlist
      COMMON/WORK1/  R2(NMAX), jlist(nmax)
      integer nshells, tens
      real *8 delta
      parameter (nshells = 24,tens=8, delta = 0.3010/2)
      REAL*8  C(3),FLAGR(nshells), v2lagr(nshells)
      REAL*8  aniso(nshells)
      real * 8 v2, vt2, vr2
      real*8 dlagr(nshells), v2sum, msum, vr2sum, vt2sum
      integer nlagr(nshells), ifirst
      integer i,j, np, jmin, im, k
      real*8 rrmin,volume, vr, rsq
      
*      DATA FLAGR/-1.9,-1.7,-1.5,-1.3,-1.1,-.9,-.7,-.5,-.3,-.1/
*
c      write(6,*)'(lagr) calculate flagr...'
      do i = nshells-tens+1,nshells
         flagr(i) = log10(0.9+(i-nshells)*0.1d0)
      enddo
      do i = 1, nshells-tens
         j = nshells - tens+1 - i
         flagr(i) = -j*delta + log10(1.0d0-tens*0.1d0)
      enddo
      do i = 1, nshells
         nlagr(i) = MAX(NINT(10.0**FLAGR(i)*(n-ifirst+1)),1)
      enddo
c      write(6,*)'flagr', flagr
c      write(6,*)'nlagr', nlagr
*
*       Set square radii for NON-BH particles.
      ifirst = nbh + 1
      NP = 0
      DO 10 I = ifirst, n
          NP = NP + 1
          R2(NP) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &                                  (X(3,I) - C(3))**2
          JLIST(NP) = I
c         write(6,*) np, jlist(np), r2(np)
   10 CONTINUE
c      write(6,*)'After loop 10'
*
*       Sort square distances with respect to the centre C.
      CALL SORT1(NP,R2,JLIST)
c      do i = 1, np
c         write(6,*) i, jlist(i), r2(i)
c      enddo
*
c     calculate temperature and density
      do i=1,nshells
c         write(6,*)'shell loop ', i
         v2sum = 0
         vr2sum = 0
         vt2sum = 0
         msum = 0
         if(i .gt. 1) then
            jmin = nlagr(i-1)+1
            rrmin = r2(nlagr(i-1))
         else
            jmin = 1
            rrmin =  0
         endif
c         write(6,*)'jmin, nlagr = ', jmin, nlagr(i)
         do j =jmin, nlagr(i)
            
            im = jlist(j)
            msum = msum + body(im)
c            if(msum .le. 0.0d0) write(6,*) '<',msum, j, im, body(im)
            v2 = xdot(1,im)**2+xdot(2,im)**2+xdot(3,im)**2
            vr = xdot(1,im)*(x(1,im)-c(1)) +
     $           xdot(2,im)*(x(2,im)-c(2)) +
     $           xdot(3,im)*(x(3,im)-c(3))
            rsq = (X(1,IM) - C(1))**2 + (X(2,IM) - C(2))**2 +
     &           (X(3,IM) - C(3))**2
c            write(6,*) im, rsq
            vr2 = vr*vr/rsq
            vt2 = v2 - vr2
                           
            v2sum = v2sum + body(im)*v2
            vr2sum = vr2sum + body(im)*vr2
            vt2sum = vt2sum + body(im)*vt2

         enddo
         if(msum .gt. 0.0) then
            v2lagr(i) = v2sum/msum
            aniso(i) = 2 - vt2sum/vr2sum
            if(aniso(i) .gt. 9.99) aniso(i) = 9.99
            if(aniso(i) .lt. -9.99) aniso(i) = -9.99
         else
            v2lagr(i) = 0
            aniso(i) = 0
         endif
         
         volume = 4*3.141592*(r2(nlagr(i))**1.5 - rrmin**1.5)
     $        /0.33333
         dlagr(i) =  1d99
         if(volume .ne. 0.0 .and. msum .ne. 0.0) then
            dlagr(i) = msum / volume
         endif
c         write(6,*)'end shell loop ', i, msum, v2sum, volume
      enddo
*
      WRITE (6,40)
     &     (.5*LOG10(R2(MAX(NINT(10.0**FLAGR(K)*NP),1))),
     $     K=nshells,1,-1)
 40   FORMAT ('  LAGR:  ',(50F7.3))
      WRITE (6,41)(sqrt(v2lagr(k)),  K=nshells,1,-1)
 41   FORMAT ('  SIGMA:  ',(50F7.3))
      WRITE (6,42)(log10(dlagr(k)),  K=nshells,1,-1)
 42   FORMAT ('  RHOR:  ',(50F7.3))
      WRITE (6,43)(aniso(k),  K=nshells,1,-1)
 43   FORMAT (' ANISO:  ',(50F7.3))
*

*
      RETURN
*
      END
