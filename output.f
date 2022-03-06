


      SUBROUTINE OUTPUT(mode)
*
*
*       Energy check and output.
*       ------------------------
*
      INCLUDE 'common1.h'
      integer i,k, mode
      real * 8 zkin, pot, de, q, tc, de1, etaold
      integer idump
      data idump /0/
*
*
      mode = 0
*       Obtain the total kinetic & potential energy at current time.
      CALL ENERGY(ZKIN,POT)
*
*       Form approximate half-mass radius.
      RSCALE = 0.5*ZMASS**2/POT
*
*       Set total energy and define standard crossing time.
      IF (TIME.GT.0.0) BE(2) = BE(3)
      BE(3) = ZKIN - POT
      TCR = ZMASS**2.5/(2.0*ABS(BE(3)))**1.5
*
*       Obtain the relative energy error (save initial energy).
      IF (TIME.EQ.0.0D0) THEN
          DE = 0.0
          BE(1) = BE(3)
      ELSE
          DE = (BE(3) - BE(2))/MAX(ZKIN,POT)
*       Relative energy error with respect to kinetic energy.
      END IF
      DE1 = (BE(3) - BE(1))/BE(1)
*
      Q = ZKIN/POT
      TC = TIME/TCR
      WRITE (6,10)  TIME, Q, NSTEPN, NBLOCK, nstepbh, DE,DE1, BE(3)
   10 FORMAT (/,' T =',F9.3,'  Q =',F5.2,'  STEPS =',3I11,/
     &	'  DE= ',1p,e10.3,' ',e10.3,0p,  '  E =',F9.5)
      call diag
      CALL FLUSH(6)
*       Check optional error control. BEFORE CREATING SNAPSHOT
      IF (KZ(2).GT.0.AND.ABS(DE).GT.80.0*QE) THEN
         if(mod(kz(2),8)/4 .ne. 0) then
c           restarting ....
            etaold =  eta
            call mydump(0,2)
            eta = etaold * 0.5
            write(6,*)'Restart, time, new eta = ', time,eta
            do i = 1, n
               if(t0(i) + 0.5*step(i) .ge. time)then
                  step(i) = 0.5*step(i)
               endif
            enddo
            mode = 1
            return
         else
            
            WRITE  (6,90)
 90         FORMAT (/,9X,'CALCULATIONS HALTED * * *')
            STOP
         endif
      END IF
      if(mod(kz(2),4)/2 .ne. 0) then
c        adjust ETA
         if(ABS(DE) .gt. QE)then
            eta = eta*sqrt(QE/abs(de))
            write(6,*)'new eta = ', eta
         else if(ABS(DE) .lt. 0.1*qe .and. eta .lt. 0.04)then
            eta = eta * 1.1
            write(6,*)'new eta = ', eta
         endif
      endif
            
*
*

*       Check optional output of single bodies & binaries.
      CALL BODIES
*
      TNEXT = TNEXT + DELTAT
*
*       Restore single precision velocities for coordinate predictions.
      DO 80 I = 1,N
          DO 70 K = 1,3
              XDOT(K,I) = X0DOT(K,I)
   70     CONTINUE
   80 CONTINUE
      if(kz(3) .gt. 0) then
          if(mod(idump,kz(3)) .eq. 0) then
              call mydump(1,2)
          endif
c          if(mod(idump,kz(3)) .eq. 1) then
c              call mydump(1,3)
c          endif
          idump = idump + 1
      endif
*
      RETURN
*
      END


      subroutine diag
      include 'common1.h'
      integer i,k,k1,k2
      integer ndim
      parameter (ndim=3)
      REAL cm(ndim), mtot, cmv(ndim)
      REAL am(ndim)
      real * 8 c(3), mbh
      do 10 k = 1,ndim
        cm(k)=0.0
        cmv(k)=0.0
 10   continue
      mtot=0.0
      do 20 k=1,ndim		
        do 30 i=1,n
          cm(k)=cm(k)+body(i)*x(k,i)
          cmv(k)=cmv(k)+body(i)*xdot(k,i)
 30     continue
 20   continue
      do 40 i=1,n
        mtot=mtot+body(i)
40    continue
      do 70 k=1,ndim		
        cm(k)=cm(k)/mtot
        cmv(k)=cmv(k)/mtot
 70   continue
      do 170 k = 1,ndim
         am(k) = 0
 170  continue
      do 180 k = 1, ndim
         k1 = mod(k,3)+1
         k2 = mod(k+1,3)+1
         do 190 i = 1,n
            am(k) = am(k)
     $           + body(i)*((x(k1,i)-cm(k1))*(xdot(k2,i)-cmv(k2))
     $           - (x(k2,i)-cm(k2))*(xdot(k1,i)-cmv(k1)))
 190     continue
 180  continue
      write(6,601)'CM ', cm
      write(6,601)'CMV', cmv
      write(6,601)'AM ', am
 601  format(' ',a3, ' : ', 3g15.6)
      call outbinary
      if(mod(kz(7),2) .ne. 0) then
c        print the lagrad infor around geometric center
         do k = 1, 3
            c(k) = 0
         enddo
         call lagr(c)
      endif
      if(mod(kz(7),4)/2 .ne. 0 .and. nbh .gt. 0) then
c        print the lagrad infor around BH cm
         do k = 1, 3
            c(k) = 0
         enddo
         mbh = 0
         do i = 1, nbh
            mbh = mbh + body(i)
            do k = 1, 3
               c(k) = c(k) + body(i)*x(k,i)
            enddo
         enddo
         do k = 1, 3
            c(k) = c(k)/mbh
         enddo
         call lagr(c)
      endif
      end
