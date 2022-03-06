      REAL*8 FUNCTION TSTEP(F,FDOT,F2DOT,F3DOT,ETA)
*
*
*       General time-step criterion.
*       ----------------------------
*
      REAL*8  F(3),FDOT(3),F2DOT(3),F3DOT(3), eta
      REAL*8 f2, fdot2, f2dot2, f3dot2, tstep2
c     automatic f2, fdot2, f2dot2, f3dot2, tstep2
*
*       Obtain new integration step using composite expression.
*       STEP = (ETA*(F*F2DOT + FDOT**2)/(FDOT*F3DOT + F2DOT**2))**0.5.
*
      F2 = F(1)**2 + F(2)**2 + F(3)**2
      FDOT2 = FDOT(1)**2 + FDOT(2)**2 + FDOT(3)**2
      F2DOT2 = F2DOT(1)**2 + F2DOT(2)**2 + F2DOT(3)**2
      F3DOT2 = F3DOT(1)**2 + F3DOT(2)**2 + F3DOT(3)**2
*
c - standard scheme
c      TSTEP = (SQRT(F2*F2DOT2) + FDOT2)/(SQRT(FDOT2*F3DOT2) + F2DOT2)
c      TSTEP = SQRT(ETA*TSTEP)

c - returns dt^2
      TSTEP = (SQRT(F2*F2DOT2) + FDOT2)/(SQRT(FDOT2*F3DOT2) + F2DOT2)
      TSTEP = ETA*TSTEP
c - test for unstable shrinkage of the timestep
      tstep2 = eta*f2/fdot2 * 0.001
      if(tstep .lt. tstep2) then
         tstep = tstep2
      endif

c     - returns (roughly) dt^4 
c      TSTEP = ((F2*F2DOT2) + FDOT2**2)/((FDOT2*F3DOT2) + F2DOT2**2)
c      TSTEP = (ETA*eta*TSTEP)
*
      RETURN
*
      END
