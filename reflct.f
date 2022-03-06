      SUBROUTINE REFLCT(I,RI2)
*
*
*       Boundary reflection.
*       --------------------
*
      INCLUDE 'common4.h'
*
*
      RDOT = X(1,I)*XDOT(1,I) + X(2,I)*XDOT(2,I) + X(3,I)*XDOT(3,I)
*
*       Avoid multiple reflections for inward moving particles.
      IF (RDOT.LT.0.0) GO TO 50
*
*       Obtain radial velocity assuming small density centre displacement.
      RI = SQRT(RI2)
      RDOT = RDOT/RI
*
*       Form radial velocity for mirror reflection.
      RDOT = -2.0*RDOT
*
*       Set reflected velocity at the boundary and initialize X0DOT.
      DO 10 K = 1,3
          XDOT(K,I) = XDOT(K,I) + RDOT*X(K,I)/RI
          X0DOT(K,I) = XDOT(K,I)
   10 CONTINUE
*
*       Predict current coordinates & velocities for all particles.
      T0(I) = TIME
      CALL XVPRED(IFIRST,NTOT)
*
*       Initialize force polynomial & time-step.
      CALL FPOLY1(I,I,0)
*
      RI2 = -1.0
*       Return of negative argument indicates new polynomials & steps.
*
      NBREF = NBREF + 1
*
   50 RETURN
*
      END
