      SUBROUTINE ENERGY(ZKIN,POT)
*
*
*       Total energy.
*       -------------
*
      INCLUDE 'common6.h'
*
*       Calculate the potential energy.
      IFIRST = 1
      NTOT = N
      ZKIN = 0.D00
      POT = 0.0
      RAV = 0.D0
      PAV = 0.D0
      PMAV = 0.D0
*
      DO 20 I = 1,NTOT
      JMIN = I + 1
      POTJ = 0.D00
      POTI = 0.D00
*       POTI contains potential at particles position to be stored later (R.Sp.)
*
      DO 30 J = 1,N
      IF (J.EQ.I) GO TO 30
          A1 = X(1,I) - X(1,J)
          A2 = X(2,I) - X(2,J)
          A3 = X(3,I) - X(3,J)
      A4 = BODY(J)/DSQRT (A1*A1 + A2*A2 + A3*A3)
      POTI = POTI - A4
*  also J.LT.N?
      IF(J.GE.JMIN)THEN
      POTJ = POTJ + A4
      RAV = RAV + DSQRT (A1*A1 + A2*A2 + A3*A3)
      PAV = PAV + A4
      PMAV = PMAV + BODY(I)*A4
      END IF
   30 CONTINUE
*       Store potential in shared vector first (R.Sp.)
      PHIDBL(I) = POTI
      POT = POT + BODY(I)*POTJ
   20 CONTINUE
*
      XNN = DBLE(N*(N-1))/2.D0
      RAV = RAV/XNN
      PAV = PAV/XNN
      PMAV = PMAV/XNN
      PRINT*,' RAV,PAV,PMAV=',RAV,PAV,PMAV
*       Sum the kinetic energy (include c.m. bodies but not components).
      DO 40 I = IFIRST,NTOT
          ZKIN = ZKIN + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 +
     &                                          XDOT(3,I)**2)
   40 CONTINUE
*
      ZKIN = 0.5D0*ZKIN
*
      RETURN
*
      END
