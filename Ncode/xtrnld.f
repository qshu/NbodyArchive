      SUBROUTINE XTRNLD(I1,I2,KCASE)
*
*
*       External force & derivatives.
*       -----------------------------
*
      INCLUDE 'common6.h'
*
*
*       See whether to include the standard or point-mass tidal force.
      IF (KZ(14).LE.2.AND.KCASE.EQ.1) THEN
*       Include tidal force & first derivative (I1 = I2 for single body).
          DO 10 I = I1,I2
              FI(1,I) = FI(1,I) + TIDAL(4)*XDOT(2,I)
              FI(2,I) = FI(2,I) - TIDAL(4)*XDOT(1,I)
              FR(1,I) = FR(1,I) + TIDAL(1)*X(1,I)
              FR(3,I) = FR(3,I) + TIDAL(3)*X(3,I)
              D1(1,I) = D1(1,I) + TIDAL(4)*(FI(2,I) + FR(2,I))
              D1(2,I) = D1(2,I) - TIDAL(4)*(FI(1,I) + FR(1,I))
              D1R(1,I) = D1R(1,I) + TIDAL(1)*XDOT(1,I)
              D1R(3,I) = D1R(3,I) + TIDAL(3)*XDOT(3,I)
   10     CONTINUE
      END IF
*
      IF (KZ(14).LE.2.AND.KCASE.EQ.2) THEN
*       Add the second and third derivatives due to the tidal field.
          DO 20 I = I1,I2
              D2(1,I) = D2(1,I) + TIDAL(4)*FDOT(2,I)
              D2(2,I) = D2(2,I) - TIDAL(4)*FDOT(1,I)
              D2R(1,I) = D2R(1,I) + TIDAL(1)*F(1,I)
              D2R(3,I) = D2R(3,I) + TIDAL(3)*F(3,I)
              D3(1,I) = D3(1,I) + TIDAL(4)*(D2(2,I) + D2R(2,I))
              D3(2,I) = D3(2,I) - TIDAL(4)*(D2(1,I) + D2R(1,I))
              D3R(1,I) = D3R(1,I) + TIDAL(1)*FDOT(1,I)
              D3R(3,I) = D3R(3,I) + TIDAL(3)*FDOT(3,I)
   20     CONTINUE
      END IF
*
      RETURN
*
      END
