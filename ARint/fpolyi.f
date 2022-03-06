      SUBROUTINE FPOLYI(I)
*
*
*       Initialization of higher derivatives.
*       -------------------------------------
*
      INCLUDE 'common6.h'
*
*
      DO 10 K = 1,3
          D2(K,I) = 0.0
          D3(K,I) = 0.0
          D2R(K,I) = 0.0
          D3R(K,I) = 0.0
   10 CONTINUE
*
      STEP(I) = SMAX/16.0D0
      STEPR(I) = SMAX/16.0D0
      T0(I) = TBLOCK
      T0R(I) = TBLOCK
*
      CALL STEPS(I,I)
*
      RETURN
*
      END
