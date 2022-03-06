      SUBROUTINE INCLIN(X,V,XCM,VCM,ALPHA)
*
*
*       Inclination of hierarchical system.
*       -----------------------------------
*
      REAL*8  A12,A22,A1A2,FAC,ALPHA
      REAL*8  XCM(3),VCM(3),X(3,3),V(3,3),A0(3),A2(3)
*
*
*       Define inner binary I1 & I2 and outer component I3.
      I1 = 1
      I2 = 2
      I3 = 3
*
*       Evaluate scalar products of angular momenta.
      A12 = 0.0
      A22 = 0.0
      A1A2 = 0.0
      DO 20 K = 1,3
          K1 = K + 1
          IF (K1.GT.3) K1 = 1
          K2 = K1 + 1
          IF (K2.GT.3) K2 = 1
          A0(K) = (X(K1,I1) - X(K1,I2))*(V(K2,I1) - V(K2,I2))
     &          - (X(K2,I1) - X(K2,I2))*(V(K1,I1) - V(K1,I2))
          A2(K) = (X(K1,I3) - XCM(K1))*(V(K2,I3) - VCM(K2))
     &          - (X(K2,I3) - XCM(K2))*(V(K1,I3) - VCM(K1))
          A12 = A12 + A0(K)**2
          A22 = A22 + A2(K)**2
          A1A2 = A1A2 + A0(K)*A2(K)
   20 CONTINUE
*
*       Determine inclination in radians within both limits.
      IF (A12.EQ.0.0D0.OR.A22.EQ.0.0D0) THEN
          FAC = 0.0
      WRITE (89,30)  A12, A22
   30 FORMAT (' INCLIN ZERO    A12 A22 ',1P,2E10.2)
      CALL FLUSH(89)
      ELSE
          FAC = A1A2/SQRT(A12*A22)
      END IF
      IF (FAC.GE.1.0d0) THEN
          ALPHA = 0d0
      ELSE IF (FAC.LE.-1.0d0) THEN
          ALPHA = ACOS(-1.0d0)
      ELSE
          ALPHA = ACOS(FAC)
      ENDIF
*
      RETURN
*
      END
