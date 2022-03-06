      SUBROUTINE XVPREDX(I,NNB)
*
*
*       Prediction of coordinates & velocities.
*       ---------------------------------------
*
      INCLUDE 'common6.h'
*
      S = TIME - T0(I)
      IF (S.EQ.0.d0) THEN
          X(1,I) = X0(1,I)
          X(2,I) = X0(2,I)
          X(3,I) = X0(3,I)
          XDOT(1,I) = X0DOT(1,I)
          XDOT(2,I) = X0DOT(2,I)
          XDOT(3,I) = X0DOT(3,I)
          GO TO 40
      END IF
*
      IF (NNB.EQ.0.AND.I.LE.N) THEN
*       Adopt low order for standard single particles (I <= N & NNB = 0).
          S1 = 1.5*S
          S2 = 2.0*S
          X(1,I) = ((FDOT(1,I)*S + F(1,I))*S + X0DOT(1,I))*S + X0(1,I)
          X(2,I) = ((FDOT(2,I)*S + F(2,I))*S + X0DOT(2,I))*S + X0(2,I)
          X(3,I) = ((FDOT(3,I)*S + F(3,I))*S + X0DOT(3,I))*S + X0(3,I)
          XDOT(1,I) = (FDOT(1,I)*S1 + F(1,I))*S2 + X0DOT(1,I)
          XDOT(2,I) = (FDOT(2,I)*S1 + F(2,I))*S2 + X0DOT(2,I)
          XDOT(3,I) = (FDOT(3,I)*S1 + F(3,I))*S2 + X0DOT(3,I)
      ELSE
*       Predict coordinates & velocities of body #I to order F3DOT.
          IF (BODY(I).EQ.0.0D0) GO TO 40
          A1 = 0.05d0*S
          A2 = 0.25*S
          A3 = T0(I) - T0R(I)
          A4 = 0.5*A3
*
          DO 30 K = 1,3
              FK = ((ONE6*D3R(K,I)*A3 + 0.5*D2R(K,I))*A3 + D1R(K,I))*A3
     &                                              + D0R(K,I) + D0(K,I)
              F1DOTK = (D3R(K,I)*A4 + D2R(K,I))*A3 + D1R(K,I) + D1(K,I)
              F2DOTK = 0.5*(D3R(K,I)*A3 + D2R(K,I) + D2(K,I))
              F3DOTK = ONE6*(D3R(K,I) + D3(K,I))
              X(K,I) = ((((F3DOTK*A1 + ONE12*F2DOTK)*S +
     &                   ONE6*F1DOTK)*S + 0.5*FK)*S + X0DOT(K,I))*S +
     &                   X0(K,I)
              XDOT(K,I)  = (((F3DOTK*A2 + ONE3*F2DOTK)*S +
     &                            0.5*F1DOTK)*S + FK)*S + X0DOT(K,I)
   30     CONTINUE
      END IF
*
   40 RETURN
*
      END
