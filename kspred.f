      SUBROUTINE KSPRED(IPAIR,I1,I,BODYIN,UI,UIDOT,XI)
*
*
*       Prediction for KS regularization.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      REAL  W1,W2,W3,RCRIT2
      REAL*8  UI(4),UIDOT(4),XI(6)
*
*       Add body #I to the perturber list for prediction.
      NNB2 = LIST(1,I1) + 2
      LIST(NNB2,I1) = I
*
*       Predict coordinates of perturbers & current c.m. to order FDOT.
      DO 10 L = 2,NNB2
          J = LIST(L,I1)
          S = TIME - T0(J)
          X(1,J) = ((FDOT(1,J)*S + F(1,J))*S + X0DOT(1,J))*S + X0(1,J)
          X(2,J) = ((FDOT(2,J)*S + F(2,J))*S + X0DOT(2,J))*S + X0(2,J)
          X(3,J) = ((FDOT(3,J)*S + F(3,J))*S + X0DOT(3,J))*S + X0(3,J)
   10 CONTINUE
*
      IF (GAMMA(IPAIR).LT.0.0001.AND.STEP(I).GT.SMIN) GO TO 20
*
*       Improve prediction of perturbers & c.m. for moderate perturbations.
      RCRIT2 = 400.0*R(IPAIR)**2
      DO 15 L = 2,NNB2
          J = LIST(L,I1)
          W1 = X(1,J) - X(1,I)
          W2 = X(2,J) - X(2,I)
          W3 = X(3,J) - X(3,I)
*       Skip distant perturbers (RIJ > 20*R) but include c.m. itself.
          IF (W1*W1 + W2*W2 + W3*W3.LT.RCRIT2) THEN
              S = TIME - T0(J)
              S1 = 0.2*S
              S4 = 0.04166666666667D0*(S*S)**2
              S5 = TIME - T0R(I)
*
              DO 12 K = 1,3
                  F2DOTK = D3R(K,J)*S5 + D2R(K,J) + D3(K,J)*S + D2(K,J)
                  X(K,J) = F2DOTK*S4 + X(K,J)
   12         CONTINUE
          END IF
   15 CONTINUE
*
*       Start integration of regularized motion.
   20 DTU = DTAU(IPAIR)
      T1PR = T0U(IPAIR) - T1U(IPAIR)
      T2PR = T0U(IPAIR) - T2U(IPAIR)
      T12PR = T1PR + T2PR
      DT06 = 0.6D0*DTU
      DT19 = ONE9*DTU
      DT12 = ONE12*DTU
      DT34 = 0.75D0*DTU
      DT32 = 1.5D0*DTU
      DT20 = 2.0D0*DTU
*
*       Check for stabilization of binaries (skip H > 0 or GAMMA > GMAX).
      IF (H(IPAIR).LT.0.0.AND.GAMMA(IPAIR).LT.GMAX) THEN
          A2 = 2.0D0*(UDOT(1,IPAIR)**2 + UDOT(2,IPAIR)**2 +
     &                UDOT(3,IPAIR)**2 + UDOT(4,IPAIR)**2) -
     &                BODY(I) - H(IPAIR)*R(IPAIR)
*       Include the stabilization term in the predicted force only.
          A3 = 0.2D0*A2*BODYIN/DTU
      ELSE
          A3 = 0.0D0
      END IF
*
*       Predict U, UDOT & R to order FUDOT3 and H to order HDOT4.
      DO 30 K = 1,4
          FSTAB = FU(K,IPAIR) - A3*UDOT(K,IPAIR)
          F2DOTK = D3U(K,IPAIR)*T12PR + D2U(K,IPAIR)
          UI(K) = ((((D3U(K,IPAIR)*DT06 + F2DOTK)*DT12 +
     &                                FUDOT(K,IPAIR))*DTU + FSTAB)*DTU +
     &                                  UDOT(K,IPAIR))*DTU + U0(K,IPAIR)
          UIDOT(K) = (((D3U(K,IPAIR)*DT34 + F2DOTK)*DT19 +
     &                FUDOT(K,IPAIR))*DT32 + FSTAB)*DT20 + UDOT(K,IPAIR)
   30 CONTINUE
*
      R(IPAIR) = UI(1)**2 + UI(2)**2 + UI(3)**2 + UI(4)**2
*
*       Predict H from Taylor series derivatives without factorials.
      HDOT2 = (D3HDOT(IPAIR)*T2PR + D2HDOT(IPAIR))*T1PR + D1HDOT(IPAIR)
      HDOT3 = D3HDOT(IPAIR)*T12PR + D2HDOT(IPAIR)
      H(IPAIR) = (((D3HDOT(IPAIR)*DT34 + HDOT3)*ONE3*DTU +
     &                    0.5D0*HDOT2)*DTU + HDOT(IPAIR))*DTU + H(IPAIR)
*
*       Form relative coordinates obtained from explicit KS transformation.
      Q1 = UI(1)**2 - UI(2)**2 - UI(3)**2 + UI(4)**2
      Q2 = UI(1)*UI(2) - UI(3)*UI(4)
      Q3 = UI(1)*UI(3) + UI(2)*UI(4)
      Q2 = Q2 + Q2
      Q3 = Q3 + Q3
*
*       Assign global coordinates of regularized components.
      A2 = BODY(I1+1)*BODYIN
      XI(1) = X(1,I) + A2*Q1
      XI(2) = X(2,I) + A2*Q2
      XI(3) = X(3,I) + A2*Q3
      XI(4) = XI(1) - Q1
      XI(5) = XI(2) - Q2
      XI(6) = XI(3) - Q3
*
      RETURN
*
      END
