      SUBROUTINE KSCORR(IPAIR,UI,UIDOT,FP,FD,TDOT,TIME0)
*
*
*       Corrector for KS motion.
*       -----------------------
*
      INCLUDE 'common6.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
      PARAMETER (ONE30=1.0/30.0D0, ONE42=1.0/42.0D0, ONE56=1.0/56.0D0,
     &           ONE72=1.0/72.0D0, ONE90=1.0/90.0D0, ONE110=1.0/110.0D0)
      REAL*8  UI(4),UIDOT(4),FP(6),FD(6),TDOT(11)
      REAL*8  A1(3,4),A(8),U2(4),U3(4),U4(4),U5(4),FP1(4),FD1(4)
      REAL*8  F0(4),F1(4),F2(4),F3(4)
      REAL*8  RDOT(3)
*
*
*       Convert from physical to regularized derivative using T' = R.
      RI = UI(1)**2 + UI(2)**2 + UI(3)**2 + UI(4)**2
      DO 1 K = 1,3
          FD(K) = RI*FD(K)
    1 CONTINUE
*
*       Include KS slow-down factor in the perturbation if ZMOD > 1.
      IF (KZ(26).GT.0) THEN
          IMOD = KSLOW(IPAIR)
          IF (IMOD.GT.1) THEN
              ZMOD = FLOAT(ISLOW(IMOD))
              DO 5 K = 1,3
                  FP(K) = ZMOD*FP(K)
                  FD(K) = ZMOD*FD(K)
    5         CONTINUE
          END IF
      END IF
*
*       Predict H to order HDOT4.
      DTU = DTAU(IPAIR)
      DT1 = 1.0/DTU
      DT2 = DT1**2
      DT12 = 12.0*DT2
*
*       Activate perturbation evaluation after one iteration above 1D-04.
      ITP = 0
*     IF (GAMMA(IPAIR).GT.1.0D-04) ITP = 1      ! Suppressed (Ryabov 5/18).
*
      W2 = 0.5*H(IPAIR) ! -w^2
      Z = W2*DTU**2     ! -(w*t)^2
      C5 = 0.05D0*(1.0 + ONE42*Z*(1.0 + ONE72*Z*(1.0 + ONE110*Z)))
      C4 =  ONE12*(1.0 + ONE30*Z*(1.0 + ONE56*Z*(1.0 + ONE90*Z)))
      C3 = ONE6 * (1.0 + Z*C5)
      C2 = 0.5 * (1.0 + Z*C4)
      C1 = 1.0 + Z*C3
      C0 = 1.0 + Z*C2
*
*       Perform one iteration without re-evaluating perturbations.
      DO 30 ITER = 1,2
*
*       Obtain new transformation matrix.
          CALL MATRIX(UI,A1)
*
*       Form twice regularized force and half first derivative of H.
          HD = 0.0D0
          TD2 = 0.0D0
          DO 10 K = 1,4
              A(K) = A1(1,K)*FP(1) + A1(2,K)*FP(2) + A1(3,K)*FP(3)
              A(K+4) = A1(1,K)*FD(1) + A1(2,K)*FD(2) + A1(3,K)*FD(3)
              HD = HD + UIDOT(K)*A(K)
              TD2 = TD2 + UI(K)*UIDOT(K)
   10     CONTINUE
*
*       Set regularized velocity matrix (Levi-Civita matrix not required).
          CALL MATRIX(UIDOT,A1)
*
*       Include the whole (L*F)' term in explicit derivatives of FU & H.
          DO 20 K = 1,4
              A(K+4) = A(K+4) + A1(1,K)*FP(1) + A1(2,K)*FP(2) +
     &                                                     A1(3,K)*FP(3)
              FP1(K) = 0.5*RI*A(K)
              FD1(K) = TD2*A(K) + 0.5*RI*A(K+4) + HD*UI(K)
              F0(K) = FP0(K,IPAIR)
              F1(K) = FD0(K,IPAIR)
              F3(K) = (F1(K) + FD1(K) - 2.0*(FP1(K)-F0(K))*DT1)*DT2
              F2(K) = 0.5*(FD1(K)-F1(K))*DT1 - 1.5*DTU*F3(K)
*       Save corrected U and UDOT in local variables.
              UI(K) = C0*U0(K,IPAIR) + DTU*(C1*UDOT(K,IPAIR) +
     &                                 DTU*(C2*F0(K) + DTU*(C3*F1(K) +
     &                                 DTU*(C4*F2(K) + C5*F3(K)*DTU))))
              UIDOT(K) = W2*DTU*C1*U0(K,IPAIR) + C0*UDOT(K,IPAIR) +
     &                          DTU*(C1*F0(K) + DTU*(C2*F1(K) +
     &                          DTU*(2.0*C3*F2(K) + 3.0*C4*F3(K)*DTU)))
   20     CONTINUE
*
*       Set improved separation.
          RI = UI(1)**2 + UI(2)**2 + UI(3)**2 + UI(4)**2
*       Check for extra perturbation calculation.
          IF (ITER.EQ.ITP) THEN
*       Note only basic KS transformation is needed here.
              CALL KSTRAN2(UI,UIDOT,Q1,Q2,Q3,RDOT)
              I = N + IPAIR
              BODYIN = 1.0/BODY(I)
              I1 = 2*IPAIR - 1
              NNB0 = LIST(1,I1)
*       Reverse index to distinguish KSCORR call (no PN energy correction).
              I1 = -I1
*       Adopt KSINT procedure based on local variables for next iteration.
              CALL KSPERT2(I1,I,NNB0,BODYIN,Q1,Q2,Q3,RDOT,FP,FD,TIME0)
              DO 28 K = 1,3
                  FD(K) = RI*FD(K)
   28         CONTINUE
          END IF
   30 CONTINUE
*
*       Form U2,...,U5 and copy corrected values to common.
      HD = 0.0D0
      HD2 = 0.0D0
      DO 35 K = 1,4
          B0 = W2*(C0*U0(K,IPAIR) + DTU*C1*UDOT(K,IPAIR)) +
     &                                           C0*F0(K) + DTU*C1*F1(K)
          B1 = W2*(W2*DTU*C1*U0(K,IPAIR) + C0*UDOT(K,IPAIR) +
     &                       DTU*C1*F0(K)) + C0*F1(K) + 2.0*DTU*C1*F2(K)
          U2(K) = B0 + 2.0*DTU**2*(C2*F2(K) + 3.0*C3*F3(K)*DTU)
          U3(K) = B1 + 6.0*DTU**2*C2*F3(K)
          U4(K) = W2*B0 + 2.0*(C0*F2(K) + 3.0*C1*F3(K)*DTU)
          U5(K) = W2*B1 + 6.0*C0*F3(K)
          U(K,IPAIR) = UI(K)
          U0(K,IPAIR) = UI(K)
          UDOT(K,IPAIR) = UIDOT(K)
          FU(K,IPAIR) = 0.5*U2(K)
          FUDOT(K,IPAIR) = ONE6*U3(K)
          FUDOT2(K,IPAIR) = U4(K)
          FUDOT3(K,IPAIR) = U5(K)
          FP0(K,IPAIR) = FP1(K)
          FD0(K,IPAIR) = FD1(K)
*       Note: A(K) correspond to previous UI/UIDOT estimation
          HD = HD + UIDOT(K)*A(K)
          HD2 = HD2 + U2(K)*A(K) + UIDOT(K)*A(K+4)
   35 CONTINUE
      HD = 2.0D0*HD
      HD2 = 2.0D0*HD2
      HD3 = (HD2 - HDOT2(IPAIR))/DTU
      HD4 = (0.5*(HD2 + HDOT2(IPAIR)) + (HDOT(IPAIR) -HD)/DTU)*DT12
      HPRED = H(IPAIR) + (0.5*(HD + HDOT(IPAIR)) +
     &                   ONE12*(HDOT2(IPAIR) - HD2)*DTU)*DTU
*
*       Save new energy and four derivatives.
      H(IPAIR) = HPRED
      HDOT(IPAIR) = HD
      HDOT2(IPAIR) = HD2
      HDOT3(IPAIR) = HD3 + 0.5*HD4*DTU
      HDOT4(IPAIR) = HD4
*
*       Form scalar terms for time derivatives (1/2 of TDOT4 & 1/4 of TDOT5).
      TD2 = 0.0D0
      TD3 = 0.0D0
      TDOT4 = 0.0D0
      TDOT5 = 0.0D0
      TDOT6 = 0.0D0
      TDOT7 = 0.0D0
      TDOT8 = 0.0D0
      TDOT9 = 0.0D0
      TDOT10 = 0.0D0
      TDOT11 = 0.0D0
      DO 40 K = 1,4
          TD2 = TD2 + UI(K)*UIDOT(K)
          TD3 = TD3 + UIDOT(K)**2 + UI(K)*U2(K)
          TDOT4 = TDOT4 + UI(K)*U3(K) + 3.0D0*UIDOT(K)*U2(K)
          TDOT5 = TDOT5 + U4(K)*UI(K) +
     &                             4.0D0*U3(K)*UIDOT(K) + 3.0D0*U2(K)**2
          TDOT6 = TDOT6 + U5(K)*UI(K) +
     &                             5.0*U4(K)*UIDOT(K) + 10.0*U3(K)*U2(K)
          TDOT7 = TDOT7 + 12.0*U5(K)*UIDOT(K) +
     &                                  30.0*U4(K)*U2(K) + 20.0*U3(K)**2
          TDOT8 = TDOT8 + 42.0*U5(K)*U2(K) + 70.0*U3(K)*U4(K)
          TDOT9 = TDOT9 + 112.0*U5(K)*U3(K) + 70.0*U4(K)**2
          TDOT10 = TDOT10 + U5(K)*U4(K)
          TDOT11 = TDOT11 + U5(K)**2
   40 CONTINUE
*
*       Save R and second & third time derivatives.
      R(IPAIR) = RI
      TDOT2(IPAIR) = 2.0D0*TD2
      TDOT3(IPAIR) = 2.0D0*TD3
*
      TDOT(1)  = RI
      TDOT(2)  = TDOT2(IPAIR)
      TDOT(3)  = TDOT3(IPAIR)
      TDOT(4)  = 2.0*TDOT4
      TDOT(5)  = 2.0*TDOT5
      TDOT(6)  = 2.0*TDOT6
      TDOT(7)  = TDOT7
      TDOT(8)  = TDOT8
      TDOT(9)  = TDOT9
      TDOT(10) = 252.0*TDOT10
      TDOT(11) = 252.0*TDOT11
*
      RETURN
*
      END
