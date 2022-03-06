      SUBROUTINE KSCORR(IPAIR,UI,UIDOT,FP,TD2,TDOT4,TDOT5)
*
*
*       Corrector for KS motion.
*       -----------------------
*
      INCLUDE 'common6.h'
      COMMON/SLOW/  RANGE,ISLOW(10)
      REAL*8  UI(4),UIDOT(4),FP(6),FREG(4),A1(3,4)
*
*
*       Set new transformation matrix.
      A1(1,1) =  UI(1)
      A1(1,2) = -UI(2)
      A1(1,3) = -UI(3)
      A1(1,4) =  UI(4)
      A1(2,1) =  UI(2)
      A1(2,2) =  UI(1)
      A1(2,3) = -UI(4)
      A1(2,4) = -UI(3)
      A1(3,1) =  UI(3)
      A1(3,2) =  UI(4)
      A1(3,3) =  UI(1)
      A1(3,4) =  UI(2)
*
*       Set time intervals for new divided differences and update the times.
      DTU = DTAU(IPAIR)
      DTUIN = 1.0/DTU
      DT06 = 0.6D0*DTU
      DT12 = ONE12*DTU
      T1PR = T0U(IPAIR) - T1U(IPAIR)
      T2PR = T0U(IPAIR) - T2U(IPAIR)
      T12PR = T1PR + T2PR
      DT1 = TAU - T1U(IPAIR)
      DT2 = TAU - T2U(IPAIR)
      DT3 = TAU - T3U(IPAIR)
      T3PR = T0U(IPAIR) - T3U(IPAIR)
      S2 = T1PR*T2PR
      S3 = S2*T3PR
      S4 = S2 + T3PR*T12PR
      S5 = T12PR + T3PR
      S6 = (((0.6666666666667D0*DTU + S5)*DT06 + S4)*DT12 + ONE6*S3)*DTU
      S7 = ((0.2D0*DTU + 0.25D0*S5)*DTU + ONE3*S4)*DTU + 0.5D0*S3
*
      T3U(IPAIR) = T2U(IPAIR)
      T2U(IPAIR) = T1U(IPAIR)
      T1U(IPAIR) = T0U(IPAIR)
      T0U(IPAIR) = TAU
      A2 = 1.0/DT1
      A3 = 1.0/DT2
      A4 = DTU*DTU/DT3
      TEMP1 = 0.5D0*H(IPAIR)
      TEMP2 = 0.5D0*R(IPAIR)
      HD = 0.0D0
*
*       Include KS slow-down factor in the perturbation if ZMOD > 1.
      IF (KZ(26).GT.0) THEN
          IMOD = KSLOW(IPAIR)
          IF (IMOD.GT.1) THEN
              ZMOD = FLOAT(ISLOW(IMOD))
              DO 5 K = 1,3
                  FP(K) = ZMOD*FP(K)
    5         CONTINUE
          END IF
      END IF
*
*       Set new force & differences and include semi-iteration for U & UDOT.
      DO 10 K = 1,4
          PERTK = A1(1,K)*FP(1) + A1(2,K)*FP(2) + A1(3,K)*FP(3)
*       Reversed indices in matrix A1 represent transpose matrix A1T.
          FREG(K) = TEMP1*UI(K) + TEMP2*PERTK
*       New regularized force.
          D1UK = (FREG(K) - 2.0D0*FU(K,IPAIR))*DTUIN
          D2UK = (D1UK - D1U(K,IPAIR))*A2
          D3UK = (D2UK - D2U(K,IPAIR))*A3
          F4DOTK = (D3UK - D3U(K,IPAIR))*A4
          D1U(K,IPAIR) = D1UK
          D2U(K,IPAIR) = D2UK
          D3U(K,IPAIR) = D3UK
          UI(K) = F4DOTK*S6 + UI(K)
          U0(K,IPAIR) = UI(K)
          U(K,IPAIR) = UI(K)
          HD = HD + UIDOT(K)*PERTK
          UIDOT(K) = F4DOTK*S7 + UIDOT(K)
          UDOT(K,IPAIR) = UIDOT(K)
   10 CONTINUE
*
*       Form scalar terms for time derivatives & HDOT and set FU & FUDOT.
      TD2 = 0.0D0
      TD3 = 0.0D0
      S1 = DTU + DT1
      TDOT4 = 0.0D0
      TDOT5 = 0.0D0
*
      DO 20 K = 1,4
          TD2 = TD2 + UI(K)*UIDOT(K)
          TD3 = TD3 + UIDOT(K)**2 + UI(K)*FREG(K)
          FUDOTK = (D3U(K,IPAIR)*DT1 + D2U(K,IPAIR))*DTU + D1U(K,IPAIR)
          F2DOTK = D3U(K,IPAIR)*S1 + D2U(K,IPAIR)
*       One-half the second force derivative.
          TDOT4 = TDOT4 + UI(K)*FUDOTK + 3.0D0*UIDOT(K)*FREG(K)
          TDOT5 = TDOT5 + UI(K)*F2DOTK + 2.0D0*UIDOT(K)*FUDOTK +
     &                                                  1.5D0*FREG(K)**2
*       One-half the fourth derivative and one-quarter the fifth derivative.
          FU(K,IPAIR) = 0.5D0*FREG(K)
          FUDOT(K,IPAIR) = ONE6*FUDOTK
*       Half the regularized force and sixth the first derivative.
   20 CONTINUE
*
*       Set second & third time derivatives (also needed in RESOLV).
      TDOT2(IPAIR) = 2.0D0*TD2
      TDOT3(IPAIR) = 2.0D0*TD3
*
      HD = 2.0D0*HD
      D1HD = (HD - HDOT(IPAIR))*DTUIN
      D2HD = (D1HD - D1HDOT(IPAIR))*A2
      D3HD = (D2HD - D2HDOT(IPAIR))*A3
      HDOT5 = (D3HD - D3HDOT(IPAIR))*A4
*       Set first derivative of binding energy & higher differences.
      HDOT(IPAIR) = HD
      D1HDOT(IPAIR) = D1HD
      D2HDOT(IPAIR) = D2HD
      D3HDOT(IPAIR) = D3HD
*
*       Include semi-iteration for the binding energy and set new distance.
      H(IPAIR) = HDOT5*S7 + H(IPAIR)
      R(IPAIR) = UI(1)**2 + UI(2)**2 + UI(3)**2 + UI(4)**2
*
      RETURN
*
      END
