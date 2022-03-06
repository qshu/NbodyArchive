      SUBROUTINE RESOLVT(IPAIR,TIME1,KCASE)
*
*
*       Transformation of KS variables.
*       -------------------------------
*       KCASE=0 updates X() and XDOT()
*       KCASE=1 updates XKS() and XDOTKS()
*
      INCLUDE 'common6.h'
      INCLUDE 'kscmn6.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
      REAL*8  Q(3),RDOT(3),UI(4),V(4),A1(3,4)
      REAL*8  FF(4),FD(4),FD2(4),FD3(4)
*
*
      I1 = 2*IPAIR - 1
      I2 = 2*IPAIR
*
      IF (LKSINT(IPAIR)) THEN
          NNB = LISTX(1,I1)
          T0I = T0X(I1)
*       Set current values of regularized coordinates & velocities.
          DO 1 K = 1,4
              UI(K) = U0X(K,IPAIR)
              V(K) = UDOTX(K,IPAIR)
    1     CONTINUE
*       Copy binding energy for routine ENERGY.
          HT1 = HX(IPAIR)
      ELSE
          NNB = LIST(1,I1)
          T0I = T0(I1)
          DO 2 K = 1,4
              UI(K) = U0(K,IPAIR)
              V(K) = UDOT(K,IPAIR)
    2     CONTINUE
          HT1 = H(IPAIR)
      END IF
*
*
*       Ensure appropriate prediction (second pair in binary collision).
      IF (NNB.EQ.0.OR.T0I.EQ.TIME1) GO TO 4
*
      IF (LKSINT(IPAIR)) THEN
          RI = RX(IPAIR)
          KSLOWI = KSLOWX(IPAIR)
          TD2 = TDOT2X(IPAIR)
          TD3 = TDOT3X(IPAIR)
          DTI = DTAUX(IPAIR)
          DO 12 K = 1,4
              FF(K) = FUX(K,IPAIR)
              FD(K) = FUDOTX(K,IPAIR)
              FD2(K) = FUDOT2X(K,IPAIR)
              FD3(K) = FUDOT3X(K,IPAIR)
   12     CONTINUE
          HD = HDOTX(IPAIR)
          HD2 = HDOT2X(IPAIR)
          HD3 = HDOT3X(IPAIR)
          HD4 = HDOT4X(IPAIR)
      ELSE
          RI = R(IPAIR)
          KSLOWI = KSLOW(IPAIR)
          TD2 = TDOT2(IPAIR)
          TD3 = TDOT3(IPAIR)
          DTI = DTAU(IPAIR)
          DO 14 K = 1,4
              FF(K) = FU(K,IPAIR)
              FD(K) = FUDOT(K,IPAIR)
              FD2(K) = FUDOT2(K,IPAIR)
              FD3(K) = FUDOT3(K,IPAIR)
   14     CONTINUE
          HD = HDOT(IPAIR)
          HD2 = HDOT2(IPAIR)
          HD3 = HDOT3(IPAIR)
          HD4 = HDOT4(IPAIR)
      END IF
*
*       Predict U & UDOT to order FUDOT3 using third-order time inversion.
      A2 = 1.0/RI
      A3 = A2*(TIME1 - T0I)
*
*       See whether the time interval should be modified by KSLOW procedure.
      IF (KSLOWI.GT.1) A3 = A3/FLOAT(ISLOW(KSLOWI))
*
*       Expand regularized interval to third order.
      A4 = 3.0D0*TD2**2*A2 - TD3
      DTU = ((ONE6*A4*A3 - 0.5D0*TD2)*A2*A3 + 1.0)*A3
*       Apply safety test near small pericentre or for unperturbed motion.
      IF (DTU.GT.DTI) DTU = 0.8*DTI
*
      IF (DTU.GT.0.0D0) THEN
          DTU1 = 0.2D0*DTU
          DTU2 = DTU/24.0D0
          DO 32 K = 1,4
              UI(K) = ((((FD3(K)*DTU1 + FD2(K))*DTU2 + FD(K))*DTU +
     &                                    FF(K))*DTU + V(K))*DTU + UI(K)
              V(K) = (((FD3(K)*DTU2 + ONE6*FD2(K))*DTU +
     &                        3.0D0*FD(K))*DTU + 2.0D0*FF(K))*DTU + V(K)
   32     CONTINUE
*
*       Predict current binding energy per unit mass for routine ADJUST.
          HT1 = (((HD4*DTU2 + ONE6*HD3)*DTU + 0.5D0*HD2)*DTU +
     &                                                     HD)*DTU + HT1
      END IF
*
*       Form current transformation matrix and two-body separation.
    4 CALL MATRIX(UI,A1)
      RI = UI(1)**2 + UI(2)**2 + UI(3)**2 + UI(4)**2
*
*       Obtain relative coordinates & velocities.
      DO 7 J = 1,3
          Q(J) = 0.0D0
          RDOT(J) = 0.0D0
          DO 6 K = 1,4
              Q(J)    = Q(J)    + A1(J,K)*UI(K)
              RDOT(J) = RDOT(J) + A1(J,K)*V(K)
    6     CONTINUE
          RDOT(J) = 2.0D0*RDOT(J)/RI
    7 CONTINUE
*
      I = N + IPAIR
*
*       Set coordinates & velocities of KS components.
      IF (KCASE.EQ.0) THEN
          DO 15 K = 1,3
              X(K,I1) = X(K,I) + BODY(I2)*Q(K)/BODY(I)
              X(K,I2) = X(K,I) - BODY(I1)*Q(K)/BODY(I)
              XDOT(K,I1) = XDOT(K,I) + BODY(I2)*RDOT(K)/BODY(I)
              XDOT(K,I2) = XDOT(K,I) - BODY(I1)*RDOT(K)/BODY(I)
   15     CONTINUE
      ELSE
          DO 16 K = 1,3
              XKS(K,I1) = X(K,I) + BODY(I2)*Q(K)/BODY(I)
              XKS(K,I2) = X(K,I) - BODY(I1)*Q(K)/BODY(I)
              XDOTKS(K,I1) = XDOT(K,I) + BODY(I2)*RDOT(K)/BODY(I)
              XDOTKS(K,I2) = XDOT(K,I) - BODY(I1)*RDOT(K)/BODY(I)
   16     CONTINUE
      END IF
*
      RETURN
*
      END
