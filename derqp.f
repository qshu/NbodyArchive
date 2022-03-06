      SUBROUTINE DERQP(Q,CMX,ENERGY,P,CMV,CHTIME,DQ,DX,DE,DP,DV,DT)
*
*       Derivatives of chain variables.
*       -------------------------------
*
      INCLUDE 'COMMON1.CH'
*     INCLUDE 'COMMON2.CH'
      COMMON/CPERT/  RGRAV,GPERT,IPERT
      COMMON/CALLS/  TPR,ICALL,NFN,NREG
      REAL*8  Q(1),CMX(1),P(1),CMV(1)
      REAL*8  DQ(1),DX(3),DP(1),DV(3)
      REAL*8  W(NMX4),AK(NMX4),DK(NMX),FNC(NMX3),FXTNL(NMX3)
      REAL*8  TP(NMX4),TQ(NMX4),UQ(NMX4),FAUX(4),AUX(-1:1),XAUX(3)
      REAL*8  FCM(3)
*
*       Mikkola & Aarseth 1990 eqs. (77) -> (80).
      UC=0.0
      RSUM=0.0
*       (77)
      DO I=1,N-1
      KS = 4*(I-1)
      K1 = KS + 1
      K2 = KS + 2
      K3 = KS + 3
      K4 = KS + 4
*       Obtain W = L?T(Q)P ; L?T = transpose of L times P.
      W(K1)=(Q(K1)*P(K1)-Q(K2)*P(K2)-Q(K3)*P(K3)+Q(K4)*P(K4))
      W(K2)=(Q(K2)*P(K1)+Q(K1)*P(K2)-Q(K4)*P(K3)-Q(K3)*P(K4))
      W(K3)=(Q(K3)*P(K1)+Q(K4)*P(K2)+Q(K1)*P(K3)+Q(K2)*P(K4))
      W(K4)=(Q(K4)*P(K1)-Q(K3)*P(K2)+Q(K2)*P(K3)-Q(K1)*P(K4))
      RIJL=Q(K1)**2+Q(K2)**2+Q(K3)**2+Q(K4)**2
*       Evaluate RSUM for decision-making.
      RSUM=RSUM+RIJL
      RINV(I)=1./RIJL
      A=.5D0*RINV(I)
      UC=UC+MKK(I)*RINV(I)
      DO K=1,4
      W(KS+K)=A*W(KS+K)
      END DO
      END DO
      LRI=N-1
*       (78)
      TKIN=0.0
      DO I=1,N-1
      J1=-1
      J2=+1
      IF(I.EQ.1)J1=0
      IF(I.EQ.N-1)J2=0
      AUX(-1)=0.5D0*TK1(I)
      AUX( 0)=TKK(I)
      AUX(+1)=0.5*TK1(I+1)
      L=4*(I-1)
      DK(I)=0.0
      DO K=1,4
      AA=0.0
      DO J=J1,J2
      LJ=L+4*J
      AA=AA+AUX(J)*W(LJ+K)
      END DO
      AK(L+K)=AA
*       (79)
      DK(I)=DK(I)+AA*W(L+K)
      END DO
*       (80)
      TKIN=TKIN+DK(I)
      END DO
*
*       Obtain physical coordinates.
      DO K=1,3
      XI(K)=0.0
      END DO
*
      DO I=1,N-1
      L=3*(I-1)
      KS=4*(I-1)
      XC(L+1)=Q(KS+1)**2-Q(KS+2)**2-Q(KS+3)**2+Q(KS+4)**2
      XC(L+2)=2.D0*(Q(KS+1)*Q(KS+2)-Q(KS+3)*Q(KS+4))
      XC(L+3)=2.D0*(Q(KS+1)*Q(KS+3)+Q(KS+2)*Q(KS+4))
      DO K=1,3
      XI(L+3+K)=XI(L+K)+XC(L+K)
      END DO
      END DO
*       External force (FXTNL = W', not the 'usual' force!).
      IF (IPERT.GT.0) THEN
*     CALL CHPERT(1)
      CALL XTF(FXTNL,FCM,CMX,CHTIME)
      IF (GPERT.LT.1.0D-06) THEN
      IPERT=0
      ELSE
*     CALL CHPERT(1)
      IPERT=1
      END IF
      END IF
*
*	Non-chained contribution.
      UNC=0.0
      DO I=1,3*(N-1)
      FNC(I)=FXTNL(I)
      END DO
      DO I=1,N-2
      LI=3*(I-1)
      DO J=I+2,N
      LJ=3*(J-1)
      RIJ2=0.0
      IF(J.GT.I+2)THEN
      DO K=1,3
      XAUX(K)=XI(LJ+K)-XI(LI+K)
      RIJ2=RIJ2+XAUX(K)**2
      END DO
      ELSE
      DO K=1,3
      XAUX(K)=XC(LI+K)+XC(LI+K+3)
      RIJ2=RIJ2+XAUX(K)**2
      END DO
      END IF
      RIJ2INV=1./RIJ2
*       Store the inverse distances.
      LRI=LRI+1
      RINV(LRI)=SQRT(RIJ2INV)
      FM=MIJ(I,J)*RINV(LRI)
      UNC=UNC+FM
      FM=FM*RIJ2INV
*       Fij attraction.
      DO K=1,3
      FAUX(K)=-FM*XAUX(K)
      END DO
*       Add the contribution to interactions depending on Rij.
      DO IK=I,J-1
      L=3*(IK-1)
      DO K=1,3
      FNC(L+K)=FNC(L+K)+FAUX(K)
      END DO
      END DO
      END DO
      END DO
*
*	Evaluate UQ & TP.
      DO I=1,N-1
      L1=3*(I-1)+1
      KS=4*(I-1)
      KS1=KS+1
      CALL QFORCE(Q(KS1),FNC(L1),UQ(KS1))
      CALL MATRIX(Q(KS1),AK(KS1),TP(KS1))
*       The * operation (85).
      AK(KS+4)=-AK(KS+4)
      CALL MATRIX(P(KS1),AK(KS1),TQ(KS1))
*
      DO K=1,4
      UQ(KS+K)=UQ(KS+K)-2.0D0*MKK(I)*Q(KS+K)*RINV(I)**2
      TQ(KS+K)=TQ(KS+K)-4.D0*DK(I)*Q(KS+K)
      END DO
      END DO
*	NOTE: The division by R above (in TP & TQ) is delayed.
*
*	Proceed to final evaluation of derivatives (90)->(94).
      UPOT=UC+UNC
      G=1./(TKIN+UPOT)
      H=TKIN-UPOT
      GAMMA=(H-ENERGY)*G
*
      GT= (1.-GAMMA)*G
      GU=-(1.+GAMMA)*G
*
      DO I=1,N-1
      KS=4*(I-1)
C       Apply the division by R here (to TP & TQ).
C       NOTE: TP & TQ never get 'correct' values. In fact TP=R*TPtrue.
      GToverR=GT*RINV(I)
      DO K=1,4
      DQ(KS+K)=GToverR*TP(KS+K)
      DP(KS+K)=-GToverR*TQ(KS+K)-GU*UQ(KS+K)
      END DO
      END DO
      DT=G
      DO K=1,3
      DX(K)=CMV(K)*G
      DV(K)=FCM(K)*G
      END DO
*	Evaluate E'.
      DE=0.0
      DO I=1,N-1
      L1=3*(I-1)+1
      KS=4*(I-1)
      KS1=KS+1
      CALL QFORCE(Q(KS1),FXTNL(L1),FAUX)
      DO K=1,4
      DE=DE+DQ(KS+K)*FAUX(K)
      END DO
      END DO
*
      IF (ICALL.GT.0) THEN
          TPR = G
          ICALL = 0
      END IF
      NFN = NFN + 1
*     IF (NFN.LE.300) THEN
*     WRITE (6,20)  DT,DE,ENERGY,CHTIME
*  20 FORMAT (' DERQP:   DT DE ENERGY CHTIME ',1P,2E10.1,0P,F10.6,F9.5)
*     END IF
*
      RETURN
      END
