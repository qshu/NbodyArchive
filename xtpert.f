      SUBROUTINE XTPERT(ACC,XCM,CHTIME)
*
*
*       External perturbations on chain members.
*       ----------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  M,MASS,MC,MIJ,MKK
      PARAMETER  (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &            NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XJ(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX,1)
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      REAL*8  ACC(NMX3),XCM(3),FIRR(3)
*
*
*       Predict current coordinates of perturbers & chain members.
      TIME0 = TIME
      ISUB = ISYS(5)
      TIME = CHTIME + T0S(ISUB)
      CALL XCPRED(1)
      TIME = TIME0
*
*       Consider each chain component in turn.
      NPC = LISTC(1,1) + 1
      DO 10 I = 1,NN
*       Copy the global chain coordinates to scalars.
          XI = XC(1,I)
          YI = XC(2,I)
          ZI = XC(3,I)
*
          DO 2 K = 1,3
              FIRR(K) = 0.0D0
    2     CONTINUE
*
*       Sum perturber contributions over each chain component.
          DO 5 L = 2,NPC
              J = LISTC(L,1)
              A1 = X(1,J) - XI
              A2 = X(2,J) - YI
              A3 = X(3,J) - ZI
              RIJ2 = A1*A1 + A2*A2 + A3*A3
              A6 = BODY(J)/(RIJ2*SQRT(RIJ2))
              FIRR(1) = FIRR(1) + A1*A6
              FIRR(2) = FIRR(2) + A2*A6
              FIRR(3) = FIRR(3) + A3*A6
    5     CONTINUE
*
*       Copy the perturbations.
          IK = 3*(I - 1)
          ACC(IK+1) = FIRR(1)
          ACC(IK+2) = FIRR(2)
          ACC(IK+3) = FIRR(3)
   10 CONTINUE
*
      RETURN
      END
