      SUBROUTINE KSINTR(I1,KCASE,TIME1)
*
*
*       Regularized integration back in time.
*       TODO: conside what checks to exlude.
*       ------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'kscmn6.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
      REAL*8  UI(4),UIDOT(4),FP(6),FD(6),RIDOT(3),TDOT(11)
*
*
*       Set second component (= 2), pair index (= 1) & c.m. index.
      I2 = I1 + 1
      IPAIR = KVEC(I1)
      I = N + IPAIR
      JPHASE = 0
*
*       Define perturber membership & inverse c.m. mass.
      NNB0 = LIST(1,I1)
      BODYIN = 1.0/BODY(I)
*
*       Perform KS prediction of U & UDOT.
      DTU = DTAU(IPAIR)
      CALL KSPRED(IPAIR,I,BODYIN,DTU,UI,UIDOT,Q1,Q2,Q3,RIDOT)
*
*       Obtain the perturbing force & derivative.
      CALL KSPERT2(I1,I,NNB0,BODYIN,Q1,Q2,Q3,RIDOT,FP,FD,TIME1)
*
*       Save old radial velocity & relative perturbation and set new GAMMA.
      RDOT = TDOT2(IPAIR)
      GI = SQRT(FP(1)**2 + FP(2)**2 + FP(3)**2)*R(IPAIR)**2*BODYIN
      GAMMA(IPAIR) = GI
*
*       Apply the Hermite corrector.
      CALL KSCORR(IPAIR,UI,UIDOT,FP,FD,TDOT,TIME1)
      TD2 = 0.5*TDOT(2)
*
*       Increase regularization time-step counter and update the time.
      NSTEPU = NSTEPU + 1
      T0(I1) = TIME1
*
*       Define useful scalars.
      RI = R(IPAIR)
      HI = H(IPAIR)
      SEMI = -0.5*BODY(I)/HI
*
*       Choose basic regularized step using binding energy.
      W1 = 0.5/ABS(HI)
      W2 = RMIN*BODYIN
      W1 = MIN(W1,W2)
      W2 = SQRT(W1)
*
*       Set new regularized step and convert to physical time units.
      DTU = 4.0*ETAU*W2
      IF (GI.GT.1.0D-02.AND.BODY(I).GT.10.0*BODYM) DTU = 0.5*DTU
*
*       Include convergence criterion DH = H'*DTU + H''*DTU**2/2 = 0.001*|H|.
      IF (GI.GT.1.0D-05) THEN
          DH = 1.0E-03*MAX(ABS(HI),0.1D0)
          XF = 2.0*DH/ABS(HDOT2(IPAIR))
          YF = HDOT(IPAIR)/HDOT2(IPAIR)
          DTU2 = SQRT(XF + YF**2) - ABS(YF)
          DTU = MIN(DTU2,DTU)
      END IF
*
*       Reduce the step for increasing PN near pericentre (GI still small).
      IF (KSTAR(I1) + KSTAR(I2).EQ.28.AND.CLIGHT.GT.0.0) THEN
          IF (RI.LT.0.1*SEMI) THEN
              ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
              DW = 3.0*TWOPI*BODY(I)/(SEMI*CLIGHT**2*(1.0 - ECC2))
*             ECC = SQRT(ECC2)
*             PM = SEMI*(1.0 - ECC)
              IF (DW.GT.1.0D-04) DTU = 0.5*DTU
              IF (DW.GT.5.0D-04) DTU = 0.5*DTU
              IF (DW.GT.1.0D-03) DTU = 0.5*DTU
*             WRITE (6,23)  TIME1, ECC, PM, DTU, DW, SEMI
*  23         FORMAT (' REDUCE DTU    T E PM DTU DW A ',
*    &                                F12.6,F9.5,1P,3E10.2,E14.6)
          END IF
      END IF
*
*       Convert to physical time units.
!     ITER = 0
!  25 STEP(I1) = (((((ZZ*TDOT6*DTU + 0.2D0*TDOT5)*DTU + 0.5D0*TDOT4)*DTU
!    &                     + TDOT3(IPAIR))*ONE6*DTU + TD2)*DTU + RI)*DTU
      S = 0.0
      DO 26 K=11,1,-1
          S = (S + TDOT(K))*DTU/K
   26 END DO
      STEP(I1) = S
      DTAU(IPAIR) = DTU
*
*       See whether the KS slow-down procedure is activated.
      IMOD = KSLOW(IPAIR)
      IF (IMOD.GT.1) THEN
          ZMOD = FLOAT(ISLOW(IMOD))
          STEP(I1) = ZMOD*STEP(I1)
      END IF
*
      RETURN
*
      END
