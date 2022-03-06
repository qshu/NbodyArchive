      SUBROUTINE KSINT(I1,KCASE)
*
*
*       Regularized integration.
*       ------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'kscmn6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &               BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &               ZP(NTMAX),ES(NTMAX),CZ(2,NTMAX),IOSC(NTMAX),
     &               NAMEC(NTMAX)
      COMMON/SLOW0/  RANGE,ISLOW(10)
      COMMON/CONNECT/  TIME_CH
!     PARAMETER (ZZ=1.0/120.0D0)
      REAL*8  COEFF(11)
      PARAMETER (C3=1.0/3.0D0,C6=1.0/6.0D0,C7=1.0/7.0D0,C9=1.0/9.0D0,
     &           C11=1.0/11.0D0)
      DATA COEFF/1D0, 0.5D0, C3, 0.25D0, 0.2D0, C6, C7, 0.125D0, C9,
     &           0.1D0, C11/
      REAL*8  UI(4),UIDOT(4),XI(6),VI(6),FP(6),FD(6),RIDOT(3)
      REAL*8  EI0(3),EI(3),XX(3,3),VV(3,3),TDOT(11)
      LOGICAL IQ,ICHAIN
      SAVE  EI0
      DATA  EI0 /3*0.0D0/
*
*
      CALL XINIT(I1)
*
*       Set second component (= 2), pair index (= 1) & c.m. index.
      I2 = I1 + 1
      IPAIR = KVEC(I1)
      I = N + IPAIR
      JPHASE = 0
      ITERM = 0
      ICHAIN = .FALSE.
*
*       Define perturber membership & inverse c.m. mass.
      NNB0 = LIST(1,I1)
      BODYIN = 1.0/BODY(I)
*
*       Check for further unperturbed motion.
      IF (NNB0.EQ.0.AND.H(IPAIR).LT.0.0) THEN
          CALL UNPERT(IPAIR,KCASE)
          IF (IPHASE.LT.0) GO TO 100     ! Exit on unperturbed coalescence.
          JPHASE = ISTAT(KCASE)
*
*       Update any unperturbed relativistic KS binary.
          IF (KZ(11).NE.0.AND.LIST(1,I1).EQ.0) THEN
              CALL BRAKE4(I1,I2,KCASE,DW)
              IF (ISTAT(KCASE).LT.0) GO TO 100
          END IF
 
*      Try re-initialize chain WD/BH system after dormant KS (#11 only).
          IF (KZ(11).NE.0.AND.NCH.EQ.0) THEN   !  note absence of LIST(1,I1).
              IF (MIN(KSTAR(I1),KSTAR(I2)).GE.10) THEN
                  SEMI = -0.5*BODY(I)/H(IPAIR)
                  NP = LIST(1,I1)
                  FMAX = 0.0
                  JCL = 0
*      Determine closest perturber based on maximum force (NP may be zero).
                  DO 215 L = 2,NP+1
                      J = LIST(L,I1)
                      RIJ2 = 0.0
                      RD = 0.0
                      DO 210 K = 1,3
                          RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
                          RD = RD +(X(K,I)-X(K,J))*(XDOT(K,I)-XDOT(K,J))
  210                 CONTINUE
                      FJ = BODY(J)/RIJ2
                      IF (FJ.GT.FMAX) THEN
                          FMAX = FJ
                          JCL = J
                          RX2 = RIJ2
                          RDX = RD
                      END IF
  215             CONTINUE
                  JCLOSEX(IPAIR) = JCL
*       Skip outward motion or separation > 5*RMIN but allow NCH = 2.
                  IF (JCL.GT.0) THEN
                      RP = SQRT(RX2)
                      IF (RDX.GE.0.0.OR.RP.GT.5.0*RMIN.OR.
     &                GAMMA(IPAIR).LT.1.0D-05) GO TO 100
                      IF (NAME(JCL).LT.0) GO TO 100
                      WRITE (6,220)  TIME+TOFF, JCL, KSTAR(I1),
     &                               KSTAR(I2), LIST(1,I1),GAMMA(IPAIR),
     &                               SEMI, R(IPAIR)
  220                 FORMAT (' ACTIVATE CHAIN    T JCL K* NP G A R ',
     &                                            F9.1,I7,3I4,1P,3E10.2)
                      KSPAIRX(IPAIR) = IPAIR
                  ELSE IF (SEMI.LT.0.1*RMIN.AND.
     &                                GAMMA(IPAIR).GT.1.0D-03) THEN
                      WRITE (6,222)  TIME+TOFF, SEMI, R(IPAIR),
     &                               GAMMA(IPAIR)
  222                 FORMAT (' ENFORCE CHAIN    T A R G ',
     &                                           F12.4,1P,3E10.2)
                      KSPAIRX(IPAIR) = IPAIR
                  ELSE
                      GO TO 100
                  END IF
*
*       Restore unperturbed motion from BRAKE4 (NP = 0 fixes some problem).
                  IF (GAMMA(IPAIR).LT.1.0D-10) THEN
                      JCLOSEX(IPAIR) = 0
                      JCL = 0
                      LIST(1,I1) = 0
                  END IF
                  KS2 = 0
*       Include case of binary as dominant perturber.
                  IF (JCL.GT.N) THEN
                      KS2 = JCL - N
                      JCOMPX(IPAIR) = JCL
                      JP = JCL - N
                      WRITE (6,223)  KSPAIRX(IPAIR), KS2, JCL, GAMMA(JP)
  223                 FORMAT (' BINARY PERT    KSP KS2 JCLOSE GJP ',
     &                                         2I4,I7,1P,E10.2)
                  ELSE
*       Avoid JCOMP > N & JCLOSE < N for spurious CALL KSTERM in DELAY.
                      JCOMPX(IPAIR) = JCL
                  END IF
                  JPHASE = 8
                  IPHASEX(IPAIR) = 8
                  KS2X(IPAIR) = KS2
                  JCMAXX(IPAIR) = 0
              END IF
          END IF
          GO TO 100
      END IF
*
*       Update perturber list, but keep original time step.
      IF (T0(I1).EQ.T0(I)) THEN
*       Update perturbers list, but keep original time step
          STEPI1 = STEP(I1)
          CALL KSLIST(IPAIR)
          STEP(I1) = STEPI1
*
          IF (NNB0.GT.0.AND.LIST(1,I1).EQ.0) THEN
*       Add a perturber to continue integration.
              LIST(1,I1) = 1
              IF (LIST(1,I).GT.0) THEN
                  LIST(2,I1) = LIST(2,I)
              ELSE
                  LIST(2,I1) = IFIRST
              END IF
          END IF
          NNB0 = LIST(1,I1)
      END IF
*
*       Perform KS prediction of U & UDOT.
      DTU = DTAU(IPAIR)
*       Perform KS prediction of U & UDOT.
      DTU = DTAU(IPAIR)
      CALL KSPRED(IPAIR,I,BODYIN,DTU,UI,UIDOT,Q1,Q2,Q3,RIDOT)
*
*       Obtain the perturbing force & derivative.
      CALL KSPERT2(I1,I,NNB0,BODYIN,Q1,Q2,Q3,RIDOT,FP,FD,TIME)
*
*       Save old radial velocity & relative perturbation and set new GAMMA.
      RDOT = TDOT2(IPAIR)
      GI = SQRT(FP(1)**2 + FP(2)**2 + FP(3)**2)*R(IPAIR)**2*BODYIN
      GAMMA(IPAIR) = GI
*
*       Apply the Hermite corrector.
      CALL KSCORR(IPAIR,UI,UIDOT,FP,FD,TDOT,TIME)
      TD2 = 0.5*TDOT(2)
      IF (IPHASEX(IPAIR).LT.0) GO TO 90     ! KSPERT2 termination for PN.
*
*       Increase regularization time-step counter and update the time.
      NSTEPU = NSTEPU + 1
      T0PREV = T0(I1)
      T0(I1) = TIME
*
*       Define useful scalars.
      RI = R(IPAIR)
      HI = H(IPAIR)
      SEMI = -0.5*BODY(I)/HI
*
*       Initialize termination indicator and check for large perturbation.
      IQ = .FALSE.
      RMASS = 4.0*BODY(I1)*BODY(I2)*BODYIN**2
      IF (GI.LT.0.01*RMASS.OR.NAME(I).LT.0) THEN
          JCOMPX(IPAIR) = 0
          GO TO 20
      END IF
      IF (GI.GT.0.05) GO TO 2
      CALL FLYBY(I,ITERM)
      IF (ITERM.EQ.0.AND.KSTAR(I).LT.0) THEN
          IQ = .TRUE.
      END IF
      IF (ITERM.EQ.1) THEN
*       Delay chain regularization search until near end of block-step.
          IF (TIME + STEP(I1).GT.TBLOCK.AND.KZ(15).GT.0) THEN
              CALL IMPACT(I,JPHASE)
              IF (IPHASEX(IPAIR).GT.0) GO TO 100
          END IF
      ELSE IF (ITERM.EQ.2) THEN
          IQ = .TRUE.
          GO TO 20
      ELSE
          GO TO 20
      END IF
*
*       Find the dominant body for large perturbations.
    2 S = 4.0*STEP(I)
      FMAX = BODY(I)/RI**2
*       Initialize JCOMP for prediction and optional diagnostics in KSTERM.
      JC = 0
      DO 10 L = NNB0+1,2,-1
          J = LIST(L,I1)
*       Only search bodies within twice the c.m. time-step.
          IF (STEP(J).GT.S) GO TO 10
*       Compare strong perturber and either component with current pair.
          DO 5 K = I1,I2
              RIJ2 = (X(1,J) - X(1,K))**2 + (X(2,J) - X(2,K))**2 +
     &                                      (X(3,J) - X(3,K))**2
              IF (BODY(J) + BODY(K).GT.RIJ2*FMAX) THEN
                  JC = J
                  GO TO 11
              END IF
    5     CONTINUE
   10 CONTINUE
   11 JCOMPX(IPAIR) = JC
*
*       Switch to CHAIN regularization on large perturbation.
      IF (GI.GT.1.0.AND.GI.LT.2.0) THEN
          ICHAIN = .TRUE.          ! Allow choice of ARC and termination.
          GO TO 84
      END IF
      IF (GI.GT.0.5) THEN
          IQ = .TRUE.
      END IF
*     IF (HI.GT.0.0.AND.RI.GT.RMIN) GO TO 90
*
*       Try chain regularization if strong perturber forms dominant pair.
      IF (JC.GT.0.OR.GI.GT.0.05) THEN
*       Check optional binary search.
*         IF (KZ(4).GT.0) THEN
*             DGAM = GI - G0
*             K = KZ(4)
*             CALL EVOLVE(IPAIR,K)
*         END IF
*       Avoid termination inside SEMI (loss of energy accuracy in KSINIT).
          IF (IQ.AND.GI.GT.0.05) THEN
              IF (RI.LT.SEMI.OR.TD2.LT.0.0) IQ = .FALSE.
              IF (SEMI.GT.3.0*RMIN.OR.GI.GT.2.0) IQ = .TRUE.
*       Note failed chain with R' < 0 & IQ = .true. leads to KS switching.
          END IF
*       Removal of old GO TO 84 here which misses time-step updating.
      END IF
*
*       Check termination of hyperbolic encounter (R > R0 or R > 2*RMIN).
   20 IF (HI.GT.0.0D0.AND.NAME(I).GT.0) THEN
          IF ((RI.GT.R0(IPAIR).AND.GI.GT.GMAX).OR.RI.GT.2.0*RMIN.OR.
     &        (GI.GT.0.1.AND.TD2.GT.0.0)) THEN
*       Skip termination delay in case of velocity kick (cf. routine KSTERM).
              IF (HI.LT.100.0.OR.GI.GT.0.1.OR.RI.GT.5.0*RMIN) THEN
                  IQ = .TRUE.
              END IF
          END IF
      END IF
*
*       Include hyperbolic termination outside 2*RMIN for small pert.
      IF (HI.GT.0.0.AND.TD2.GT.0.0.AND.(RI.GT.2.0*RMIN.OR.
     &                        (RI.GT.RMIN.AND.GI.GT.0.0001))) GO TO 90
*
*       Choose basic regularized step using binding energy.
      W1 = 0.5/ABS(HI)
      W2 = RMIN*BODYIN
      W1 = MIN(W1,W2)
      W2 = SQRT(W1)
*
*       Set new regularized step and convert to physical time units.
      DTU = 4.0*ETAU*W2
*       Employ safety reduction near eccentric apocentre (4 is too much).
      IF (RI.GT.1.9*SEMI) DTU = 0.5*DTU
!     IF (RI.LT.0.1d0*SEMI) DTU = DTU * (RI/(0.1d0*SEMI))**0.2d0
!     IF (RI.LT.0.1d0*SEMI.AND.TD2.LT.0.0) DTU = 0.5*DTU
!     IF (HI.LT.-1.0) DTU = DTU * (-HI)**(-0.125d0)
      IF (GI.GT.1.0D-02.AND.BODY(I).GT.10.0*BODYM) DTU = 0.5*DTU
!     DTU = MAX(DTU, 1D-10/RI)
*
*       Include convergence criterion DH = H'*DTU + H''*DTU**2/2 = 0.001*|H|.
      IF (GI.GT.1.0D-04) THEN
          DH = 1.0E-03*MAX(ABS(HI),0.1D0)
          XF = 2.0*DH/ABS(HDOT2(IPAIR))
          YF = HDOT(IPAIR)/HDOT2(IPAIR)
          DTU2 = SQRT(XF + YF**2) - ABS(YF)
          DTU = MIN(DTU2,DTU)
      END IF
*
*       Reduce the step for increasing PN near pericentre (GI still small).
      IF (KSTAR(I1) + KSTAR(I2).EQ.27.AND.CLIGHT.GT.0.0) THEN
          IF (RI.LT.0.1*SEMI) THEN
              ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
              DW = 3.0*TWOPI*BODY(I)/(SEMI*CLIGHT**2*(1.0 - ECC2))
*             ECC = SQRT(ECC2)
*             PM = SEMI*(1.0 - ECC)
              IF (DW.GT.1.0D-04) DTU = 0.5*DTU
              IF (DW.GT.5.0D-04) DTU = 0.5*DTU
              IF (DW.GT.1.0D-03) DTU = 0.5*DTU
*             WRITE (6,23)  TIME, ECC, PM, DTU, DW, SEMI
*  23         FORMAT (' REDUCE DTU    T E PM DTU DW A ',
*    &                                F12.6,F9.5,1P,3E10.2,E14.6)
          ELSE IF (SEMI.GT.3.0*RMIN) THEN
*       Include reduction for wide BH binary to catch low-mass perturber.
              DTU = 0.5*DTU
          END IF
      END IF
*
*       Convert DTU to physical time units.
!     ITER = 0
!  25 STEP(I1) = (((((ZZ*TDOT6*DTU + 0.2D0*TDOT5)*DTU + 0.5D0*TDOT4)*DTU
!    &                     + TDOT3(IPAIR))*ONE6*DTU + TD2)*DTU + RI)*DTU
      S = 0.0
      DO 26 K=11,1,-1
          S = (S + TDOT(K))*DTU/K
   26 CONTINUE
      STEP(I1) = S
!*      t(tau) is required to increase monotonically, i.e. dt/dtau > 0
!     DTDTU = ((((0.05D0*TDOT6*DTU + TDOT5)*DTU + 2.0D0*TDOT4)*DTU
!    &              + 3.0D0*TDOT3(IPAIR))*ONE6*DTU + 2.0D0*TD2)*DTU + RI
!     IF (STEP(I1).LT.1.0D-12.OR.DTDTU.LE.0.0) THEN
      IF (STEP(I1).LT.1.0D-12) THEN
          WRITE (6,30)  NAME(I1), KSLOW(IPAIR), HI, RI, DTU, STEP(I1),GI
!         WRITE (6,*) RI/SEMI
!         WRITE (6,*) X(1:3,I1), X(1:3,I2), X(1:3,N+IPAIR)
!         WRITE (6,*) ZZ,TDOT6,TDOT5,TDOT4,TDOT3(IPAIR),TD2,RI,R(IPAIR),
!    &                R0(IPAIR),IPAIR,STEP(I1),DTDTU,0.001*RI
   30     FORMAT (' KSINT SMALL STEP    NM KSL H R DTU S1 G ',
     &                                  I7,I4,1P,5E10.2)
!         DTU = 0.5*DTU
!         ITER = ITER + 1
!         IF (ITER.LT.10) GO TO 25
          STOP
      END IF
      DTAU(IPAIR) = DTU
*
*       See whether the KS slow-down procedure is activated.
      IMOD = KSLOW(IPAIR)
      IF (IMOD.GT.1) THEN
          ZMOD = FLOAT(ISLOW(IMOD))
          STEP(I1) = ZMOD*STEP(I1)
      END IF
*
*       Truncate STEP to quantized value restricted to STEP(I).
      STEP(I1) = MIN(STEP(I1), STEP(I))
      CALL STEPK(STEP(I1),DTN)
*
*       Perform iteration by Brent's method..
      DTU = DTU*DTN/STEP(I1)
      CALL BRENT(TDOT,IPAIR,DTU)
      DTAU(IPAIR) = DTU
      STEP(I1) = DTN
*
*       Check diagnostics print option.
      IF (KZ(10).GE.3) THEN
          WRITE (6,40) IPAIR, TIME, H(IPAIR), RI, DTAU(IPAIR), GI,
     &                  STEP(I1), IMOD, LIST(1,I1)
   40     FORMAT (3X,'KS MOTION',I4,2F10.4,F8.4,1P,3E10.2,0P,2I4)
      END IF
*
*       Employ special termination criterion in merger case.
      IF (NAME(I).LT.0) THEN
*       Set termination flag for perturbed hierarchy (no chain possible).
          IF (GI.GT.0.5) IQ = .TRUE.
*       Terminate if not enough perturbers or GI > 0.01.
          IF (HI.LT.0.0) THEN
*             ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
*             A0 = SEMI*(1.0 + SQRT(ECC2))/RI
*       Replace eccentricity calculation with typical value.
*             A0 = 1.5*SEMI/RI
*             GA = GI*A0*A0*A0
*             IF (GA.GT.0.25.AND.RI.GT.SEMI) IQ = .TRUE.
              IF (RI.GT.30*RMIN.AND.NNB0.GT.0.8*LIST(1,I)) IQ = .TRUE.
*       Delay termination for massive system with small perturbation.
              IF (GI.LT.3.0D-03.AND.BODY(I).GT.20.0*BODYM) IQ = .FALSE.

              IF (GI.GT.0.02.AND.RI.GT.RMIN) IQ = .TRUE.
*             IF (GI.GT.0.01) IQ = .TRUE.
*       Include extra condition for planet case.
              IF (MIN(BODY(I1),BODY(I2)).LT.0.05*BODYM) THEN
                  IF (GI.GT.2.0D-04) IQ = .TRUE.
              END IF
              IF (.NOT.IQ) GO TO 60
*       Include extra condition for planet case.
              IF (MIN(BODY(I1),BODY(I2)).LT.0.05*BODYM) THEN
                  IF (GI.GT.2.0D-04) IQ = .TRUE.
              END IF
          ELSE
              IF (TD2.GT.0.0.AND.(GI.GT.GMAX.OR.RI.GT.RMIN)) IQ = .TRUE.
          END IF
          IF (.NOT.IQ) GO TO 60
      END IF
*
      ICHAIN = .FALSE.
      IF (KZ(24).GT.0.AND.BODY(I).GT.50.0*BODYM) THEN
          IF (SEMI.LT.RMIN.AND.RI.LT.0.1*SEMI) ICHAIN = .TRUE.
          IF (ICHAIN) GO TO 84
      END IF
*
*       Include algorithm for switch to PN treatment in ARC.
      IF (NCH.EQ.0.AND.KSTAR(I1) + KSTAR(I2).EQ.27
     &                          .AND.CLIGHT.GT.0.0) THEN
          ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
          DW = 3.0*TWOPI*BODY(I)/(SEMI*CLIGHT**2*(1.0 - ECC2))
*       Adopt medium PN condition and activate chain at end of block-step.
          IF (DW.GT.1.1D-03.AND.STEP(I1).GT.TBLOCK-TIME) THEN
*       Copy the NEW CHAIN procedure below.
              IPHASEX(IPAIR) = 8
              JPHASE = 8
*             JCLOSEX(IPAIR) = 0
              JCLOSEX(IPAIR) = -1  ! Negative flag ensures NCH = 2 in SETSYS.
              JCOMPX(IPAIR) = 0
              JCMAXX(IPAIR) = 0
              KSPAIRX(IPAIR) = IPAIR
              KS2X(IPAIR) = 0
*       Select closest perturber above GI > 1.0D-03.
              IF (GI.GT.1.0D-03) THEN
                  RX2 = 4.0*RMIN2
                  NNB0 = LIST(1,I1)    ! Note possible change in unpert.f.
                  JX = 0
                  DO 39 L = 2,NNB0+1
                      J = LIST(L,I1)
                      RIJ2 = 0.0
                      DO 38 K = 1,3
                          RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
   38                 CONTINUE
                      IF (RIJ2.LT.RX2) THEN
                          RX2 = RIJ2
                          JX = J
                      END IF
   39             CONTINUE
*       Accept one close encounter inside 2*RMIN.
                  IF (RX2.LT.4.0*RMIN2) THEN
                      JCLOSEX(IPAIR) = JX
                      JCOMPX(IPAIR) = JX
                      JCLOSEX(IPAIR) = -2    ! Consistent with JCLOSE < 0.
                      JCOMPX(IPAIR) = -3
                  END IF
              END IF
              WRITE (6,42)  TOFF+TIME, TBLOCK-TIME, STEP(I1), RI, DW, GI
   42         FORMAT (' ARC SWITCH    T TB-T S1 R DW G',F13.7,1P,5E10.2)
              GO TO 100
          END IF
      END IF
*
*       Delay termination until end of block for large perturbations.
      IF (IQ) THEN
          DTR = TBLOCK - TIME
*         WRITE (6,45)  IPAIR, TTOT, GI, RI, DTR, STEP(I1)
*  45     FORMAT (' TERM TEST    KS T G R DTR DT  ',
*    &                           I4,F10.4,F7.3,1P,E10.2,2E9.1)
*       See whether chain test is needed.
          IF (DTR.LT.STEP(I1)) THEN
              IF (GI.GT.0.02) GO TO 84
              GO TO 100
          END IF
      END IF
      IF (IQ.AND.GI.GT.0.5) GO TO 90
*
*       Check standard termination criterion (suppress on IQ = .true.).
      IF (IQ) GO TO 60
      RRC = 2.0*RMIN
      IF (HI.GE.0.0D0) RRC = MIN(RRC,R0(IPAIR))
      IF (RI.GT.RRC) THEN
*     IF ((RI.GT.R0(IPAIR).OR.RI.GT.2.0*RMIN).AND..NOT.IQ) THEN
*       Include termination for rare tidal capture starting at pericentre.
          IF (KSTAR(I).LT.0.AND.RI.GT.5.0*RMIN) GO TO 90
*       Impose a limit using size of neighbour sphere (100*R > 0.80*RS).
          IF (RI.GT.8.0D-03*RS(I).AND.GI.GT.GMAX) GO TO 90
*       See whether termination can be delayed for sufficient perturbers.
          IF (BODY(I)*SMU.GT.20.0.AND.GI.LT.0.03) GO TO 60
          IF (NNB0.LT.0.80*LIST(1,I).AND.GI.LT.0.1) GO TO 60
*       Check updating of R0 for newly hardened binary orbit.
          IF (HI.LT.-ECLOSE) THEN
              R0(IPAIR) = MAX(RMIN,2.0D0*SEMI)
              R0(IPAIR) = MIN(R0(IPAIR),5.0*RMIN)
              GO TO 70
          END IF
*       Avoid repeated terminations for large distance (other limits above).
          IF (NNB0.GT.2.AND.LIST(1,I).GT.2) GO TO 90
      END IF
      IF (.NOT.IQ.AND.RI.GT.0.1D0*RS(I)) GO TO 90
*
*       End integration cycle for hyperbolic motion.
   60 IF (HI.GE.0.0D0) THEN
          IF (RDOT*TD2.LT.0.0D0) THEN
*       Determine pericentre for hyperbolic two-body motion.
              ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
              QPERI = SEMI*(1.0D0 - SQRT(ECC2))
*             DMIN2 = MIN(DMIN2,QPERI)
*
*       Check optional tidal interaction or stellar collision.
              IF (KZ(19).GE.3.AND.KSTAR(I).LT.10) THEN
                  VINF = SQRT(2.0*HI)*VSTAR
                  KS1 = KSTAR(I1)
                  KS2 = KSTAR(I2)
                  RX1 = MAX(RADIUS(I1),RADIUS(I2))
*       Determine maximum periastron factor for capture (VINF in km/sec).
                  IF (KZ(27).LE.2) THEN
                      RFAC = RPMAX2(RADIUS(I1),RADIUS(I2),BODY(I1),
     &                              BODY(I2),KS1,KS2,VINF)
                      RCAP = RFAC*RX1
                  ELSE
                      DV = SQRT(2.0*HI)
*       Note that Quinlan & Shapiro function returns actual distance.
                      RCAP = RPMAX(BODY(I1),BODY(I2),VSTAR,DV,QPERI)
                  END IF
                  IF (QPERI.LT.5.0*RX1) THEN
                      WRITE (54,54)  TTOT, NAME(I1), NAME(I2), KS1,
     &                              KS2, VINF, RCAP*SU, RX1*SU, QPERI*SU
   54                 FORMAT (' CLOSE   T NAM K* VINF RCAP RX QP  ',
     &                                  F7.1,2I6,2I4,F6.2,3F6.1)
                  END IF
                  IF (QPERI.LT.RCAP) THEN
                      J1 = I1
                      IF (RADIUS(I2).GT.RADIUS(I1)) J1 = I2
                      FAC = 0.5*BODY(I)/BODY(J1)
*       Set possible BH index and check disruption condition (& #43) first.
                      J2 = 2*IPAIR + 1 - J1
                      IF (KZ(43).GE.2.AND.KSTAR(J2).EQ.14) THEN
                          RCOLL = (BODY(J2)/BODY(J1))**0.3333*RADIUS(J1)
                      ELSE IF (KZ(27).LE.2) THEN
*       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
                          RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
                      ELSE IF (CLIGHT.GT.0.0) THEN
                          RCOLL = 6.0*BODY(I)/CLIGHT**2
                      ELSE
                          RCOLL = 0.0
                      END IF
                      WRITE (55,58)  TTOT, IPAIR, NAME(I1), NAME(I2),
     &                               KS1, KS2, KSTAR(I), VINF
   58                 FORMAT (' RPMAX:    T KS NAM K* VINF ',
     &                                    F7.1,I5,2I6,3I4,F7.2)
                      WRITE (55,59)  SQRT(ECC2), HI, R(IPAIR), SEMI,
     &                               QPERI, BODY(I1), BODY(I2),
     &                               BODY(I)*ZMBAR
   59                 FORMAT (' E H R A QP BODY MT ',
     &                          F9.5,1P,6E10.2,0P,F6.1)
                      RI2 = 0.0
                      VI2 = 0.0
                      DO 61 K = 1,3
                          RI2 = RI2 + (X(K,I) - CMR(K))**2
                          VI2 = VI2 + XDOT(K,I)**2
   61                 CONTINUE
                      WRITE (55,62)  SQRT(RI2)/RC, SQRT(VI2)*VSTAR,
     &                               RHOD, RADIUS(I1)*SU, RADIUS(I2)*SU,
     &                               RCAP, RADIUS(J1)/QPERI, RCOLL/QPERI
   62                 FORMAT (' r/RC V* <C> R* RCAP R1/QP RCOLL/QP ',
     &                          2F5.1,3F6.1,3F5.1)
                      CALL FLUSH(55)
                      IF (QPERI.LT.RCOLL) THEN
*       Obtain KS variables at pericentre before merging into one body.
                          CALL KSPERI(IPAIR)
                          KSPAIRX(IPAIR) = IPAIR
                          JPHASE = -1
                          IQCOLL = -2
*                         CALL CMBODY(QPERI,2)
                      ELSE IF (KSTAR(I).GE.0.AND.KZ(27).GT.0) THEN
                          CALL KSTIDE(IPAIR,KCASE,QPERI)
                      END IF
                  END IF
*       Check options for artificial collisions.
              ELSE IF (KZ(27).EQ.-1.AND.KZ(13).LT.0) THEN
                  RFAC = 2.0
                  IF (QPERI.LT.RFAC*MAX(RADIUS(I1),RADIUS(I2))) THEN
                      J1 = I1
                      IF (RADIUS(I2).GT.RADIUS(I1)) J1 = I2
                      FAC = 0.5*BODY(I)/BODY(J1)
*                     IF (KZ(27).LE.2) THEN
*       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
                          RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
*                     ELSE
*                         RCOLL = 6.0*BODY(I)/CLIGHT**2
*                     END IF
                      IF (QPERI.LT.RCOLL) THEN
                          CALL TOUCH(IPAIR,I1,I2,RCOLL)
                      END IF
                  END IF
              END IF
          END IF
          GO TO 100
      END IF
*
*       Check optional output of periapse angle on unit #49.
      IF (KZ(29).GE.0) GO TO 70
*
*       Evaluate Runge-Lenz eccentricity vector from physical variables.
      CALL KSRES(IPAIR,I1,I2,0.0D0)
      VI2 = 0.0
      RD = 0.0
      RR = 0.0
      JC = JCOMPX(IPAIR)
      DO 63 K = 1,3
          XI(K) = X(K,I) - X(K,JC)
          VI(K) = XDOT(K,I) - XDOT(K,JC)
          RR = RR + XI(K)**2
          VI2 = VI2 + VI(K)**2
          RD = RD + XI(K)*VI(K)
   63 CONTINUE
      RR = SQRT(RR)
      EI2 = 0.0
      DO 64 K = 1,3
          EI(K) = (VI2*XI(K) - RD*VI(K))/BODY(NTOT) - XI(K)/RR
          EI2 = EI2 + EI(K)**2
   64 CONTINUE
*       Save initial values for error plot.
      IF (EI0(1).EQ.0.0D0) THEN
          DO 65 K = 1,3
              EI0(K) = EI(K)
   65     CONTINUE
      ELSE
*       Construct cross product to yield deviation (in degrees).
          ERR = (EI(1)*EI0(2) - EI(2)*EI0(1))/EI2
          ERR = 360.0*ASIN(ERR)/TWOPI
          WRITE (49,66)  TIME/TWOPI, ERR
   66     FORMAT (' ',F12.6,1P,E12.4)
      END IF
*
*       Determine new perturbers for binary at apocentre turning point.
   70 IF (RDOT*TD2.GE.0.0D0) THEN
          IF (ICHAIN.AND.RI.GT.SEMI) GO TO 84
          GO TO 100
      END IF
*
*       Check minimum two-body separation just after pericentre.
      IF (RDOT.LT.0.0D0) THEN
*       Obtain pericentre by Mikkola's algorithm (GAMMA < 0.001).
          IF (GI.LT.0.001) THEN
              CALL PERI(UI,UIDOT,RI,BODY(I1),BODY(I2),QPERI)
          ELSE
              QPERI = RI
          END IF
*         DMIN2 = MIN(DMIN2,QPERI)
*
*       Check optional tidal interaction or stellar collision (skip merger).
          IF (KZ(19).GE.3.AND.KSTAR(I).LE.10.AND.NAME(I).GT.0) THEN
              RFAC = 5.0
              IF (KZ(27).LE.2) THEN
                  IF (KZ(27).EQ.1) RFAC = 4.0
                  RX1 = RFAC*MAX(RADIUS(I1),RADIUS(I2))
              ELSE
                  RX1 = RPMIN(BODY(I1),BODY(I2),VSTAR,HI,QPERI)
              END IF
              IF (QPERI.LT.RX1) THEN
                  J1 = I1
                  IF (RADIUS(I2).GT.RADIUS(I1)) J1 = I2
                  FAC = 0.5*BODY(I)/BODY(J1)
*       Set possible BH index and check disruption condition (& #43) first.
                  J2 = 2*IPAIR + 1 - J1
                  IF (KZ(43).GE.2.AND.KSTAR(J2).EQ.14) THEN
                      RCOLL = (BODY(J2)/BODY(J1))**0.3333*RADIUS(J1)
                  ELSE IF (KZ(27).LE.2) THEN
*       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
                      RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
                  ELSE IF (CLIGHT.GT.0.0) THEN
                      RCOLL = 6.0*BODY(I)/CLIGHT**2
                  ELSE
                      RCOLL = 0.0
                  END IF
                  IF (QPERI.LT.RCOLL) THEN
*       Obtain KS variables at pericentre before merging into one body.
                      CALL KSPERI(IPAIR)
                      KSPAIRX(IPAIR) = IPAIR
                      JPHASE = -1
                      IQCOLL = -2
*                     CALL CMBODY(QPERI,2)
                  ELSE IF (KSTAR(I).GE.0) THEN
*       Distinguish between sequential, standard and GR circularization.
                      IF (KZ(27).EQ.1) THEN
                          ICIRC = 1
                          TC = 0.0
                      ELSE IF (KZ(27).EQ.2.AND.KSTAR(I).LT.10) THEN
                          ECC2 = (1.0 - RI/SEMI)**2 +
     &                                    TDOT2(IPAIR)**2/(BODY(I)*SEMI)
                          ECC = SQRT(ECC2)
                          ICIRC = 0
                          CALL TCIRC(QPERI,ECC,I1,I2,ICIRC,TC)
                      ELSE
                          ICIRC = 1
                          TC = 0.0
                      END IF
                      IF (KSTAR(I).GE.10) ICIRC = 0
*       Skip tidal effects for circularization time above 100 Myr (07/08).
                      IF (ICIRC.GT.0.AND.KZ(27).GT.0.AND.
     &                    TC.LT.100.0) THEN
                          CALL KSTIDE(IPAIR,KCASE,QPERI)
                      END IF
                  END IF
              END IF
*       Check for perturbed spiral or chaos case (skip collision).
              IF (KSTAR(I).EQ.-2.AND.IPHASEX(IPAIR).EQ.0) THEN
                  CALL SPIRAL(IPAIR)
              ELSE IF (KSTAR(I).EQ.-1.AND.IPHASEX(IPAIR).EQ.0) THEN
                  CALL KSTIDE(IPAIR,KCASE,QPERI)
              END IF
*       Check options for artificial collisions.
          ELSE IF (KZ(27).EQ.-1.AND.KZ(13).LT.0) THEN
              RFAC = 2.0
              IF (QPERI.LT.RFAC*MAX(RADIUS(I1),RADIUS(I2))) THEN
                  J1 = I1
                  IF (RADIUS(I2).GT.RADIUS(I1)) J1 = I2
                  FAC = 0.5*BODY(I)/BODY(J1)
*       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
                  RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
                  IF (QPERI.LT.RCOLL) THEN
                      CALL TOUCH(IPAIR,I1,I2,RCOLL)
                  END IF
              END IF
          END IF
          GO TO 100
      END IF
*
*       Save maximum separation of persistent binary (skip parallel case).
*     RMAX = MAX(RMAX,RI)
*
*       Check binary reference radius or merger stability criterion.
      IF (NAME(I).GT.0) THEN
*       Update termination length scale in case of initial soft binary.
          EB = BODY(I1)*BODY(I2)*HI*BODYIN
          IF (EB.LT.EBH) R0(IPAIR) = MAX(RMIN,2.0*SEMI)
      ELSE
          ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
          ECC = SQRT(ECC2)
          RP = SEMI*(1.0 - ECC)
*       Find merger index.
          IM = 0
          DO 72 K = 1,NMERGE
              IF (NAMEM(K).EQ.NAME(I)) IM = K
   72     CONTINUE
*       Exclude inner planets from the general stability test.
          IF (MIN(BODY(I1),BODY(I2)).LT.0.05*BODYM) THEN
              IF (RP.LT.R0(IPAIR)) GO TO 90
          END IF
*       Include optional Kozai diagnostics (large EMAX) & circularization.
          IF (KZ(42).GT.0) THEN
              CALL KOZAI(IPAIR,IM,ECC,SEMI,ITERM)
*       Check merger termination and perform circularization or collision.
              IF (ITERM.GT.0) THEN
                  JPHASE = 7
                  KSPAIRX(IPAIR) = IPAIR
                  KSPAIR = IPAIR
                  CALL RESET
                  IF (ITERM.GT.1) CALL CIRC(KCASE)
                  GO TO 100
              END IF
          END IF
*       Assess the stability inside critical pericentre (safe estimate).
          SEMI0 = -0.5*BODY(I1)/HM(IM)
          IF (RP.LT.6.0*SEMI0) THEN
*       Note: assessment needs to use same eccentricity as for acceptance.
              CALL ASSESS(IPAIR,IM,ECC,SEMI,ITERM)
              IF (ITERM.GT.0) THEN
                  INSTAB = INSTAB + 1
                  GO TO 90
              END IF
          END IF
*       Check possible eccentricity modulation or t_circ update.
          IF (IM.GT.0.AND.(TIME.GT.TMDIS(IM).OR.
     &        TMDIS(IM).GT.1.0D+06)) THEN
              IF (KZ(27).EQ.2) THEN
                  CALL ECCMOD(I,ITERM)
                  IF (ITERM.GT.0) THEN
*                     WRITE (6,76)  RP, R0(IPAIR)
*  76                 FORMAT (' ECCMOD TERM    RP R0 ',1P,2E10.2)
                      GO TO 90
                  END IF
*       Consider both inner and possible outer circularizing binary.
                  DO 78 K = 1,NCHAOS
                      IF (NAMEC(K).EQ.NZERO - NAMEM(IM).AND.
     &                    KSTARM(IM).EQ.-2) THEN
*       Update unperturbed binary if T - TOSC > 10 Myr (cf. IMPACT & DECIDE).
                          IF ((TIME - TOSC(K))*TSTAR.GT.10.0) GO TO 90
                      END IF
                      IF (NAMEC(K).EQ.NAMEG(IM).AND.
     &                    KSTARM(IM).EQ.-2) THEN
                          IF ((TIME - TOSC(K))*TSTAR.GT.10.0) GO TO 90
                      END IF
*       Note: perturbed binary is treated if pericentre before next IMPACT.
   78             CONTINUE
              END IF
          END IF
      END IF
*
*       See whether KS slow-down procedure should be (re)-checked (no Chaos).
      IF (KZ(26).GT.0.AND.KSTAR(I).GE.0) THEN
          KMOD = RANGE*GMIN/MAX(GI,1.0D-10)
          IF (KMOD.GT.1.OR.IMOD.GT.1) THEN
              CALL KSMOD(IPAIR,KMOD)
              IF (KMOD.LT.0) GO TO 100
              GO TO 80
          END IF
      END IF
*
*       Set approximate value of next period.
      TK = TWOPI*SEMI*SQRT(SEMI*BODYIN)
      IF (IMOD.GT.1) THEN
          TK = ZMOD*TK
      END IF
*
*       Use old perturber list if next apocentre is before the c.m. step.
      IF (TIME + TK.LT.T0(I) + STEP(I)) THEN
          GO TO 100
      END IF
*
*       Select new perturbers (adopt unperturbed period if none found).
   80 CALL KSLIST(IPAIR)
*
*       Check rectification of chaotic spiral at start of unperturbed motion.
      IF (KSTAR(I).EQ.-2.AND.LIST(1,I1).EQ.0) THEN
          DMR = 0.D0
          CALL CHRECT(IPAIR,DMR)
          IF (IPHASEX(IPAIR).LT.0) JPHASE = -1
          IF (IPHASEX(IPAIR).LT.0) GO TO 100
      ELSE
*         CALL KSRECT(IPAIR)
      CONTINUE
      END IF
*
*       See whether a massive subsystem can be selected for ARC.
   84 IF (NCH.EQ.0.AND.SEMI.LT.5.0*RMIN.AND.NAME(I).GE.0.AND.
     &    (GI.GT.0.02.OR.ICHAIN).AND.
     &    (KZ(30).NE.0.AND.KZ(30).NE.-2)) THEN
*         IF (RI.GT.SEMI) GO TO 88   ! Requires investigation.
*
*       Check optional BH condition (prevents mass-loss complications).
          IF (KZ(11).LE.-2) THEN
              IF (KSTAR(I1).NE.14.OR.KSTAR(I2).NE.14) GO TO 88
          END IF
*
*       Search the perturbers (singles and binaries) inside 10*SEMI.
          RSEP = MAX(SEMI,RMIN)
          RCR2 = 100.0*RSEP**2
          NNB1 = LIST(1,I1) + 1
          JCL = 0
          JCL2 = 0
          FMAX = 0.0
          FMAX2 = 0.0
          JP = 0
          DO 85 L = 2,NNB1
              J = LIST(L,I1)
              IF (J.GT.N) THEN
                  JP = J - N
                  IF (LIST(1,2*JP-1).EQ.0) GO TO 85
              END IF
              RD = (X(1,I)-X(1,J))*(XDOT(1,I)-XDOT(1,J)) +
     &             (X(2,I)-X(2,J))*(XDOT(2,I)-XDOT(2,J)) +
     &             (X(3,I)-X(3,J))*(XDOT(3,I)-XDOT(3,J))
*             RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
*    &                                      (X(3,I) - X(3,J))**2
*     WRITE (6,444) J, NAME(J), SEMI, SQRT(RIJ2), RD, GI
* 444 FORMAT (' SEARCH    J NM A RIJ RD G ',2I6,1P,5E10.2)
*             IF (RD.GE.0.0.AND.GI.LT.0.02) GO TO 85
              RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                      (X(3,I) - X(3,J))**2
*       Save possible strong perturber just outside limit.
              IF (RIJ2.GT.RCR2) THEN
                 IF (RIJ2.LT.4.0*RCR2) THEN
                      JP = J
                      RP2 = RIJ2
                      GO TO 85
                  ELSE
                      GO TO 85
                  END IF
              END IF
              FIJ = (BODY(I) + BODY(J))/RIJ2
              IF (FIJ.GT.FMAX) THEN
                  FMAX2 = FMAX
                  JCL2 = JCL
                  FMAX = FIJ
                  JCL = J
                  RRD = RD
                  RJ2 = RIJ2
                  VIJ2 = (XDOT(1,I) - XDOT(1,J))**2 +
     &                   (XDOT(2,I) - XDOT(2,J))**2 +
     &                   (XDOT(3,I) - XDOT(3,J))**2
              ELSE IF (NAME(J).GE.0.AND.FIJ.GT.FMAX2) THEN
                  FMAX2 = FIJ
                  JCL2 = J
              END IF
   85     CONTINUE
          JCLOSEX(IPAIR) = JCL
          JCMAXX(IPAIR) = JCL2
*
*       Quit on dominant external perturber.
          IF (JCL.GT.0.AND.JP.GT.0) THEN
              JLIST(1) = JP
              CALL FPERT(I,JCL,1,JLIST,PERT)
              GP = PERT*RP2/(BODY(I) + BODY(JCL))
              IF (GP.GT.0.10) GO TO 88
          END IF
*
*       Exclude zero or hierarchical perturber and also large mass ratio.
          IF (JCL.LT.IFIRST) GO TO 88     ! Exclude JCL < IFIRST & = 0.
          IF (NAME(JCL).LE.0) GO TO 88
          IF (BODY(JCL).GT.50.0*BODY(I)) GO TO 88
*       Include delay time to avoid repeat events (DANGER: do not use GPERT).
          IF (TTOT.GE.TIME_CH) THEN
              IF (NAME(I).LT.0) WRITE (6,83)  I, NAME(I), NAME(JCL)
   83         FORMAT (' NEGATIVE ARC    I NM NMJ ',3I7)
              RX1 = SQRT(RJ2)
              SEMI1 = 2.0/RX1 - VIJ2/(BODY(I) + BODY(JCL))
              SEMI1 = 1.0/SEMI1
              ECC2 = (1.0-RX1/SEMI1)**2 +
     &                           RRD**2/((BODY(I) + BODY(JCL))*SEMI1)
              ECC1 = SQRT(ECC2)
              PMIN = SEMI1*(1.0 - ECC1)
*       Increase PMIN for hyperbolic perturber (also JCLOSE > N).
              IF (ECC1.GT.1.0) THEN
                  PMIN = 0.25*RX1      ! Allows more distant perturber.
                  QST = -1.0
              ELSE
*        Apply stability test with actual ECC and inclination.
                  ECC2 = (1.0 - RI/SEMI)**2 + 4.0*TD2**2/(BODY(I)*SEMI)
                  ECC = SQRT(ECC2)
                  DO 660 K = 1,3
                      XX(K,1) = X(K,I1)
                      XX(K,2) = X(K,I2)
                      XX(K,3) = X(K,JCL)
                      VV(K,1) = XDOT(K,I1)
                      VV(K,2) = XDOT(K,I2)
                      VV(K,3) = XDOT(K,JCL)
  660             CONTINUE
                  IF (XDOT(1,I1).EQ.XDOT(1,I2)) THEN
                      CALL RESOLV(IPAIR,1)
                  END IF
                  CALL INCLIN(XX,VV,X(1,I),XDOT(1,I),ANGLE)
                  QST = QSTAB(ECC,ECC1,ANGLE,BODY(I1),BODY(I2),
     &                                                BODY(JCL))
*       Abandon ARC for stable hierarchy but not hyperbolic binary.
                  IF (PMIN.GT.0.9*QST*SEMI.AND.SEMI.GT.0) THEN
*                     WRITE (6,666)  NAME(I1), ECC1, PMIN, QST*SEMI
* 666                 FORMAT (' NOCHAIN!   NM E1 PM PCR ',
*    &                                     I6,F8.3,1P,2E10.2)
                      GO TO 88
                  END IF
              END IF
*       Delay until end of the block-step.
              IF (TBLOCK-TIME.GT.STEP(I1)) GO TO 100
*       Limit energy of triple system (< 5*EBH) using radial velocity.
*             ZMU = BODY(I)*BODY(JCL)/(BODY(I) + BODY(JCL))
*             EBT = EB + ZMU*(0.5*(RD/RX1)**2-(BODY(I)+BODY(JCL)/RX1))
*             IF (EBT.GT.5.0*EBH) GO TO 88
*
              EB = BODY(I1)*BODY(I2)*HI*BODYIN
*             EORB = -0.5*BODY(I)*BODY(JCL)/SEMI1
              WRITE (6,86)  TTOT, NAME(JCL), LIST(1,I1),
     &                      STEP(JCL), SEMI, RX1, EB, GI, QST
   86         FORMAT (' NEW CHAIN   T NMJ NP STEPJ A RIJ EB G QST ',
     &                              F9.3,I7,I4,1P,6E10.2)
*       Set next new chain time to avoid escaper being absorbed.
              TIME_CH = TTOT + 0.001
*       Switch JCLOSE and JCMAX if the latter is a KS binary and JCLOSE < N.
              IF (JCMAXX(IPAIR).GT.N.AND.JCL.LE.N) THEN
                  JSAVE = JCMAXX(IPAIR)
                  JCMAXX(IPAIR) = JCL
                  JCL = JSAVE
                  JCLOSEX(IPAIR) = JCL
              END IF
              IF (JCL.GT.N.AND.JCMAXX(IPAIR).GT.N) THEN
                  WRITE (6,87)  I, JCL, JCMAXX(IPAIR)
   87             FORMAT (' SIX-BODY CHAIN   ICM JCL JCX ',3I7)
              END IF
              JCOMPX(IPAIR) = JCL
              KSPAIRX(IPAIR) = IPAIR
              IPHASEX(IPAIR) = 8
              EBCH0X(IPAIR) = EB
*       Distinguish between case of single and binary intruder(s).
              IF (JCL.LE.N) THEN
                  KS2 = 0
              ELSE
                  KS2 = JCL - N
*       Adopt termination sequence: I (IPAIR), JCLOSE, JCMAX (DELAY, SUBSYS).
              END IF
*       Specify index for ARchain initialization (cf. SUBINT).
              ISTAT(KCASE) = 8
*       Initiate chain regularization directly (B-B or B-S: see IMPACT).
              KS2X(IPAIR) = KS2
              GO TO 100
          ELSE
*       Reset indicators to zero after unsuccessful test.
              JCLOSEX(IPAIR) = 0
              JCMAXX(IPAIR) = 0
          END IF
      END IF
*
*       Include KS termination after failed chain test GO TO 84.
   88 IF (IQ) GO TO 90
*
*       Check optional search criterion for multiple encounter or merger.
      IF (KZ(15).GT.0.AND.STEP(I).LT.DTMIN) THEN
!          write(6,*) 'ksint try impact',NAME(I1),NAME(I2)
          IF (TTOT.GT.TIME_CH) CALL IMPACT(I,JPHASE)
      END IF
      GO TO 100
*
*       Terminate regularization of current pair (IPAIR set in KSPAIR).
   90 KSPAIRX(IPAIR) = IPAIR
*       Set indicator for calling KSTERM in MAIN (permits phase overlay).
      JPHASE = 2
*       Check case of hierarchical binary.
      IF (NAME(I).LT.0) JPHASE = 7
*
*       Save activity index unless already non-zero.
  100 IF (ISTAT(KCASE).EQ.0) ISTAT(KCASE) = JPHASE
*
      CALL XUPDATE(I1)
*
      RETURN
*
      END
