      SUBROUTINE CMBODY(ENERGY,NSYS)
*
*
*       Formation of c.m. body by collision.
*       ------------------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX4=4*NMX)
      COMMON/CLOSE/  RIJ4(4,4),RCOLL4,QPERI4,SIZE4(4),ECOLL4,IP(4)
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIK(NMX,NMX),SIZE(NMX),ASYNC,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/EBSAVE/  EBS
      REAL*8  CM(6),LUMS(10),TSCLS(20),GB(10)
      REAL*8  M01,M02,M03,M1,M2,M3,LUM,MC
      CHARACTER*8  WHICH1
*
*
*       Distinguish between chain and triple or quad case (ICH > 0 or = 0).
      IF (IPHASE.EQ.9) THEN
          ICH = 1
      ELSE
*       Activate collision indicator (otherwise done in CHTERM).
          ICH = 0
          IPHASE = 9
      END IF
*
*       Specify global indices of subsystem (membership: NSYS = 2 - 5).
      IF (NSYS.EQ.2) THEN
*
*       Define discrete time for prediction & new polynomials (T <= TBLOCK).
          TIME2 = TIME - TPREV
          DT = 0.1*STEP(N+KSPAIR)
          CALL STEPK(DT,DTN)
          TIME = TPREV + INT((TIME2 + DT)/DTN)*DTN
          TIME = MIN(TBLOCK,TIME)
          TIME0 = TIME
*
*       Save binding energy and terminate KS pair.
          EB = BODY(2*KSPAIR-1)*BODY(2*KSPAIR)*H(KSPAIR)/BODY(N+KSPAIR)
          I = N + KSPAIR
          JCLOSE = 0
*
*       Check for hierarchical configuration.
          NP1 = LIST(1,2*KSPAIR-1) + 1
          DO 5 L = 2,NP1
              J = LIST(L,2*KSPAIR-1)
              RIJ2 = 0.0
              VIJ2 = 0.0
              RDOT = 0.0
              DO 2 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
                  VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,J))**2
                  RDOT = (X(K,I) - X(K,J))*(XDOT(K,I) - XDOT(K,J))
    2         CONTINUE
              RIP = SQRT(RIJ2)
              A1 = 2.0/RIP - VIJ2/(BODY(I) + BODY(J))
              A1 = 1.0/A1
              IF (1.0/A1.GT.0.5/RMIN) THEN
                  ECC2 = (1.0 - RIP/A1)**2 +
     &                                  RDOT**2/(A1*(BODY(I) + BODY(J)))
                  RP = A1*(1.0 - SQRT(ECC2))
                  A0 = -0.5*BODY(I)/H(KSPAIR)
                  ECC = 1.0 - R(KSPAIR)/A0
                  RA = A0*(1.0 + ECC)
                  SR = RP/RA
                  WRITE (6,4)  KSPAIR, NAME(J), H(KSPAIR), ECC, A0, A1,
     &                         RP, ECC, SR
    4             FORMAT (' HIERARCHY:   KS NMJ H E A0 A1 RP E1 SR',
     &                              2I6,F7.0,F8.4,1P,3E9.1,0PF6.2,F6.1)
              END IF
*       Select closest single body inside 0.5*RMIN as KS component.
              IF (RIP.LT.0.5*RMIN.AND.J.LE.N) THEN
                  IF (JCLOSE.GT.0) THEN
                      IF (RIP.GT.RIP0) GO TO 5
                      JCLOSE = J
                      RIP0 = RIP
                  ELSE
                      JCLOSE = J
                      RIP0 = RIP
                  END IF
              END IF
    5     CONTINUE
*
*       Search for evidence of recent regularization.
          NAM1 = NAME(2*KSPAIR-1)
          NAM2 = NAME(2*KSPAIR)
          NNB = LISTD(1)
          DO 7 K = 2,NNB+1
              IF (LISTD(K).EQ.NAM1.OR.LISTD(K).EQ.NAM2) THEN
                  WRITE (6,6)  NAM1, NAM2, LISTD(K), K
    6             FORMAT (' KS REMNANT:    NAM LISTD K  ',3I6,I4)
              END IF
    7     CONTINUE
*
*       Update body #JCLOSE to current time for new KS with combined c.m.
          IF (JCLOSE.GT.0) THEN
              CALL XVPRED(JCLOSE,-1)
              T0(JCLOSE) = TIME
              DO 8 K = 1,3
                  X0DOT(K,JCLOSE) = XDOT(K,JCLOSE)
                  X0(K,JCLOSE) = X(K,JCLOSE)
    8         CONTINUE
          END IF
*
*       Ensure orbit is at pericentre (perturbed hyperbolic case is OK).
          SEMI = -0.5*BODY(N+KSPAIR)/H(KSPAIR)
          IF (R(KSPAIR).GT.SEMI.AND.SEMI.GT.0.0) THEN
              CALL KSAPO(KSPAIR)
              CALL KSPERI(KSPAIR)
*       Restore quantized time to avoid small STEP (KSTERM needs T0 = TIME).
              TIME = TIME0
              T0(2*KSPAIR-1) = TIME
          END IF
*
*       Save collision distance and VINF (km/sec).
          RCOLL = R(KSPAIR)
          VINF = 0.0
          IF (H(KSPAIR).GT.0.0) VINF = SQRT(2.0*H(KSPAIR))*VSTAR
          ECC = 1.0 - R(KSPAIR)/SEMI
*
*       Terminate KS pair and set relevant indices for collision treatment.
          CALL KSTERM
          I1 = 2*NPAIRS + 1
          I2 = I1 + 1
          I3 = 0
          ICOMP = I1
          DMIN2 = MIN(DMIN2,RCOLL)
          WHICH1 = ' BINARY '
      ELSE
          I1 = JLIST(1)
          I2 = JLIST(2)
          I3 = JLIST(3)
          I4 = JLIST(4)
          I5 = JLIST(5)
*       Ignore case of three-body system here (JLIST(4) = 0).
          VINF = 0.0
          ECC = 1.0 + 2.0*EBS*DMINC/(BODY(I1)*BODY(I2))
          ECC = MAX(ECC,0.0D0)
          IF (EBS.GT.0.0) THEN
              HI = EBS*(BODY(I1) + BODY(I2))/(BODY(I1)*BODY(I2))
              VINF = SQRT(2.0*HI)*VSTAR
          END IF
      END IF
*
*       Set new quantized time (note: restore if problems in CHTERM).
*     TIME = TBLOCK
*
*       Form global c.m. coordinates & velocities from body #I1 & I2.
      ZM = BODY(I1) + BODY(I2)
      DO 10 K = 1,3
          CM(K) = (BODY(I1)*X(K,I1) + BODY(I2)*X(K,I2))/ZM
          CM(K+3) = (BODY(I1)*XDOT(K,I1) + BODY(I2)*XDOT(K,I2))/ZM
   10 CONTINUE
*
*       Change stellar evolution parameters if mass loss option is active.
      IF (KZ(19).GE.3) THEN
          M1 = BODY(I1)*ZMBAR
          M2 = BODY(I2)*ZMBAR
          RS1 = RADIUS(I1)*SU
          RS2 = RADIUS(I2)*SU
          TPHYS = (TIME + TOFF)*TSTAR
*
*       Obtain new stellar parameters based on complete mixing.
          IF (MAX(KSTAR(I1),KSTAR(I2)).LE.2) THEN
*
*       Determine evolution time scales for first star.
              M01 = BODY0(I1)*ZMBAR
              KW = KSTAR(I1)
              CALL STAR(KW,M01,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
              AGE1 = TPHYS - EPOCH(I1)
              TMS1 = TM
*
*       Obtain time scales for second star.
              M02 = BODY0(I2)*ZMBAR
              KW = KSTAR(I2)
              CALL STAR(KW,M02,M2,TM,TN,TSCLS,LUMS,GB,ZPARS)
              AGE2 = TPHYS - EPOCH(I2)
              TMS2 = TM
*
*       Specify new age based on complete mixing.
              M3 = M1 + M2
              M03 = M01 + M02
              KW = 1
*       Determine consistent stellar type.
              CALL STAR(KW,M03,M3,TM,TN,TSCLS,LUMS,GB,ZPARS)
              TMS3 = TM
              AGE3 = 0.1*TMS3*(AGE1*M1/TMS1 + AGE2*M2/TMS2)/M3
              CALL HRDIAG(M03,AGE3,M3,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                    RM,LUM,KW,MC,RCC)
*       Update initial mass, epoch & evolution times.
              BODY0(I1) = M03/ZMBAR
              EPOCH(I1) = TPHYS - AGE3
              TEV(I1) = TIME
              TEV0(I1) = TEV(I1)
              KSTAR(I1) = KW
              RADIUS(I1) = RM/SU
              IF (KW.LE.1.AND.TM.LT.TPHYS) THEN
                  NBS = NBS + 1
              END IF
          ELSE
*       Re-define new star with zero age on main sequence for simplicity.
              KSTAR(I1) = 1
              BODY0(I1) = ZM
              EPOCH(I1) = TPHYS
              TEV(I1) = TIME + 0.01/TSTAR
              TEV0(I1) = TEV(I1)
              RADIUS(I1) = MAX(RADIUS(I1),RADIUS(I2))
          END IF
*
          WRITE (6,15)  NAME(I1), KSTAR(I1), TEV(I1)*TSTAR, M1, M2,
     &                  RS1, RS2, RADIUS(I1)*SU
   15     FORMAT (' NEW STAR:    NAME KW TEV M1 M2 R1 R2 R* ',
     &                           I6,I4,F7.1,2F6.1,F7.1,2F6.1)
          TEV(I2) = 1.0E+10
      END IF
*
*       Copy all members of neighbour or perturber list.
      NNB = LIST(1,I1)
      DO 20 L = 1,NNB
          JPERT(L) = LIST(L+1,I1)
   20 CONTINUE
*
*       Evaluate potential energy with respect to colliding bodies.
      IF (NSYS.EQ.2) THEN
          JLIST(1) = I1
          JLIST(2) = I2
*       Replace second old KS component temporarily by arbitrary body.
          JPERT(1) = N
          IF (I2.EQ.N) JPERT(1) = N + 1
          CALL NBPOT(2,NNB,POT1)
      ELSE
*       Obtain differential effect on #I1 & #I2 due to closest member