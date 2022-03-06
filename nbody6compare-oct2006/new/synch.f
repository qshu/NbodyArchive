      SUBROUTINE SYNCH(IPAIR)
*
*
*       Spin synchronization of circularized orbit.
*       -------------------------------------------
*
*       Rational function approximations of solution to the Hut
*       evolution equation with spin. Ref: A & A 99, 126, eqn (A15).
*
      INCLUDE 'common6.h'
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &               BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &               RP(NTMAX),ES(NTMAX),CM(2,NTMAX),IOSC(NTMAX),
     &               NAMEC(NTMAX)
      common/spins/angmom0,rg2(2),m21,r21,semi0,C1,C2,C3,C4,semi
      REAL*8  WW(3),QQ(3),W(2),Q(2),AT0(2),M21,WG(2),QG(2),WSCALE(2),
     &        QSCALE(2),A(2),B(2),C(6)
      REAL*8  M0,M1,MC1,MC2,CORERD
      EXTERNAL CORERD
      DATA  WW  /2.119,3.113,8.175/
      DATA  QQ  /0.4909,0.4219,0.2372/
      DATA  A  /6.306505,-7.297806/
      DATA  B  /32.17211,13.01598/
      DATA  C  /5.101417,24.71539,-9.627739,1.733964,
     &                            -2.314374,-4.127795/
      SAVE  ICOUNT,ITRY,ISYNCH,ITER,ITIME,ECCM
      DATA  ICOUNT,ITRY,ISYNCH,ITER,ITIME,ECCM /0,0,0,0,0,0.002D0/
*
*
*       Define c.m. & KS indices and search current names for chaos index.
      I = N + IPAIR
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      IC = 0
      DO 2 K = 1,NCHAOS
          IF (NAMEC(K).EQ.NAME(I)) IC = K
    2 CONTINUE
*
*       Exit if binary is not identified but include updating of small ECC.
      IF (IC.EQ.0) THEN
          JPAIR = -IPAIR
          CALL TRFLOW(JPAIR,DTR)
          M1 = BODY(I1)*SMU
          CALL TRDOT(I1,DTM,M1)
          M1 = BODY(I2)*SMU
          CALL TRDOT(I2,DTM2,M1)
          DTM = MIN(DTM,DTM2)
          DTM = MAX(DTM,0.1D0)
          TEV(I) = TIME + 0.1*MIN(DTM,DTR)
          A0 = -0.5*BODY(I)/H(IPAIR)
          ECC2 = (1.0 - R(IPAIR)/A0)**2 + TDOT2(IPAIR)**2/(A0*BODY(I))
          ECC = SQRT(ECC2)
          ITIME = ITIME + 1
          IF (ITIME.LT.100) THEN
          WRITE (6,3) NAME(I1), NAME(I), KSTAR(I), NCHAOS, ECC, DTM, DTR
    3     FORMAT (' MISSING!   SYNCH    NM K* NCH ECC DTM DTR ',
     &                                  2I6,2I4,F8.4,1P,2E10.2)
          END IF
*
*       Check enforcement of circularization for small TCIRC.
          IF (ECC.LT.0.10.AND.LIST(1,I1).EQ.0) THEN
              ICIRC = -1
              QPERI = A0*(1.0 - ECC)
              CALL TCIRC(QPERI,ECC,I1,I2,ICIRC,TC)
              IF (TC.LT.100.0) THEN
                  CALL DEFORM(IPAIR,ECC,ECCM)
                  ECC = ECCM
              END IF
          END IF
*
*       Initialize chaos elements after possible common envelope or E < 0.1.
          IF (ECC.LE.0.003) THEN
              NCHAOS = NCHAOS + 1
              IC = NCHAOS
              NAMEC(IC) = NAME(I)
              RP(IC) = A0*(1.0 - ECC)
              ES(IC) = ECC
              TOSC(IC) = TIME
              IF (NCHAOS.GT.NTMAX) THEN
                  WRITE (6,4)  NAME(I1), NCHAOS, ECC, A0*(1.0 - ECC)
    4             FORMAT (' FATAL ERROR!    SYNCH    NM NCH E QP ',
     &                                               I6,I4,F8.4,1P,E9.1)
                  STOP
              END IF
      WRITE (6,444)  (NAMEC(K),K=1,NCHAOS)
  444 FORMAT (' SYNCH!!   NMC  ',12I5)
      STOP
          END IF
          GO TO 100
      END IF
*
*       Check Roche overflow time during small intervals.
      IF (TIME - TOSC(IC).LT.1.0) THEN
          JPAIR = -IPAIR
          CALL TRFLOW(JPAIR,DTR)
          IF (DTR.LT.STEP(I)) THEN
              TEV(I) = TIME + STEP(I)
              GO TO 100
          END IF
      END IF
*
*       Copy spiral parameters and set semi-major axis & eccentricity.
      TIME0 = TOSC(IC)
      RP0 = RP(IC)
      ES0 = ES(IC)
      HI = H(IPAIR)
      A0 = -0.5*BODY(I)/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/A0)**2 + TDOT2(IPAIR)**2/(A0*BODY(I))
      ECC = SQRT(ECC2)
      ECC0 = ECC
*
*       Rectify elements after possible common envelope process.
      QPERI = A0*(1.0 - ECC)
      IF (ES0.GT.0.01.OR.ABS(QPERI-RP0).GT.0.01*QPERI) THEN
          DT = TIME - TOSC(IC)
          WRITE (6,5)  NAME(I1), KSTAR(I1),KSTAR(I2),KSTAR(I), ES0,
     &                 ECC, RP0, QPERI, DT
    5     FORMAT (' WARNING!    SYNCH    NM K* ES0 E RP0 QP DT  ',
     &                                   I6,3I4,2F8.4,1P,3E10.2)
          RP0 = QPERI
          ES0 = ECC
      END IF
*
*       Set standard binary or current elements on significant departure.
      IF (ECC.GT.0.01.OR.ES0.GT.0.01) THEN
          ECC = SQRT(ECC2)
          IF (ECC.GT.0.01) THEN
              KSTAR(I) = 0
              WRITE (6,8)  NAME(I1), LIST(1,I1), ES0, ECC, A0*SU
    8         FORMAT (' NON-CIRCULAR    SYNCH    NM NP ES0 ECC A ',
     &                                           I6,I4,2F8.4,1P,E10.2)
              NAMI = -I
              CALL SPIRAL(NAMI)
              GO TO 100
          ELSE
              ES0 = ECC
              RP0 = A0*(1.0 - ECC)
          END IF
      END IF
*
*       Specify index J1 as biggest radius to be used with AT0(1).
      IF (RADIUS(I1).GE.RADIUS(I2)) THEN
          J1 = I1
          J2 = I2
      ELSE
          J1 = I2
          J2 = I1
      END IF
*
*       Define oscillation period (dimensionless time) and damping constants.
      XN = 0.0
      QD = 0.0
      DO 10 K = 1,2
          IK = I1 + K - 1
          IF (K.EQ.1) THEN
              IK = J1
          ELSE
              IK = J2
          END IF
*       Specify polytropic index for each star (n = 3, 2 or 3/2).
          IF (KSTAR(IK).EQ.3.OR.KSTAR(IK).EQ.5.OR.
     &        KSTAR(IK).EQ.6.OR.KSTAR(IK).EQ.9) THEN
              CALL GIANT(IPAIR,IK,WG,QG,WSCALE,QSCALE,XN,QL)
              W(K) = WG(1)
              Q(K) = QG(1)
	      rg2(k)= 0.1*(1.0 - CM(K,IC)/BODY(IK))
              QD = QL
          ELSE
              QL = 1.0E+04
              IP = 3
              IF (KSTAR(IK).GE.3) IP = 2
              IF (KSTAR(IK).EQ.4.OR.KSTAR(IK).EQ.7) IP = 3
              IF (KSTAR(IK).EQ.8) IP = 3
              IF (KSTAR(IK).EQ.0) IP = 1
              W(K) = WW(IP)
              Q(K) = QQ(IP)
              IF (KSTAR(IK).LE.2.OR.KSTAR(IK).EQ.7) THEN
                  rg2(k) = 0.1
	      ELSE IF (KSTAR(IK).EQ.4) THEN
                  CM(K,IC) = MIN(0.89D0*BODY(IK),CM(K,IC))
                  rg2(k)= 0.1*(1.0 - CM(K,IC)/BODY(IK))
              ELSE
	          rg2(k)= 0.21
              END IF
          END IF
          TL = TWOPI*RADIUS(IK)*SQRT(RADIUS(IK)/BODY(IK)/W(K))
          IF (KSTAR(IK).GE.11) QL = 1.0D+10
          AT0(K) = 1.0/(QL*TL)
   10 CONTINUE
*
*       Form mass, radius & pericentre ratio.
      IF (RADIUS(I1).GE.RADIUS(I2)) THEN
          M21 = BODY(I2)/BODY(I1)
          R21 = RADIUS(I2)/RADIUS(I1)
	  RP1 = RP(IC)/RADIUS(I1)
	  rad = radius(i1)
      ELSE
	  M21 = BODY(I1)/BODY(I2)
          R21 = RADIUS(I1)/RADIUS(I2)
	  RP1 = RP(IC)/RADIUS(I2)
	  rad = radius(i2)
      END IF
*
*       Define initial angular momentum from the scaled semi-major axis.
      semi0 = RP1/(1.0 - ES0)
*
*       Form the initial mean motion in N-body units.
      omeq = sqrt(BODY(I)/(rad*semi0)**3)
*
*       Convert from angular momentum to omega (denoted spin1 & spin2).
      IF (KSTAR(J1).LE.2.OR.(KSTAR(J1).GE.7.AND.KSTAR(J1).NE.9)) THEN
          SPIN1 = SPIN(J1)/(rg2(1)*BODY(J1)*RADIUS(J1)**2)
      ELSE
          KW = KSTAR(J1)
          M0 = BODY0(J1)*SMU
          MC1 = CM(1,IC)*SMU
          IF (MC1.LE.0.0D0.OR.MC1.GT.M0) THEN
              MC1 = 0.3 + 0.1*(KSTAR(J1) - 3)
              IF(KW.EQ.9) MC1 = MIN(0.3D0,0.95*M0)
              CM(1,IC) = MC1/ZMBAR
          END IF
          ZDUM = 2.0D0
          RC1 = CORERD(KW,MC1,M0,ZDUM)/SU
          SPIN1 = SPIN(J1)/(rg2(1)*BODY(J1)*RADIUS(J1)**2 +
     &                      0.21*MC1/SMU*RC1**2)
      END IF
      IF (KSTAR(J2).LE.2.OR.(KSTAR(J2).GE.7.AND.KSTAR(J2).NE.9)) THEN
          spin2 = SPIN(J2)/(rg2(2)*BODY(J2)*RADIUS(J2)**2)
      ELSE
          KW = KSTAR(J2)
          M0 = BODY0(J2)*SMU
          MC2 = CM(2,IC)*SMU
          IF (MC2.LE.0.0D0.OR.MC2.GT.M0) THEN
              MC2 = 0.3 + 0.1*(KSTAR(J1) - 3)
              IF(KW.EQ.9) MC2 = MIN(0.3D0,0.95*M0)
              CM(2,IC) = MC2/ZMBAR
          END IF
          ZDUM = 2.0D0
          RC2 = CORERD(KW,MC2,M0,ZDUM)/SU
          spin2 = SPIN(J2)/(rg2(2)*BODY(J2)*RADIUS(J2)**2 +
     &                      0.21*MC2/SMU*RC2**2)
      END IF
*
*       Set synchronous spin for degenerate stars.
      DS10 = (spin1 - omeq)/omeq
      DS20 = (spin2 - omeq)/omeq
      IF (KSTAR(J1).GE.11) spin1 = omeq
      IF (KSTAR(J2).GE.11) spin2 = omeq
      DS1 = (spin1 - omeq)/omeq
      DS2 = (spin2 - omeq)/omeq
*
*       Skip integration if both spins are synchronous.
      IF (ABS(DS1).LT.0.01.AND.ABS(DS2).LT.0.01) THEN
          IF (ABS(DS20).GT.0.01) THEN
              WRITE (6,20)  KSTAR(I1), KSTAR(I2), A0*SU, DS1, DS20, OMEQ
   20         FORMAT (' SYNCH    K* A DS1 DS2 omeq ',2I4,F9.3,1P,3E10.2)
          END IF
          ISYNCH = 1
          GO TO 50
      END IF
*
*       Scale the spins by mean motion and define angular momentum.
      spin10=spin1/omeq
      spin20=spin2/omeq
      angmom0=(m21/(1+m21))*semi0**2*sqrt(1-es0**2)+rg2(1)*spin10+
     &         m21*r21**2*rg2(2)*spin20
*
*	Evaluate damping coefficients (Mardling & SJA, M.N. 321, 398, 2001).
      cf = 54.0*twopi/5.0
      C3 = (cf/9.0)*(AT0(1)*(Q(1)/W(1))**2*M21**2)/rg2(1)/semi0**6
      C4 = (cf/9.0)*(AT0(2)*(Q(2)/W(2))**2/M21**2)*R21**6/rg2(2)/
     &      semi0**6
*
*       Skip on zero interval.
      IF (time - time0.LE.0.0D0) GO TO 70
*
*       Obtain dominant terms of the derivatives (cf. routine HUT/DERIV2).
      udot1=C3
      udot2=C4
*
*       Choose the step from smallest time-scale (! time0 < time possible).
      taux = max(abs(1.0/udot1),abs(1.0/udot2))
      nstep = 1 + 100.0*sqrt(ABS(time - time0)/taux)
      nstep = min(nstep,100)
*       Increase number of steps on slow primary rotation.
      IF (DS1.LT.-0.2) nstep = nstep + 10
      DT = TIME - TIME0
*
*       Reduce integration interval for slow evolution or after merger.
      IF (DT.GT.1.0) THEN
          JPAIR = -IPAIR
          CALL TRFLOW(JPAIR,DTR)
          M1 = BODY(J1)*SMU
          CALL TRDOT(J1,DTM,M1)
          DTM = MAX(DTM,0.1D0)
          DT = 0.1*MIN(DTR,DTM,10.0D0)
          IF (ABS(DS2).GT.0.05) THEN
              DT = 0.1*DT
              nstep = nstep + 50
          END IF
          TK = DAYS*A0*SQRT(A0/BODY(I))
          WRITE (6,30)  NAME(J1), DS1, DS2, DT, TK
   30     FORMAT (' REDUCE    SYNCH    NM DS1 DS2 DT P ',I6,1P,4E10.2)
      END IF
*
   40 dtau1 = ABS(DT)/float(nstep)
*
*       Integrate equations for angular velocities only.
      call hut2(spin10,spin20,spin1,spin2,nstep,dtau1)
*
*       Re-scale the semi-major axis and angular velocities to N-body units.
      semi = rad*semi0*semi
      spin1 = omeq*spin1
      spin2 = omeq*spin2
      ecc = es0
      omeq = SQRT(BODY(I)/SEMI**3)
      DS1 = (spin1 - omeq)/omeq
      DS2 = (spin2 - omeq)/omeq
      NSYNCH = NSYNCH + 1
*
      IF (spin1.LT.0.0.OR.spin2.LT.0.0) THEN
          WRITE (6,42)  ITRY, nstep, DT, spin1, spin2, omeq
   42     FORMAT (' REPEAT    SYNCH    IT # DT S1 S2 omeq ',
     &                                 2I5,1P,4E10.2)
          ITRY = ITRY + 1
          nstep = nstep + 50
          IF (ITRY.LT.4) GO TO 40
          STOP
      END IF
      ITRY = 0
*
      IF ((DS10.LT.0.0.AND.DS1.GT.0.01).OR.
     &    (DS20.LT.0.0.AND.DS2.GT.0.01)) THEN
           ITER = ITER + 1
           nstep = nstep + 50
           IF (ITER.GT.1) DT = 0.5*DT
           IF (ITER.LT.4) GO TO 40
           WRITE (6,45)  nstep, DS10, DS1, DS20, DS2, DT
   45      FORMAT (' LIMIT    SYNCH    # DS10, DS1, DS20, DS2, DT ',
     &                                 I4,1P,5E10.2)
      END IF
*
      IF ((DS1.GT.0.0.AND.DS10.LT.-0.01).OR.
     &    (DS2.GT.0.0.AND.DS20.LT.-0.01)) THEN
           ITER = ITER + 1
           nstep = nstep + 50
           IF (ITER.GT.1) DT = 0.5*DT
           IF (ITER.LT.4) GO TO 40
           WRITE (6,45)  nstep, DS10, DS1, DS20, DS2, DT
      END IF
      ITER = 0
*
*       Convert back to angular momenta.
   50 IF (KSTAR(J1).LE.2.OR.(KSTAR(J1).GE.7.AND.KSTAR(J1).NE.9)) THEN
          SPIN(J1) = rg2(1)*BODY(J1)*RADIUS(J1)**2*spin1
      ELSE
          SPIN(J1) = (rg2(1)*BODY(J1)*RADIUS(J1)**2 +
     &                0.21*MC1/SMU*RC1**2)*spin1
      END IF
      IF (KSTAR(J2).LE.2.OR.(KSTAR(J2).GE.7.AND.KSTAR(J2).NE.9)) THEN
          SPIN(J2) = rg2(2)*BODY(J2)*RADIUS(J2)**2*spin2
      ELSE
          SPIN(J2) = (rg2(2)*BODY(J2)*RADIUS(J2)**2 +
     &                0.21*MC2/SMU*RC2**2)*spin2
      END IF
      IF (ISYNCH.GT.0) THEN
          ISYNCH = 0
          GO TO 70
      END IF
*
      ERROR = (SEMI - A0)/A0
*       Obtain the tidal contributions from integration.
      A0 = RP0/(1.0 - ES0)
      DH = -0.5*BODY(I)*(1.0/SEMI - 1.0/A0)
*
*       Update energy and semi-major axis.
      H(IPAIR) = H(IPAIR) + DH
      SEMI = -0.5*BODY(I)/H(IPAIR)
*
*       Set new pericentre from final elements (delay updating).
      IF (LIST(1,I1).EQ.0) THEN
          QPERI = semi*(1.0 - ecc)
      ELSE
          QPERI = RP0*(1.0 + ES0)/(1.0 + ECC)
      END IF
*
*       Form KS coordinate scaling factor from pericentre ratio.
      C1 = SQRT(QPERI/RP0)
*
*       Set KS velocity scaling from energy relation with RP0 instead of R.
      C2 = SQRT((BODY(I) + H(IPAIR)*QPERI)/(BODY(I) + HI*RP0))
*
*       Scale KS variables to yield the prescribed elements (set TDOT2 >= 0).
      R(IPAIR) = 0.0D0
      TDOT2(IPAIR) = 0.0
      DO 60 K = 1,4
          U(K,IPAIR) = C1*U(K,IPAIR)
          UDOT(K,IPAIR) = C2*UDOT(K,IPAIR)
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
          TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0*U(K,IPAIR)*UDOT(K,IPAIR)
   60 CONTINUE
      TDOT2(IPAIR) = MAX(TDOT2(IPAIR),0.0D0)
*
*       Rectify the orbit to prevent eccentricity growth.
      CALL KSRECT(IPAIR)
*
*       Perform energy correction to maintain conservation.
      ZMU = BODY(I1)*BODY(I2)/BODY(I)
      ECOLL = ECOLL + ZMU*(HI - H(IPAIR))
*
*       Initialize KS for perturbed motion.
      IF (LIST(1,I1).GT.0) THEN
          CALL RESOLV(IPAIR,1)
          IMOD = KSLOW(IPAIR)
          CALL KSPOLY(IPAIR,IMOD)
      END IF
*
      DA = (SEMI - A0)/SEMI
      DS1 = (spin1-omeq)/spin1
      DS2 = (spin2-omeq)/spin2
      TS1 = MIN(TEV(J1),TEV(J2))
      WRITE (7,65) nstep,time-time0,udot1,udot2,DS1,DS2,DA,TS1
   65 FORMAT (' HUT2    # t-t0 ud DS DA/A TM ',
     &                  I5,F9.3,1P,5E9.1,0P,F10.2)
*
*       Include magnetic or gravitational braking for small separations.
*     IF (SEMI*SU.LT.10.0) THEN
*         DT = TIME - TIME0
*         CALL BRAKE(IPAIR,DT)
*         IF (IPHASE.LT.0) GO TO 100
*     END IF
*
*       Update look-up time from radial expansion or Roche interval.
   70 JPAIR = -IPAIR
      CALL TRFLOW(JPAIR,DTR)
      M1 = BODY(J1)*SMU
      CALL TRDOT(J1,DTM,M1)
      DTM = MAX(DTM,0.1D0)
      TOSC(IC) = TIME
      TEV(I) = TIME + 0.1*MIN(DTR,DTM)
      TEV(I) = MIN(TEV(J1),TEV(J2),TEV(I))
      ICOUNT = ICOUNT + 1
*
      SEMI = -0.5*BODY(I)/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(SEMI*BODY(I))
      ECC = SQRT(ECC2)
*       Correct for small eccentricity drift.
      IF (ABS(ECC - ES0).GT.0.001*ES0) THEN
          CALL DEFORM(IPAIR,ECC,ES0)
          ECC = ES0
          QPERI = SEMI*(1.0 - ES0)
      END IF
*       Save final reference values.
      RP(IC) = QPERI
      ES(IC) = ECC
      TOSC(IC) = TIME
*
      IF (ECC.GT.0.01) THEN
          WRITE (6,80)  NAME(I1), ECC, nstep, dtau1, ERROR,C1-1.0,C2-1.0
   80     FORMAT (' DANGER!    NAM E # dtau DA/A DC1 DC2 ',
     &                         I6,F8.4,I5,1P,5E10.2)
          STOP
      END IF
*
      IF (ICOUNT.LE.300) THEN
          DT = TEV(I) - TIME
          WRITE (7,85)  NAME(J1), nstep, KSTAR(J1), KSTAR(J2),
     &                  SEMI*SU, DA, omeq, spin1, spin2, DT
   85     FORMAT (' SYNCH    NM # K* A DA/A om s1 s2 DT ',
     &                       2I8,2I4,F7.2,1P,5E10.2)
          CALL FLUSH(7)
      END IF
*
      IF (spin1.LT.0.0) THEN
          WRITE (6,90)  NAME(J1), KSTAR(J1), spin1, TIME-TIME0
   90     FORMAT (' DANGER!    NEGATIVE SPIN    NM K* s1 T-T0 ',
     &                         I6,I4,1P,2E10.2)
          STOP
      END IF
*
  100 RETURN
*
      END
