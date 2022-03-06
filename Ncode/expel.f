      SUBROUTINE EXPEL(J1,J2,ICASE)
*
*
*       Common envelope stage of interacting stars.
*       -------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &               BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &               RP(NTMAX),ES(NTMAX),CM(2,NTMAX),IOSC(NTMAX),
     &               NAMEC(NTMAX)
      REAL*8 M01,M1,R1,AJ1,M02,M2,R2,AJ2,SEP,MI(2)
      REAL*8 TSCLS(20),LUMS(10),GB(10),TM,TN,LUM1,LUM2,MC1,MC2,RCC
      REAL*8 JSPIN1,JSPIN2,MENV,RENV,K2
      LOGICAL COALS
*
*
*       Define global indices such that body #I1 is giant-like.
      IF(KSTAR(J1).GE.2.AND.KSTAR(J1).LE.9.AND.KSTAR(J1).NE.7)THEN
          I1 = J1
          I2 = J2
      ELSE
          I1 = J2
          I2 = J1
      ENDIF
*
*       Update the stars to latest previous time.
      TEV1 = MAX(TEV0(I1),TEV0(I2))
*
*       Specify basic parameters for both stars.
      M01 = BODY0(I1)*ZMBAR
      M1 = BODY(I1)*ZMBAR
      MC1 = 0.D0
      AJ1 = TEV1*TSTAR - EPOCH(I1)
      JSPIN1 = SPIN(I1)*SPNFAC
      KW1 = KSTAR(I1)
      M02 = BODY0(I2)*ZMBAR
      M2 = BODY(I2)*ZMBAR
      MC2 = 0.D0
      AJ2 = TEV1*TSTAR - EPOCH(I2)
      JSPIN2 = SPIN(I2)*SPNFAC
      KW2 = KSTAR(I2)
*
*       Save current semi-major axis.
      IPAIR = KSPAIR
      I = N + IPAIR
      SEMI0 = -0.5d0*BODY(I)/H(IPAIR)
      SEP = SEMI0*SU
      ECC2 = (1.D0-R(IPAIR)/SEMI0)**2+TDOT2(IPAIR)**2/(SEMI0*BODY(I))
      ECC = SQRT(ECC2)
      ECC0 = ECC
*
*       Perform common envelope evolution (note: SEP in SU).
      CALL COMENV(M01,M1,MC1,AJ1,JSPIN1,KW1,
     &            M02,M2,MC2,AJ2,JSPIN2,KW2,ECC,SEP,COALS)
*
*       Obtain consistent radii for the stars (skip #I2 on coalescence).
      CALL star(KW1,M01,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
      CALL hrdiag(M01,AJ1,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &            R1,LUM1,KW1,MC1,RCC,MENV,RENV,K2)
      IF(KW1.NE.KSTAR(I1))THEN
         WRITE(38,*)' EXPEL TYPE CHANGE1 ',KSTAR(I1),KW1
      ENDIF
      IF(COALS)THEN
          KW2 = KSTAR(I2)
          R2 = 0.d0
          LUM2 = 0.d0
          JSPIN2 = 0.d0
      ELSE
          CALL star(KW2,M02,M2,TM,TN,TSCLS,LUMS,GB,ZPARS)
          CALL hrdiag(M02,AJ2,M2,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                R2,LUM2,KW2,MC2,RCC,MENV,RENV,K2)
          IF(KW2.NE.KSTAR(I2))THEN
             WRITE(38,*)' EXPEL TYPE CHANGE2 ',KSTAR(I2),KW2
          ENDIF
      END IF
*
*       Skip small mass loss unless coalescence (May 2003).
      DMS = M1 + M2 - BODY(I)*SMU
      IF (ABS(DMS).LT.1.0D-10.AND..NOT.COALS) THEN
*         WRITE (77,7)  TIME, ECC, DMS, SEMI0*SU-SEP
*   7     FORMAT (' FAIL!    T ECC DMS DSEP ',F12.4,F8.4,1P,2E10.2)
*         CALL FLUSH(77)
          IPHASE = 0
          GO TO 60
      END IF
*
*       Copy new values of radius, luminosity & semi-major axis.
      RADIUS(I1) = R1/SU
      RADIUS(I2) = R2/SU
      ZLMSTY(I1) = LUM1
      ZLMSTY(I2) = LUM2
      SPIN(I1) = JSPIN1/SPNFAC
      SPIN(I2) = JSPIN2/SPNFAC
      SEMI = SEP/SU
*       Choose separation as semi-major axis in hyperbolic case.
      IF (SEMI.LT.0.0) SEMI = R(IPAIR)
*
*       Set reduced mass & binding energy and update masses & time T0(J1).
      ZMU0 = BODY(I1)*BODY(I2)/BODY(I)
      HI = H(IPAIR)
      BODY0(I1) = M01/ZMBAR
      BODY0(I2) = M02/ZMBAR
      T0(J1) = TIME
*       Add current energy (final state: coalescence or close binary).
      ECOLL = ECOLL + ZMU0*HI
*
*       Distinguish between coalescence and surviving binary.
      EPOCH(I1) = TEV1*TSTAR - AJ1
      IF(COALS)THEN
          MI(1) = M1
          MI(2) = M2
          BODY0(I2) = 0.D0
*       Check for TZ object formed by CE evolution.
          IF(KSTAR(I2).GE.13.AND.KW1.GE.13)THEN
              NTZ = NTZ + 1
              WRITE (6,5)  M1, M2
    5         FORMAT (' NEW TZ    M1 M2 ',2F7.2)
          ENDIF
          CALL COAL(IPAIR,KW1,KW2,MI)
      ELSE
*
*       Update evolution times.
          EPOCH(I2) = TEV1*TSTAR - AJ2
          TEV(I1) = TEV1
          TEV0(I1) = TEV(I1)
          TEV(I2) = TEV1
          TEV0(I2) = TEV(I2)
*       Shrink the orbit in case of no coalescence.
          BODY(I1) = M1/ZMBAR
          BODY(I2) = M2/ZMBAR
*       Form mass loss & new binary mass and re-define binding energy.
          DM = BODY(I) - (BODY(I1) + BODY(I2))
          ZMASS = ZMASS - DM
          BODY(I) = BODY(I) - DM
          H(IPAIR) = -0.5d0*BODY(I)/SEMI
*
*       Subtract final binding energy of surviving binary.
          ZMU = BODY(I1)*BODY(I2)/BODY(I)
          ECOLL = ECOLL - ZMU*H(IPAIR)
*
*       Set argument for new pericentre velocity (need to have H*A*(1-E)).
          RX = R(IPAIR)/(1.0 - ECC)
*       Modify KS variables consistent with new eccentricity and energy H.
          CALL EXPAND(IPAIR,RX)
          CALL RESOLV(IPAIR,1)
*
          WRITE (6,10)  NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2), KW1,
     &                  KW2, M1, M2, DM*ZMBAR, ECC, R1, R2, SEMI0*SU,
     &                  SEMI*SU
   10     FORMAT (' NEW CE    NAM K0* K* M1 M2 DM E R1 R2 A0 A ',
     &                        2I6,4I3,3F5.1,F8.4,4F7.1)
*
*       Copy new types and reset IPHASE (prevents CALL CHAIN from DELAY).
          KSTAR(I1) = KW1
          KSTAR(I2) = KW2
          IPHASE = 0
*
*       Obtain neighbour list for force corrections and polynomials.
          NNB = LIST(1,I)
          DO 11 L = 1,NNB+1
              ILIST(L) = LIST(L,I)
   11     CONTINUE
*
*       Ensure rare case of massless remnant will escape at next output.
         IF (KW1.EQ.15.OR.KW2.EQ.15) THEN
            CALL KSTERM
            IPHASE = -1
            I = IFIRST
            IF (BODY(I).GT.0.0D0) I = IFIRST + 1
            T0(I) = TADJ + DTADJ
            STEP(I) = 1.0D+06
            RI = SQRT(X(1,I)**2 + X(2,I)**2 + X(3,I)**2)
            VI = SQRT(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2)
            DO 12 K = 1,3
               X0(K,I) = MIN(1.0d+04+X(K,I),1000.0*RSCALE*X(K,I)/RI)
               X(K,I) = X0(K,I)
               X0DOT(K,I) = SQRT(0.004*ZMASS/RSCALE)*XDOT(K,I)/VI
               XDOT(K,I) = X0DOT(K,I)
               F(K,I) = 0.0D0
               FDOT(K,I) = 0.0D0
               D2(K,I) = 0.0D0
               D3(K,I) = 0.0D0
   12       CONTINUE
            WRITE (6,13)  NAME(I), KSTAR(I), BODY(I)*ZMBAR
   13       FORMAT (' MASSLESS GHOST    NAM K* M ',I6,I4,1P,E10.2)
*       Include special treatment for velocity kick of KS binary.
          ELSE IF (KW1.GE.13) THEN
*
*       Evaluate binary disruption velocity and introduce velocity kick.
              VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
              KSTAR(I1) = -KSTAR(I1)
              CALL KICK(IPAIR,0)
*
*       Include potential and kinetic energy corrections due to mass loss.
              CALL POTI(I,POTJ)
              ECDOT = ECDOT + DM*POTJ + 0.5d0*DM*VI2
*
*       See whether tidal terms should be included (standard or scaled).
              IF (KZ(14).GT.0) THEN
                  IF (KZ(14).LE.2) THEN
                      ECDOT = ECDOT - 0.5d0*DM*(TIDAL(1)*X(1,I)**2 +
     &                                          TIDAL(3)*X(3,I)**2)
                  ELSE
                      BODY(I) = BODY(I) + DM
                      CALL XTRNLV(I,I)
                      ECDOT = ECDOT + HT
                      BODY(I) = BODY(I) - DM
                      CALL XTRNLV(I,I)
                      ECDOT = ECDOT - HT
                  END IF
              END IF
*
*       Form new c.m. variables and update GRAPE.
              T0(I) = TIME
              DO 20 K = 1,3
                  X(K,I) = (BODY(I1)*X(K,I1) + BODY(I2)*X(K,I2))/BODY(I)
                  XDOT(K,I) = (BODY(I1)*XDOT(K,I1) +
     &                         BODY(I2)*XDOT(K,I2))/BODY(I)
                  X0(K,I) = X(K,I)
                  X0DOT(K,I) = XDOT(K,I)
   20         CONTINUE
*
*       Obtain new F & FDOT and time-steps for neighbours.
              DO 30 L = 2,NNB+1
                  J = ILIST(L)
                  DO 25 K = 1,3
                      X0DOT(K,J) = XDOT(K,J)
   25             CONTINUE
                  CALL FPOLY1(J,J,0)
                  CALL FPOLY2(J,J,0)
   30         CONTINUE
              TPREV = TIME - STEPX
*
*       Terminate KS binary and assign kick velocity to single star #I.
              I = I1 + 2*(NPAIRS - IPAIR)
              CALL KSTERM
              CALL KICK(I,1)
*       Initialize new KS polynomials after velocity kick (H > 0 is OK).
              ICOMP = IFIRST
              JCOMP = IFIRST + 1
              IPHASE = 1
              CALL KSREG
*
*       Provide diagnostic output if binary survives the kick.
              I = NTOT
              KSTAR(I) = 0
              IF (H(IPAIR).LT.0.0) THEN
                  SEMI = -0.5d0*BODY(I)/H(IPAIR)
                  WRITE (6,35)  R(IPAIR)/SEMI, SEMI*SU, BODY(I)*ZMBAR,
     &                          (XDOT(K,I)*VSTAR,K=1,3)
   35             FORMAT (' NS BINARY    R/A A M V ',3F6.2,3F7.1)
              END IF
          ELSE
*       Set new binary indicator and update Roche look-up time.
              KSTAR(I) = 0
              CALL TRFLOW(IPAIR,DTR)
              TEV(I) = TIME + DTR
*
*       Update chaos variables for routine SYNCH after circularization.
              IF (ECC.LE.0.001D0) THEN
                 IC = 0
*       Note: NAMEC may be set in routine CHAOS without call to SPIRAL.
                 DO 36 L = 1,NCHAOS
                     IF (NAME(I).EQ.NAMEC(L)) IC = L
   36            CONTINUE
*       Include safety test on no membership just in case (not sure!).
                 IF (IC.EQ.0) THEN
                     NCHAOS = NCHAOS + 1
                     IC = NCHAOS
                     NAMEC(IC) = NAME(I)
                     IF (NCHAOS.GT.NTMAX) THEN
                          WRITE (6,38)  NAME(I1), NCHAOS, ECC
   38                     FORMAT (' FATAL ERROR!    EXPEL    NM NCH E ',
     &                                              I6,I4,F8.4)
                          STOP
                      END IF
                 END IF
                 IF (IC.GT.0) THEN
                     RP(IC) = SEMI*(1.0 - ECC)
                     ES(IC) = ECC
                     TOSC(IC) = TIME
                     KSTAR(I) = 10
                 END IF
              END IF
*
*       Include prediction to current time before new KS and possible FPOLYI.
              IF (T0(I).LT.TIME.AND.DM.GT.0.0D0) THEN
                  CALL XVPRED(I,-2)
                  T0(I) = TIME
                  CALL DTCHCK(TIME,STEP(I),DTK(40))
                  DO 40 K = 1,3
                      X0(K,I) = X(K,I)
                      X0DOT(K,I) = XDOT(K,I)
   40             CONTINUE
                  IPHASE = -1
              END IF
*
*       Re-initialize the KS solution.
              CALL KSPOLY(IPAIR,1)
*
*       Correct neighbour forces and update system energy loss due to DM.
              IF (DM.GT.0.0) THEN
                  IF (DM*SMU.LT.0.05) THEN
                      CALL FICORR(I,DM)
                  ELSE
                      CALL FCORR(I,DM,0)
                  END IF
              END IF
*
              IF (DM*SMU.GT.0.1) THEN
*       Include body #I at the end (counting from location #2).
                  NNB2 = NNB + 2
                  ILIST(NNB2) = I
*
*       Obtain new F & FDOT and time-steps (note FPOLYI may over-write I).
                  DO 50 L = 2,NNB2
                      J = ILIST(L)
                      IF (L.EQ.NNB2) THEN
                          J = I
                      ELSE
                          CALL XVPRED(J,-2)
                          CALL DTCHCK(TIME,STEP(J),DTK(40))
                          DO 45 K = 1,3
                              X0DOT(K,J) = XDOT(K,J)
   45                     CONTINUE
                      END IF
                      CALL FPOLY1(J,J,0)
                      CALL FPOLY2(J,J,0)
   50             CONTINUE
              END IF
              TPREV = TIME - STEPX
          END IF
      END IF
*
*       Reverse case indicator for exit from collision routine.
   60 ICASE = -ICASE
*
      RETURN
*
      END
