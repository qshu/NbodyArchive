      SUBROUTINE CHAIN(ISUB)
*
*
*       Perturbed chain regularization.
*       -------------------------------
*
*       Method of Mikkola & Aarseth, Celestial Mechanics 47, 375.
*       .........................................................
*
      INCLUDE 'COMMON1.CH'
      INCLUDE 'COMMON2.CH'
      REAL*8  G0(3),Y(NMX8)
      LOGICAL  CHECK
      REAL*8  M,MIJ
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NMX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NMX,5),ISYS(5)
      COMMON/INTFAC/  LX,LE,LP,LV,LT,J10
      COMMON/ECHAIN/  ECH
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
      COMMON/CALLS/  TPR,ICALL,NFN,NREG
      SAVE
      EXTERNAL CHMOD
*
*
*       Main variables for chain regularization
*       ***************************************
*
*       -------------------------------------------------------------------
*       CHTIME  Local physical time (from CHAIN code).
*       CM(1-7) Coordinates, velocities & total mass of system.
*       CM(8)   Total energy of N-body system (copied in routine CHINIT).
*       CM(9)   Binary energy change (copied to CHCOLL in CHTERM).
*       ECH     Total energy of perturbed system (N-body interface).
*       ENERGY  Total energy.
*       EPS     Tolerance for DIFSY1 (1.0E-08 is recommended).
*       ICALL   Pericentre indicator (activated if R(IMIN) < EPSR2**0.5).
*       ICOLL   Collision indicator (activated in routine DERQP4).
*       IPERT   Perturbation indicator (=0: once per call; =1: every call).
*       I1-I4   Indices of final configuration (I1 & I2 is closest binary).
*       KZ27    Tidal dissipation option (not implemented yet).
*       KZ30    Diagnostic option (copied in CHINIT; full output if > 1).
*       M       Particle mass (CM(7) = sum M(I), I = 1,N).
*       NAMEC   Particle identity (initialized to global name).
*       NFN     Number of function calls.
*       NPERT   Number of perturbers (for diagnostic output).
*       NREG    Number of regularization switches.
*       NSTEP1  Number of DIFSY1 calls.
*       P       Regularized momenta.
*       Q       Regularized coordinates.
*       RINV    Inverse chain distances.
*       RCOLL   Minimum two-body separation (not activated yet).
*       RGRAV   Gravitational radius ((sum M(I)*M(J))/ABS(ENERGY)).
*       RMAXS   Maximum size of unperturbed configuration.
*       RMAXC   Maximum of unperturbed & initial size (+ 20 per cent).
*       RSUM    Sum of all chain distances.
*       STEPS   Current integration interval (set by routine INTGRT).
*       TCR     Local crossing time ((sum M(I))**2.5/ABS(2*ENERGY)**1.5).
*       TIMEC   Local physical time in scaled units (= CHTIME).
*       TMAX    Maximum integration time (based on c.m. step).
*       TS      Global physical time at each return to N-body integration.
*       X       Particle coordinates (X(1), X(2), X(3) is X, Y, Z).
*       X       Velocity components.
*       -------------------------------------------------------------------
*
*
*       Save termination indicator and check for restart.
      ITERM = ISUB
      IF (ISUB.GT.0) THEN
*       Synchronize next time interval with subsystem step.
          TMAX = TIMEC + STEPS(ISUB)
          ITER = 0
*       Choose next step from first order relation for T' (but < STEP1).
          IF (STEPS(ISUB).GT.0.0D0) THEN
              STEP = MIN(STEPS(ISUB)/TPR2,ABS(STEP1))
          ELSE
              STEP = 1.0E-06*MIN(ABS(STEP1),ABS(STEP))
          END IF
          GO TO 30
      END IF
*
*       Copy initial conditions from N-body COMMON and prepare chain.
      CALL CHINIT(ISUB)
*
*       Initialize diagnostic & COMMON variables.
      TIMEC = 0.0D0
      RCOLL = 100.0
      ICALL = 0
      ICOLL = 0
      ITER = 0
      NSTEP1 = 0
      NREG = 0
      NFN = 0
      KCASE = 0
      J10 = 10
*
*       Specify the tolerance for DIFSY1.
      EPS = 1.0E-08
*
*       Initialize subsystem time and set dummy variable for DIFSY.
      CHTIME = 0.0D0
      STIME = 0.0D0
*
*       Initialize chain regularization (new case or modified chain).
   10 NEQ = 8*N
      Y(NEQ) = CHTIME
      IF (KCASE.GT.0) THEN
          TMAX = TIMEC + STEPS(ISUB)
          ITER = 0
      END IF
*
*       Define indices for DIFSY1 calls to derivative routine.
      NC = N - 1
      LX = 4*NC + 1
      LE = 4*NC + 4
      LP = 4*NC + 5
      LV = 8*NC + 5
      LT = 8*NC + 8
*       Ensure routine XTPERT uses correct subsystem index (ISYS(5) is free).
      ISYS(5) = ISUB
*
*       Evaluate constants of the motion.
      CALL CONST(X,V,M,N,ENER0,G0,ALAG)
*
*       Select chain indices and transform from chain to KS variables.
      CALL SELECT
      CALL TRANSQ
*
*       Suppress XTRNLU & UG if external potential is absent.
*     CALL XTRNLU(X,M,N,UG)
*
*       Define total energy (including any c.m. motion & external potential).
      ENERGY = ENER0 - 0.5D0*MASS*(CMV(1)**2 + CMV(2)**2 + CMV(3)**2)
*    &               + UG
*
*       Copy whole input array (for DIFSY call).
      CALL YCOPY(Y)
*
*       Find sum of mass products and set current gravitational radius.
      SUM = 0.0D0
      DO 20 L = 1,N-1
          DO 15 K = L+1,N
              SUM = SUM + MIJ(L,K)
   15     CONTINUE
   20 CONTINUE
      RGRAV = SUM/ABS(ENERGY)
*
*       Set current crossing time and specify initial step using T' = 1/L.
      TCR = MASS**2.5/ABS(2.0D0*ENERGY)**1.5
      STEPIN = EPS**0.2*MASS**2.5/ALAG**0.5
      TPR1 = 1.0/ALAG
*
      IF (KZ30.GT.1) THEN
          WRITE (6,25)  N, NPERT, ENERGY, RSUM, RGRAV, TCR, RMAXS(ISUB)
   25     FORMAT (' NEW CHAIN   N NP E RSUM RGRAV TCR RMAXS ',
     &                          2I4,F10.5,F8.4,1P,3E9.1)
      END IF
*
*       Assign integration step (initial case or modified membership).
      IF (KCASE.EQ.0) THEN
          STEP = STEPIN
      ELSE
*       Base continuation step on remaining time interval (but <= STEP1).
          STEP = MIN(ALAG*ABS(TMAX - TIMEC),ABS(STEP1))
          KCASE = 0
*       Ensure termination if only two particles left or positive energy.
          IF (N.EQ.2.OR.ENERGY.GT.0.0) GO TO 40
          GO TO 31
      END IF
*
*       Evaluate initial binary energies.
      CALL RECOIL(0)
*
*       Perform regularized integration until termination or modification.
   30 KSTEPS = 0
      RMAXC = 3.0*RGRAV
*       Replace perturbed boundary radius with MIN(3*RGRAV,2*RSUM).
      IF (KCASE.GT.0) RMAXC = MIN(RMAXC,3.0*RGRAV,2.0*RSUM)
   31 KSTEPS = KSTEPS + 1
      STEP0 = STEP
*
*       Modify DIFSY integration order according to relative perturbation.
      JDIF = 10 - (1.0D+04*GPERT)**0.4
*     IF (JDIF.LT.J10) THEN
*         J10 = MAX(4,JDIF)
*     ELSE IF (JDIF.GT.J10) THEN
*         J10 = MIN(10,J10+1)
*     END IF
*
*       Switch on indicator for external perturbations (first function call).
      IPERT = 1
      ICALL = 1
      TLAST = CHTIME
      CALL DIFSY1(NEQ,EPS,STEP,STIME,Y)
*
*     IF (N.GT.3) THEN
*         WRITE (6,32)  (1.0/RINV(K),K=1,N-1)
*  32     FORMAT (' CHAIN:   R  ',1P,5E10.2)
*         CALL FLUSH(6)
*     END IF
*
      IF (STEP.EQ.0.0D0) THEN
          WRITE (6,*) ' Stepsize = 0!', char(7)
          STOP
      END IF
*
*       Copy current physical time and save COMMON variables.
      CHTIME = Y(NEQ)
      CALL YSAVE(Y)
      CALL TRANSX
*
*       Set last physical integration interval and remaining interval.
      DTR = TMAX - CHTIME
      STEP1 = STEP
      IF (ITER.EQ.0.AND.DTR.LT.0.0) THEN
          STEP = DTR/TPR
          ITER = ITER + 1
*         WRITE (6,33)  TPR, CHTIME-TLAST, DTR, DTR/TPR, STEP0
*  33     FORMAT (' REV:  TPR DT DTR DTR/TPR STEP0 ',1P,5E9.1)
          GO TO 31
      END IF
      ITER = ITER + 1
*
*       Define average of T' for step prediction and save current value.
      TPR2 = 0.5*(TPR1 + TPR)
      TPR1 = TPR
*       Form new step from remaining interval (add small term and use MAX).
      STEP = MAX(ABS(DTR)/TPR2,1.0D-02*STEPIN)
      STEP = MIN(STEP,2.0*ABS(STEP0),ABS(STEP1))
*
*       Check switching condition (Note: RINV not known after switch).
      ISW = 0
      CALL SWCOND(CHECK)
      IF (CHECK) THEN
          CALL SWITCH(Y)
          NREG = NREG + 1
          ISW = 1
      END IF
*
*       Check termination criteria (T > TMAX or RSUM > RMAXC).
      IF ((CHTIME.GT.TMAX).OR.(RSUM.GT.RMAXC).OR.GPERT.GT.0.05) THEN
          CALL YSAVE(Y)
          CALL TRANSX
          ECH = ENERGY
          TIMEC = CHTIME
          NSTEP1 = NSTEP1 + KSTEPS
          IF (KZ30.GT.2) THEN
              WRITE (6,34)  KSTEPS, T0S(ISUB)+TIMEC, TIMEC-TMAX, RSUM,
     &                      RMAXC, (1.0/RINV(K),K=1,N-1)
   34         FORMAT (' CHAIN:  # T DT RS RM R',I3,2F8.4,2F6.3,1P,5E9.1)
          END IF
          IF (ISW.EQ.0.AND.(RSUM.GT.RMAXC.OR.GPERT.GT.0.05)) THEN
              CALL CHMOD(ISUB,KCASE)
              IF (KCASE.GT.0) THEN
                  CALL RECOIL(1)
                  GO TO 10
              END IF
              IF (KCASE.LT.0) GO TO 35
          END IF
          IF (CHTIME.LT.TMAX.AND.STEPS(ISUB).GT.0.0D0) GO TO 31
          GO TO 35
      ELSE
*       Exit chain integration if remaining interval is < 0.
          IF (DTR.LT.0.0) GO TO 35
          IF (STEPS(ISUB).GT.0.0D0) GO TO 31
      END IF
*
*       See whether temporary or actual termination (continue if N > 4).
   35 IF (N.GT.4.OR.(KCASE.EQ.0.AND.STEPS(ISUB).GT.0.0D0)) GO TO 100
      IF (TIMEC.GT.TMAX.AND.RSUM.LT.RMAXC) THEN
         IF (STEPS(ISUB).GT.0.0D0.OR.N.GT.4) GO TO 100
      END IF
*
*       Check for dominant binary.
   40 CALL RECOIL(2)
*
*       Set zero step to define termination (just in case).
      STEPS(ISUB) = 0.0D0
*
*       Transform to global variables and begin new KS (I1 & I2), I3 & I4.
      CALL CHTERM(ISUB)
*
*       Activate termination index for routine INTGRT.
      ITERM = -1
*
*       Update current time unless termination and set subsystem index.
  100 IF (ITERM.GE.0) TS(ISUB) = T0S(ISUB) + TIMEC
      ISUB = ITERM
*
      RETURN
*
      END
