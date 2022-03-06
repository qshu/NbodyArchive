# 1 "nbody6.F"
# 1 "<built-in>"
# 1 "<command line>"
# 1 "nbody6.F"
        PROGRAM NBODY6
*
* N B O D Y 6++
* *************
*
* Regularized AC N-body code with triple & binary collisions.
* --------------------------------------------------------
*
* Hermite integration scheme with block-steps (V 4.0.0 April/99).
* ------------------------------------------------------------------
*
* Developed by Sverre Aarseth, IOA, Cambridge.
* ............................................
* Message Passing Version NBODY6++ for Massively Parallel Systems
* Developed by Rainer Spurzem, ARI, Heidelberg
*
      INCLUDE 'common6.h'
      COMMON/STSTAT/ TINIT,NIR,NIB,NRGL,NKS
      EXTERNAL MERGE
*
# 39 "nbody6.F"
*
* Initialize the timer.
      CALL CPUTIM(ttota)
*
* Read start/restart indicator & CPU time.
      IF(rank.eq.0)READ (5,*) KSTART, TCOMP, TCRITp,
     * isernb,iserreg
*
# 58 "nbody6.F"
*
      IF (KSTART.EQ.1) THEN
*
* Read input parameters, perform initial setup and obtain output.
          CPU = TCOMP
          CALL START
          call cputim(tt7)
          CALL ADJUST
          call cputim(tt8)
          ttadj = ttadj + (tt8-tt7)*60.
      ELSE
*
* Read previously saved COMMON variables from tape/disc on unit 1.
          CALL MYDUMP(0,1)
*
          IF (NDUMP.GE.3) STOP
* Safety indicator preventing repeated restarts set in routine CHECK.
          CPU = TCOMP
          CPU0 = 0.0
* Set IPHASE = -1 for new NLIST in routine INTGRT (Hermite version).
          IPHASE = -1
*
* Initialize evolution parameters which depend on metallicity.
          CALL ZCNSTS(ZMET,ZPARS)
*
* Check reading modified restart parameters (KSTART = 3, 4 or 5).
          IF (KSTART.GT.2) THEN
              CALL MODIFY(KSTART)
          END IF
*
* Open all other files.
          CALL FILE_INIT(0)
*
* If no explicit new TCRIT given just go for another TCRIT of common block.
      TTOT = TIME + TOFF
      TCRIT = TTOT + TCRIT
      if(rank.eq.0)then
      WRITE (6,10) TTOT/TCR0, TIME/TCR0, TCRIT/TCR0, TTOT, TIME, TCRIT
      WRITE (6,20) DTADJ, DELTAT, TADJ, TNEXT, TCRIT, QE
      WRITE (6,30) ETAI, ETAR, ETAU, DTMIN, RMIN, NNBOPT
   10 FORMAT (' START AT TTOT/TIME ',2F16.8,' STOP INTENDED AT ',
     & F16.8,' TCR0',/,' START AT TTOT/TIME ',2F16.8,
     & ' STOP INTENDED AT ',F16.8,' NBODY-UNITS ',/)
   20 FORMAT (/,7X,'RESTART PARAMETERS:   DTADJ =',F7.3,'  DELTAT =',
     & F7.3,'   TADJ =',F7.3,'   TNEXT =',
     & F7.3,'  TCRIT =',F7.1,'  QE =',1PE9.1)
   30 FORMAT (/,7X,'                      ETAI =',F7.3,'  ETAR =',
     & F7.3,'  ETAU =',F7.3,'  DTMIN =',1PE9.1,
     & '  RMIN =',E9.1,' NNBOPT =',I5,/)
      end if
*
      END IF
*
* (R.Sp.)Set time flag and step number flags for beginning of run
      TINIT = TTOT
      NIR = NSTEPI
      NIB = NSTEPB
      NRGL = NSTEPR
      NKS = NSTEPU
*
      call cputim(tt2)
      ttinit = ttinit + (tt2-ttota)*60.
* Advance solutions until next output or change of procedure.
    1 CONTINUE
      call cputim(tt1)
*
      CALL INTGRT
*
      call cputim(tt2)
      ttint = ttint + (tt2-tt1)*60.
*
      IF (IPHASE.EQ.1) THEN
* Prepare new KS regularization.
      call cputim(tt1)
          CALL KSREG
          CALL FLUSH(6)
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.2) THEN
* Terminate KS regularization.
      call cputim(tt1)
          CALL KSTERM
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.3) THEN
* Perform energy check & parameter adjustments and print diagnostics.
          call cputim(tt7)
          CALL ADJUST
          call cputim(tt8)
          ttadj = ttadj + (tt8-tt7)*60.
*
      ELSE IF (IPHASE.EQ.4) THEN
* Switch to unperturbed three-body regularization.
      call cputim(tt1)
          ISUB = 0
          CALL TRIPLE(ISUB)
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.5) THEN
* Switch to unperturbed four-body regularization.
      call cputim(tt1)
          ISUB = 0
          CALL QUAD(ISUB)
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
* Adopt c.m. approximation for inner binary in hierarchical triple.
      ELSE IF (IPHASE.EQ.6) THEN
      call cputim(tt1)
          CALL MERGE
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.7) THEN
* Restore old binary in hierarchical configuration.
      call cputim(tt1)
          CALL RESET
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
* Begin chain regularization.
      ELSE IF (IPHASE.EQ.8) THEN
      call cputim(tt1)
          ISUB = 0
          CALL CHAIN(ISUB)
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
      END IF
*
* Continue integration.
      GO TO 1
*
      END
