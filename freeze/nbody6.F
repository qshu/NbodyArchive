	PROGRAM NBODY6
*
*             N B O D Y 6++
*             *************
*
*       Regularized AC N-body code with triple & binary collisions.
*       --------------------------------------------------------
*
*       Hermite integration scheme with block-steps (V 4.0.0 April/99).
*       ------------------------------------------------------------------
*
*       Developed by Sverre Aarseth, IOA, Cambridge.
*       ............................................
*       Message Passing Version NBODY6++ for Massively Parallel Systems
*       Developed by Rainer Spurzem, ARI, Heidelberg
*
      INCLUDE 'common6.h'
      COMMON/STSTAT/  TINIT,NIR,NIB,NRGL,NKS
      EXTERNAL MERGE
*
#ifdef PARALLEL
#define MPIINIT 1
#else
#ifdef ENSEMBLE
#define MPIINIT 1
#else
#define MPIINIT 0
#endif
#endif

#if MPIINIT
*       Initialize MPI
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_GROUP(MPI_COMM_WORLD,group,ierr)
      CALL MPI_GROUP_SIZE(group,isize,ierr)
      CALL MPI_GROUP_RANK(group,rank,ierr)
*     PRINT*,' This is rank=',rank,' size=',isize,' group=',group
#endif
*
*       Initialize the timer.
      CALL CPUTIM(ttota)
*
*       Read start/restart indicator & CPU time.
      IF(rank.eq.0)READ (5,*)  KSTART, TCOMP, TCRITp,
     *    isernb,iserreg
*
#if MPIINIT
      CALL MPI_BCAST(isernb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iserreg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(KSTART,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(TCOMP,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(TCRITp,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
*
      iserreg = max(isize,iserreg)
      isernb = max(isize,isernb)
      IF(rank.eq.0)PRINT*,' iserreg,isernb=',iserreg,isernb
#endif
*
      IF (KSTART.EQ.1) THEN
*
*       Read input parameters, perform initial setup and obtain output.
          CPU = TCOMP
          CALL START
          call cputim(tt7)
          CALL ADJUST
          call cputim(tt8)
          ttadj = ttadj + (tt8-tt7)*60.
      ELSE
*
*       Read previously saved COMMON variables from tape/disc on unit 1.
          CALL MYDUMP(0,1)
*
          IF (NDUMP.GE.3) STOP
*       Safety indicator preventing repeated restarts set in routine CHECK.
          CPU = TCOMP
          CPU0 = 0.0
*       Set IPHASE = -1 for new NLIST in routine INTGRT (Hermite version).
          IPHASE = -1
*
*       Metallicity should be set to minimum value.
      IF(ZMET.EQ.0.D0) ZMET = 1.D-4
*
*       Initialize evolution parameters which depend on metallicity.
          CALL ZCNSTS(ZMET,ZPARS)
*
*       Check reading modified restart parameters (KSTART = 3, 4 or 5).
          IF (KSTART.GT.2) THEN
              CALL MODIFY(KSTART)
          END IF
      END IF
*
*       Initialize the file names for read and write.
      DO 100 IF = 1,9
          WRITE (FILE(IF),112) IF
 100  CONTINUE
 112  FORMAT('fort0000.',I1)
      FILE(11) = 'esc0000.11'
      FILE(15) = 'per0000.15'
#ifdef ENSEMBLE
      DO 101 IF = 1,9
      if (isize.lt.10) then
          WRITE (FILE(IF),113) rank,IF
 113  FORMAT('fort000',I1,'.',I1)
      else if (isize.ge.10.and.isize.lt.100) then
          WRITE (FILE(IF),115) rank,IF
 115  FORMAT('fort00',I2,'.',I1)
      else if (isize.ge.100.and.isize.lt.1000) then
          WRITE (FILE(IF),117) rank,IF
 117  FORMAT('fort0',I3,'.',I1)
      else if (isize.ge.1000) then
          WRITE (FILE(IF),119) rank,IF
 119  FORMAT('fort',I4,'.',I1)
      end if
 101  CONTINUE
*
      if (isize.lt.10) then
          WRITE (FILE(11),213) rank
          WRITE (FILE(15),313) rank
 213  FORMAT('esc000',I1,'.11')
 313  FORMAT('per000',I1,'.15')
      else if (isize.ge.100.and.isize.lt.1000) then
          WRITE (FILE(11),215) rank
          WRITE (FILE(15),315) rank
 215  FORMAT('esc00',I2,'.11')
 315  FORMAT('per00',I2,'.15')
      else if (isize.ge.100.and.isize.lt.1000) then
          WRITE (FILE(11),217) rank
          WRITE (FILE(15),317) rank
 217  FORMAT('esc0',I3,'.11')
 317  FORMAT('per0',I3,'.15')
      else if (isize.ge.1000) then
          WRITE (FILE(11),219) rank
          WRITE (FILE(15),319) rank
 219  FORMAT('esc',I4,'.11')
 319  FORMAT('per',I4,'.15')
      end if
#endif
#ifdef PARALLEL
      if(rank.eq.0) then
#endif
      IF (KZ(1).GT.0)
     &OPEN (UNIT=1,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(1))
      IF (KZ(2).GT.0)
     &OPEN (UNIT=2,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(2))
      IF (KZ(3).GT.0) 
     &OPEN (UNIT=3,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(3))
      IF (BK(3).GT.0)
     &OPEN (UNIT=4,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(4))
      IF (KZ(7).GE.3)
     &OPEN (UNIT=7,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(7))
      IF (KZ(8).GE.2)
     &OPEN (UNIT=9,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(9))
      IF (KZ(23).EQ.2.OR.KZ(23).EQ.4)
     &OPEN (UNIT=11,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(11))
      IF (BK(4).EQ.1)
     &OPEN (UNIT=15,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(15))
#ifdef PARALLEL
      end if
#endif


* (R.Sp.)Set time flag and step number flags for beginning of run
      TINIT = TTOT
      NIR = NSTEPI
      NIB = NSTEPB
      NRGL = NSTEPR
      NKS = NSTEPU
*
      call cputim(tt2)
      ttinit = ttinit + (tt2-ttota)*60.
*       Advance solutions until next output or change of procedure.
    1 CONTINUE
      call cputim(tt1)
*
      CALL INTGRT
*
      call cputim(tt2)
      ttint = ttint + (tt2-tt1)*60.
*
      IF (IPHASE.EQ.1) THEN
*       Prepare new KS regularization.
      call cputim(tt1)
          CALL KSREG
          CALL FLUSH(6)
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.2) THEN
*       Terminate KS regularization.
      call cputim(tt1)
          CALL KSTERM
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.3) THEN
*       Perform energy check & parameter adjustments and print diagnostics.
          call cputim(tt7)
          CALL ADJUST
          call cputim(tt8)
          ttadj = ttadj + (tt8-tt7)*60.
*
      ELSE IF (IPHASE.EQ.4) THEN
*       Switch to unperturbed three-body regularization.
      call cputim(tt1)
          ISUB = 0 
          CALL TRIPLE(ISUB)
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.5) THEN
*       Switch to unperturbed four-body regularization.
      call cputim(tt1)
          ISUB = 0
          CALL QUAD(ISUB)
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
*       Adopt c.m. approximation for inner binary in hierarchical triple.
      ELSE IF (IPHASE.EQ.6) THEN
      call cputim(tt1)
          CALL MERGE
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.7) THEN
*       Restore old binary in hierarchical configuration.
      call cputim(tt1)
          CALL RESET
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
*       Begin chain regularization.
      ELSE IF (IPHASE.EQ.8) THEN
      call cputim(tt1)
          ISUB = 0
          CALL CHAIN(ISUB)
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
      END IF
*
*       Continue integration.
      GO TO 1
*
      END
