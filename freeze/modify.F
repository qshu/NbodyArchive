      SUBROUTINE MODIFY(KSTART)
*
*
*       Parameter modification at restart.
*       ----------------------------------
*
      INCLUDE 'common6.h'
      EXTERNAL VERIFY
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
*
*       Read first, second or both lines (KSTART = 3, 4, 5).
      IF (KSTART.EQ.4) GO TO 10
*
*       Read new DTADJ, DELTAT, TADJ, TNEXT, TCRIT, QE & KZ(J) (if > 0).
      if(rank.eq.0)then
      READ (5,*)  DTA, DT, TA, TN, TC, QE1, J, K
      end if
*
#if MPIINIT
      CALL MPI_BCAST(DTA,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(DT,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(TA,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(TN,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(TC,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(QE1,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(J,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(K,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
*
*       Choose appropriate scaling factor.
      IF (KZ(35).EQ.0) THEN
          TCR1 = TCR0
      ELSE
          TCR1 = 1.0
      END IF
*
*       Set new parameters if corresponding input is non-zero.
      IF (DTA.LE.0.0) THEN
          DTA = DTADJ
      ELSE
          DTA = DTA*TCR1
      END IF
*
      IF (DT.LE.0.0) THEN
          DT = DELTAT
      ELSE
          DT = DT*TCR1
      END IF
*
      IF (TA.LE.0.0) THEN
          TADJ = MAX(TADJ - DTADJ + DTA,TIME)
      ELSE
          TADJ = MAX(TA*TCR1,TIME)
      END IF
*
      IF (TN.LE.0.0) THEN
          TNEXT = MAX(TNEXT - DELTAT + DT,TIME)
      ELSE
          TNEXT = MAX(TN*TCR1,TIME)
      END IF
*
      DTADJ = DTA
      DELTAT = DT
      IF (TC.GT.0.0) TCRIT = TIME + TC*TCR1
      IF (QE1.GT.0.0) QE = QE1
*
*       See whether any options should be changed.
      IF (J.GT.0) KZ(J) = K
*
      if(rank.eq.0)then
      WRITE (6,5)  DTADJ, DELTAT, TCRIT, QE, J, K
    5 FORMAT (///,7X,'RESTART PARAMETERS:   DTADJ =',F7.3,'  DELTAT =',
     &                            F7.3,'  TCRIT =',F7.1,'  QE =',1PE9.1,
     &                                            '  KZ(',I2,') =',I2,/)
      WRITE(6,50) TIME/TCR0, TCRIT/TCR0, TIME
   50 FORMAT (' START AT',F12.4,' REG. STOP INTENDED AT ',F12.4,' TCR0',
     &       /,' START AT TIME ',F16.8)
      end if
*
*       Read new ETAI, ETAR, ETAU, DTMIN, RMIN (IF > 0 & KSTART = 4 or 5).
   10 IF (KSTART.GE.4) THEN
          if(rank.eq.0)then
          READ (5,*)  ETA1, ETA2, ETA3, DTM, RM, NNBO
          end if
#if MPIINIT
      CALL MPI_BCAST(ETA1,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ETA2,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ETA3,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(DTM,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RM,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif
*
*       Check modification of integration parameters.
          IF (ETA1.GT.0.0) ETAI = ETA1
          IF (ETA2.GT.0.0) ETAR = ETA2
          IF (ETA3.GT.0.0) ETAU = ETA3
          IF (DTM.GT.0.0) THEN
              DTMIN = DTM
              SMIN = 2.0*DTM
          END IF
          IF (RM.GT.0.0) THEN
              RMIN = RM
              RMIN2 = RM**2
              RMIN22 = 4.0*RMIN2
          END IF
          IF (NNBO.GT.0) NNBOPT = NNBO
*
          if(rank.eq.0)then
          WRITE (6,15)  ETAI, ETAR, ETAU, DTMIN, RMIN, NNBOPT
   15     FORMAT (/,7X,'RESTART PARAMETERS:   ETAI =',F7.3,'  ETAR =',
     &                          F7.3,'  ETAU =',F7.3,'  DTMIN =',1PE9.1,
     &                               '  RMIN =',E9.1,' NNBOPT =',I5,/)
          end if
      END IF
*
*       Perform a simple validation check on main input parameters.
      CALL VERIFY
*
*       Save the new parameters on tape/disc in case a restart is needed.
      CALL MYDUMP(1,1)
*
      RETURN
*
      END
