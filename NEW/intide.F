      SUBROUTINE INTIDE
*
*
*       Input & scaling for tidal dissipation.
*       --------------------------------------
*
      INCLUDE 'common6.h'
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
*       Read parameters for tidal capture simulation.

      if(rank.eq.0)READ (5,*)  RSTAR, IMS, IEV, RMS, REV
#if MPIINIT
      CALL MPI_BCAST(RSTAR,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(IMS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(IEV,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RMS,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(REV,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif
*
*       Convert radii from S.U. to internal units.
      RSTAR = RSTAR/SU
*
      WRITE (6,10)  RSTAR, IMS, IEV, RMS, REV
   10 FORMAT (/,12X,'OLD TIDAL:    RSTAR =',1PE8.1,'  IMS =',0P,I5,
     &              '  IEV =',I4,'  RMS/R* =',F6.2,'  REV/R* =',F6.1)
*
*       Assign individual radii for main-sequence and evolved stars.
      DO 30 I = 1,N
*       Adopt a primitive scheme in case of no stellar evolution.
          IF (I.LE.IMS) THEN
              RADIUS(I) = RMS*RSTAR
          ELSE
              RADIUS(I) = REV*RSTAR
          END IF
   30 CONTINUE
*
*     WRITE (6,40)  BODY(1), BODY(N), RADIUS(1), RADIUS(N)
*  40 FORMAT (/,12X,'SCALED RADII:    M(1) =',F8.4,'  M(N) =',F8.4,
*    &                                '  R(1) =',1PE8.1,'  R(N) =',E8.1)
*
      RETURN
*
      END
