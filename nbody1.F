*
*             N B O D Y 1 H
*             *************
*
*
*       Standard N-body code with Hermite scheme.
*       -----------------------------------------
*
*       Developed by Jun Makino & Sverre Aarseth.
*       .........................................
*
*
      PROGRAM NBODY1
*
      INCLUDE 'common1.h'
      INCLUDE 'harp3.h'
      integer kstart, j, newval, imode
      real*8 tcomp
*
*
*      memory locking routine
*      call lock_memory
*       Initialize the timer.
#ifdef HARP3
      call set_debug_level(-1)
#endif
#ifdef MONOPOLIZEHARP
      call h3open
#endif
      CALL CPUTIM(TCOMP) 
*
#ifdef GRAPE6      
c IF GRAPE-6, first read number of clusters
c      READ (5,*)  g6ncluster
c      write(6,*) 'G6NCLUSTER = ', g6ncluster
      g6ncluster = 1
#endif      
*       Read start/restart indicator & CPU time.
      READ (5,*)  KSTART, TCOMP
      write(6,*)'Kstart, tcomp = ',kstart, tcomp
      IF (KSTART.EQ.1) GO TO 2
*
*       Read previously saved COMMON variables from tape/disc on unit 1.
      CALL MYDUMP(0,1)
      CPU = TCOMP
      IF (KSTART.EQ.2) GO TO 4
*
*       Read modified parameters.
      READ (5,*)  ETA, j, newval
      kz(j) = newval
      WRITE (6,1)  ETA, QE, TCRIT, j, kz(j)
    1 FORMAT (/,7X,'RESTART PARAMETERS:   ETA =',F7.3,'  QE =',F8.5,
     &     '  TCRIT =',F6.1,
     &     ' kz(', i2, ')= ', i5)
      GO TO 4
*
    2 CPU = TCOMP
*
*       Read input parameters and perform initial setup.
      CALL START
*
*       Calculate total energy and produce output.
    3 CALL OUTPUT(imode)
      call flush(6)
*
*       Check the elapsed CPU time or termination time.
      CALL CPUTIM(TCOMP)
      IF (TCOMP.GT.CPU.OR.TIME.GT.TCRIT) GO TO 6
*
    4 continue
*      memory locking routine
*      call lock_memory
      CALL INTGRT
*      memory locking routine
*      call unlock_memory
      GO TO 3
*
*       Check optional COMMON save for subsequent restart.
    6 IF (KZ(1).GT.0) CALL MYDUMP(1,1)
      WRITE (6,8)  TIME, TCOMP, KZ(1)
    8 FORMAT (/,9X,'TIME =',F7.2,'  TCOMP =',F8.2,'  KZ(1) =',I2)
*
      STOP
*
      END
