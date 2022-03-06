      SUBROUTINE STEPS(I1,I2,KCASE)
*
*
*       Initialization of time-steps & prediction variables.
*       ----------------------------------------------------
*
      INCLUDE 'common1.h'
      integer i,i1,i2,kcase,iter,k
      real * 8 fi2, fd2, dt, dtn
*
*
*       Set new steps and initialize prediction variables.
      DO 40 I = I1,I2
*
*       Obtain time-step using DT = ETA*F/FDOT (D2 & D3 not available).
          FI2 = F(1,I)**2 + F(2,I)**2 + F(3,I)**2
          FD2 = FDOT(1,I)**2 + FDOT(2,I)**2 + FDOT(3,I)**2
          DT = ETA*SQRT(FI2/FD2)
           dt=dt*0.1
          if(dt .gt. deltat) dt = deltat
*
*       Initialize the time and obtain discrete time-step if required.
          T0(I) = TIME
*
*       Convert predicted step to nearest block time-step (truncated down).
          CALL STEPK(DT,DTN)
          IF (TIME.LE.0.0D0) THEN
              STEP(I) = DTN
          ELSE 
*       Reduce step by factor 2 until commensurate with current time.
              STEP(I) = DTN
              ITER = 0
   10         IF (DMOD(TIME,STEP(I)).NE.0.0D0) THEN
                  STEP(I) = 0.5D0*STEP(I)
                  ITER = ITER + 1
                  IF (ITER.LT.10) GO TO 10
                  WRITE (6,15)  I, ITER, TIME/STEP(I), DT, STEP(I)
   15             FORMAT (' WARNING!   I ITER T/STEP DT STEP ',
     &                                 I5,I4,F16.4,1P,2E9.1)
                  IF (STEP(I).GT.DTK(32)) GO TO 10
                  STOP
              END IF
          END IF
*
*       Update array for new block times (#I copied to free location).
          TIMENW(I) = T0(I) + STEP(I)
*         IF (TIME.GT.0.0) THEN
*             WRITE (7,20)  I, NAME(I), TIME, DT, STEP(I)
*  20         FORMAT (' STEPS:   I NM TIME DT STEP ',2I5,F12.6,1P,2E9.1)
*             CALL FLUSH(7)
*         END IF
*
*       Set prediction variables.
          DO 30 K = 1,3
              X0(K,I) = X(K,I)
              X0DOT(K,I) = XDOT(K,I)
              F(K,I) = 0.5D0*F(K,I)
              FDOT(K,I) = ONE6*FDOT(K,I)
   30     CONTINUE
   40 CONTINUE
*
      RETURN
*
      END
