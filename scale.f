


      SUBROUTINE mySCALE
*
*
*       Scaling to new units.
*       ---------------------
*
      INCLUDE 'common1.h'
      integer i,k
      REAL*4  CMR(3),CMRDOT(3)
      real * 8 q, vxrot, vzrot, zkin, pot, qv, e0, etot, sx, omega
*
*
*       Read virial ratio and rotation scaling factors.
      READ (5,*)  Q, VXROT, VZROT
*
      ZMASS = 0.0
      DO 10 K = 1,3
          CMR(K) = 0.0
          CMRDOT(K) = 0.0
   10 CONTINUE
*
*       Form centre of mass displacements.
      DO 30 I = 1,N
          ZMASS = ZMASS + BODY(I)
          DO 25 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   25     CONTINUE
   30 CONTINUE
*
*       Adjust coordinates and velocities to c.m. rest frame.
      DO 40 I = 1,N
          DO 35 K = 1,3
              X(K,I) = X(K,I) - CMR(K)/ZMASS
              XDOT(K,I) = XDOT(K,I) - CMRDOT(K)/ZMASS
   35     CONTINUE
   40 CONTINUE
*
*       Scale masses to standard units of <M> = 1/N.
      DO 50 I = 1,N
          BODY(I) = BODY(I)/ZMASS
   50 CONTINUE
      ZMASS = 1.0
*
*       Obtain the total kinetic & potential energy.
      CALL ENERGY(ZKIN,POT)
      write(6,*) 'Energy retvals:', zkin, pot
*
*       Scale non-zero velocities by virial theorem ratio.
      IF (ZKIN.GT.0.0) THEN
          QV = SQRT(Q*POT/ZKIN)
          DO 55 I = 1,N
              DO 52 K = 1,3
                  XDOT(K,I) = XDOT(K,I)*QV
   52         CONTINUE
   55     CONTINUE
      END IF
*
*       Scale total energy to standard units (E = -0.25).
      E0 = -0.25
      ETOT = (Q - 1.0)*POT
      SX = E0/ETOT
*
      WRITE (6,65)  SX, ETOT, BODY(1), BODY(N), ZMASS/FLOAT(N)
   65 FORMAT (/,12X,'SCALING:   SX  =',F6.2,'  E =',1PE10.2,
     &                   '  M(1) =',E9.2,'  M(N) =',E9.2,'  <M> =',E9.2)
*
*       Scale coordinates & velocities to the new units.
      DO 70 I = 1,N
          DO 68 K = 1,3
              X(K,I) = X(K,I)/SX
              XDOT(K,I) = XDOT(K,I)*SQRT(SX)
   68     CONTINUE
   70 CONTINUE
*
*       Check whether to include rotation (VXROT = 0 in standard case). 
      IF (VXROT.GT.0.0) THEN
*
*       Set angular velocity for retrograde motion (i.e. star clusters).
          OMEGA = -SX*SQRT(ZMASS*SX)
          WRITE (6,75)  VXROT, VZROT, OMEGA
   75     FORMAT (/,12X,'VXROT =',F6.2,'  VZROT =',F6.2,
     &                                                 '  OMEGA =',F7.2)
*
*       Add solid-body rotation about Z-axis (reduce random velocities).
          DO 80 I = 1,N
              XDOT(1,I) = XDOT(1,I)*VXROT - X(2,I)*OMEGA
              XDOT(2,I) = XDOT(2,I)*VXROT + X(1,I)*OMEGA
              XDOT(3,I) = XDOT(3,I)*VZROT
   80     CONTINUE
      END IF
*
*       Set initial crossing time in scaled units.
      TCR = ZMASS**2.5/(2.0*ABS(E0))**1.5
*
*       Scale output time interval & termination time by initial TCR.
      IF (KZ(10).EQ.0) THEN
          DELTAT = DELTAT*TCR
          TCRIT = TCRIT*TCR
      END IF
*
      RETURN
*
      END
