      SUBROUTINE SCALE
*
*
*       Scaling to new units.
*       ---------------------
*
      INCLUDE 'common6.h'
*
*
*       Read virial ratio, rotation scaling factors & boundary radius.
      READ (5,*)  Q, VXROT, VZROT, RSPH2
*
      ZMASS = 0.0D0
      DO 10 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   10 CONTINUE
*
*       Form total mass and centre of mass displacements.
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
*       Scale masses to standard units of <M> = 1/N and set total mass.
      DO 50 I = 1,N
          BODY(I) = BODY(I)/ZMASS
   50 CONTINUE
      ZMASS = 1.0
*
*       Obtain the total kinetic & potential energy.
      CALL ENERGY
*
*       Use generalized virial theorem for external tidal field.
      IF (KZ(14).EQ.1) THEN
          AZ = 0.0D0
          DO 55 I = 1,N
              AZ = AZ + BODY(I)*(X(1,I)*XDOT(2,I) - X(2,I)*XDOT(1,I))
   55     CONTINUE
          ZKIN1 = ZKIN + 0.5*TIDAL(4)*AZ
          VIR = POT - 2.0*ETIDE
      ELSE
          ZKIN1 = ZKIN
          VIR = POT
      END IF
*
*       Scale non-zero velocities by virial theorem ratio.
      IF (ZKIN.GT.0.0D0) THEN
          QV = SQRT(Q*VIR/ZKIN1)
          DO 60 I = 1,N
              DO 58 K = 1,3
                  XDOT(K,I) = XDOT(K,I)*QV
   58         CONTINUE
   60     CONTINUE
      END IF
*
*       Scale total energy to standard units (E = -0.25 for Q < 1).
      E0 = -0.25
      ETOT = (Q - 1.0)*POT
      SX = E0/ETOT
*       Include case of hot system inside reflecting boundary.
      IF (KZ(29).GT.0.AND.Q.GT.1.0) THEN
          E0 = ETOT
      END IF
*       Define scaling factor (set E0 = ETOT if energy scaling not desired).
      SX = E0/ETOT
*
      WRITE (6,65)  SX, ETOT, BODY(1), BODY(N), ZMASS/FLOAT(N)
   65 FORMAT (//,12X,'SCALING:   SX =',F6.2,'  E =',1PE10.2,
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
      IF (VXROT.GT.0.0D0) THEN
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
*       Check option for writing the initial conditions on unit 10.
      IF (KZ(22).EQ.1) THEN
          DO 85 I = 1,N
              WRITE (10,84)  BODY(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3)
   84         FORMAT (1P,7E14.6)
   85     CONTINUE
      END IF
*
*       Check option for reading initial subsystems.
      IF (KZ(24).GT.0) THEN
          K = KZ(24)
          DO 90 I = 1,K
              READ (5,*)  (X(J,I),J=1,3), (XDOT(J,I),J=1,3)
   90     CONTINUE
      END IF
*
*       Set initial crossing time in scaled units.
      TCR = ZMASS**2.5/(2.0D0*ABS(E0))**1.5
      TCR0 = TCR
*
*       Obtain approximate half-mass radius after scaling.
      RSCALE = 0.5*ZMASS**2/(SX*POT)
*       Set square radius of reflecting sphere.
      RSPH2 = (RSPH2*RSCALE)**2
*       Form equilibrium rms velocity (temporarily defined as VC).
      VC = SQRT(2.0D0*ABS(E0)/ZMASS)
*
*       Check for general binary search of initial condition.
      IF (KZ(4).GT.0) THEN
          CALL EVOLVE(0,0)
      END IF
*
*       Print half-mass relaxation time & equilibrium crossing time.
      A1 = FLOAT(N)
      TRH = 4.0*TWOPI/3.0*(VC*RSCALE)**3/(15.4*ZMASS**2*LOG(A1)/A1)
      WRITE (6,95)  TRH, TCR, 2.0*RSCALE/VC
   95 FORMAT (/,12X,'TIME SCALES:   TRH =',1PE8.1,'  TCR =',E8.1,
     &                                            '  2<R>/<V> =',E8.1,/)
*
      RETURN
*
      END
