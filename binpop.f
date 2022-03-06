      SUBROUTINE BINPOP
*
*
*       Initial binary population.
*       --------------------------
*
      INCLUDE 'common6.h'
      REAL*8  XORB(2),VORB(2),XREL(3),VREL(3),PX(3),QX(3),BS(NMAX)
      REAL  RAN2
*
*
      READ (5,*)  NBIN, SEMI, ECC0, RATIO, NBGR, REDUCE, RANGE, NSKIP,
     &            IDORM
      NBIN0 = NBIN
      NBIN1 = NBIN + 1
      SEMI0 = SEMI
      WRITE (6,1)  NBIN, SEMI, ECC0, RATIO, NBGR, REDUCE, RANGE, NSKIP,
     &             IDORM
    1 FORMAT (/,7X,'BINARIES:   NBIN =',I4,'  A =',F9.6,'  E =',F6.2,
     &               '  RATIO =',F4.1,'  NBGR =',I3,'  REDUCE =',F5.2,
     &               '  RANGE =',F6.1,'  NSKIP =',I3,'  IDORM =',I2,/)
*
*       Check type of mass splitting & sampling procedure.
      IF (NSKIP.EQ.0.OR.KZ(20).GE.2) GO TO 10
      IF (RATIO.EQ.1.0) GO TO 20
*
*       Select binaries from the most massive bodies (frequency NSKIP).
      ILAST = (1 + NSKIP)*NBIN
      JSKIP = 0
      JS = 0
      JB = 1
*
*       Transfer binary masses to first NBIN locations.
      DO 6 I = 2,ILAST
          JSKIP = JSKIP + 1
*       Copy binary mass of body #I to new global location.
          IF (JSKIP.GT.NSKIP) THEN
              JSKIP = 0
              JB = JB + 1
              BODY(JB) = BODY(I)
          ELSE
*       Save next NSKIP masses of single bodies.
              JS = JS + 1
              BS(JS) = BODY(I)
          END IF
    6 CONTINUE
*
*       Restore the single bodies in subsequent locations.
      JS = 0
      DO 8 I = NBIN1,ILAST
          JS = JS + 1
          BODY(I) = BS(JS)
    8 CONTINUE
*
*       Move main variables of all single bodies.
   10 DO 15 I = N,NBIN1,-1
          J = I + NBIN
          BODY(J) = BODY(I)
          DO 12 K = 1,3
              X(K,J) = X(K,I)
              XDOT(K,J) = XDOT(K,I)
   12     CONTINUE
   15 CONTINUE
*
*       Create space for each binary component next to primary.
   20 DO 30 I = NBIN,2,-1
          J = 2*I - 1
          BODY(J) = BODY(I)
          DO 25 K = 1,3
              X(K,J) = X(K,I)
              XDOT(K,J) = XDOT(K,I)
   25     CONTINUE
   30 CONTINUE
*
*       Specify binary components from relative motion.
      IBIN = 0
      DO 60 I = 1,NBIN
*       Randomize perihelion, node & inclination.
          PI = TWOPI*RAN2(IDUM1)
          OMEGA = TWOPI*RAN2(IDUM1)
          ZI = 0.25*TWOPI*RAN2(IDUM1)
*
*       Set transformation elements (Brouwer & Clemence p. 35).
          PX(1) = COS(PI)*COS(OMEGA) - SIN(PI)*SIN(OMEGA)*COS(ZI)
          QX(1) =-SIN(PI)*COS(OMEGA) - COS(PI)*SIN(OMEGA)*COS(ZI)
          PX(2) = COS(PI)*SIN(OMEGA) + SIN(PI)*COS(OMEGA)*COS(ZI)
          QX(2) =-SIN(PI)*SIN(OMEGA) + COS(PI)*COS(OMEGA)*COS(ZI)
          PX(3) = SIN(PI)*SIN(ZI)
          QX(3) = COS(PI)*SIN(ZI)
*
*       Form component masses.
          I1 = 2*I - 1
          I2 = 2*I
          ZMBIN = BODY(I1)
*
*       Specify component masses (copy BODY0 from IMF2 or use RATIO).
          IF (KZ(20).GE.2) THEN
              BODY(I1) = BODY0(I1)
              BODY(I2) = BODY0(I2)
              ZMBIN = BODY(I1) + BODY(I2)
          ELSE IF (RATIO.EQ.1.0) THEN
              BODY(I2) = BODY(I1)
              ZMBIN = ZMBIN + BODY(I2)
          ELSE
              BODY(I1) = RATIO*ZMBIN
              BODY(I2) = (1.0 - RATIO)*ZMBIN
          END IF
*
*       Choose random (thermalized) or fixed eccentricity.
          IF (ECC0.LT.0.0) THEN
              ECC2 = RAN2(IDUM1)
              ECC = SQRT(ECC2)
          ELSE
              ECC = ECC0
          END IF
*
*       Select semi-major axis from uniform distribution in log(A).
          IF (RANGE.GT.0.0) THEN
              EXP = RAN2(IDUM1)*LOG10(RANGE)
              SEMI = SEMI0/10.0**EXP
          END IF
*
*       Specify relative motion at apocentre and sum binding energy.
          XORB(1) = SEMI*(1.0 + ECC)
          XORB(2) = 0.0
          VORB(1) = 0.0
          VORB(2) = SQRT(ZMBIN*(1.0D0 - ECC)/(SEMI*(1.0D0 + ECC)))
          EBIN0 = EBIN0 - 0.5*BODY(I1)*BODY(I2)/SEMI
*
*       Transform to relative variables.
          DO 40 K = 1,3
              XREL(K) = PX(K)*XORB(1) + QX(K)*XORB(2)
              VREL(K) = PX(K)*VORB(1) + QX(K)*VORB(2)
   40     CONTINUE
*
*       Set global variables for each component.
          DO 50 K = 1,3
              X(K,I1) = X(K,I1) + BODY(I2)*XREL(K)/ZMBIN
              X(K,I2) = X(K,I1) - XREL(K)
              XDOT(K,I1) = XDOT(K,I1) + BODY(I2)*VREL(K)/ZMBIN
              XDOT(K,I2) = XDOT(K,I1) - VREL(K)
   50     CONTINUE
*
          IBIN = IBIN + 1
          IF (IBIN.GE.NBGR) THEN
              SEMI = REDUCE*SEMI
              IBIN = 0
          END IF
   60 CONTINUE
*
*       Update the total particle number after splitting primaries.
      IF (RATIO.LT.1.0) THEN
          N = N + NBIN
          NZERO = N
          NTOT = N
          IF (NSKIP.GT.0) THEN
              WRITE (6,62)  (BODY(J),J=1,10)
              WRITE (6,64)  (BODY(J),J=2*NBIN+1,2*NBIN+10)
   62         FORMAT (/,7X,'BINARY MASSES (1-10):  ',10F9.5)
   64         FORMAT (/,7X,'SINGLE MASSES (1-10):  ',10F9.5,/)
          END IF
      END IF
*
*       Check procedure for introducing dormant binaries.
      IF (IDORM.GT.0) THEN
          DO 66 I = 1,NBIN
              I1 = 2*I - 1
              I2 = I1 + 1
              ZM = BODY(I1) + BODY(I2)
              DO 65 K = 1,3
              X(K,I) = (BODY(I1)*X(K,I1) + BODY(I2)*X(K,I2))/ZM
              XDOT(K,I) = (BODY(I1)*XDOT(K,I1) + BODY(I2)*XDOT(K,I2))/ZM
   65         CONTINUE
              BODY(I) = ZM
   66     CONTINUE
          I1 = 2*NBIN + 1
          I2 = NBIN
          DO 68 I = I1,N
              I2 = I2 + 1
              BODY(I2) = BODY(I)
              DO 67 K = 1,3
                  X(K,I2) = X(K,I)
                  XDOT(K,I2) = XDOT(K,I)
   67         CONTINUE
   68     CONTINUE
*
*       Reset particle membership and turn off binary output option.
*         N = N - NBIN
          NZERO = N
          NTOT = N
          NBIN0 = 0
          EBIN0 = 0.0
          IF (KZ(8).EQ.1) KZ(8) = 0
      END IF
*
*       Set coordinates & velocities in c.m. rest frame.
      DO 70 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   70 CONTINUE
*
      DO 80 I = 1,N
          DO 75 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   75     CONTINUE
   80 CONTINUE
*
      DO 90 I = 1,N
          DO 85 K = 1,3
              X(K,I) = X(K,I) - CMR(K)/ZMASS
              XDOT(K,I) = XDOT(K,I) - CMRDOT(K)/ZMASS
   85     CONTINUE
   90 CONTINUE
*
  100 RETURN
*
      END
