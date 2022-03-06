      SUBROUTINE DATA
*
*
*       Initial conditions.
*       -------------------
*
      INCLUDE 'common6.h'
      REAL*4  RAN2
*
*
*       Initialize the portable random number generator (range: 0 to 1).
      KDUM = -1
      RN1 = RAN2(KDUM)
*       Skip the first random numbers (IDUM1 specified at input).
      DO 1 K = 1,IDUM1
          RN1 = RAN2(KDUM)
    1 CONTINUE
*
*       Save random number sequence in COMMON for future use.
      IDUM1 = KDUM
*
*       Read mass function parameters, # primordials, Z-abundance & epoch.
      READ (5,*)  ALPHA, BODY1, BODYN, NBIN0, ZMET, EPOCH0
      IF (ZMET.LE.0.0D0) ZMET = 1.0D-04
*
*       Check option for reading initial conditions from input file.
      IF (KZ(22).GE.2) THEN
          ZMASS = 0.0
          DO 5 I = 1,N
              READ (10,*)  BODY(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3)
              ZMASS = ZMASS + BODY(I)
    5     CONTINUE
*       Include possibility of retaining the mass distribution.
          IF (KZ(22).GT.2) GO TO 30
      END IF
*
*       Include the case of equal masses (ALPHA = 1 or BODY1 = BODYN).
      IF (ALPHA.EQ.1.0.OR.BODY1.EQ.BODYN) THEN
          DO 10 I = 1,N
              BODY(I) = 1.0
   10     CONTINUE
*       Set provisional total mass (rescaled in routine SCALE).
          ZMASS = FLOAT(N)
          GO TO 40
      END IF
*
*       Choose between two realistic IMF's and standard Salpeter function.
      IF (KZ(20).EQ.1) THEN
          CALL IMF(BODY1,BODYN)
          GO TO 30
      ELSE IF (KZ(20).GE.2) THEN
          CALL IMF2(BODY1,BODYN)
          GO TO 30
      END IF
*
      WRITE (6,15)  ALPHA, BODY1, BODYN
   15 FORMAT (/,12X,'STANDARD IMF    ALPHA =',F5.2,
     &              '  BODY1 =',F5.1,'  BODYN =',F5.2)
*
*       Generate a power-law mass function with exponent ALPHA.
      ALPHA1 = ALPHA - 1.0
      FM1 = 1.0/BODY1**ALPHA1
      FMN = (FM1 - 1.0/BODYN**ALPHA1)/(FLOAT(N) - 1.0)
      ZMASS = 0.0D0
      CONST = 1.0/ALPHA1
*
*       Assign individual masses sequentially.
      DO 20 I = 1,N
          FMI = FM1 - FLOAT(I - 1)*FMN
          BODY(I) = 1.0/FMI**CONST
          ZMASS = ZMASS + BODY(I)
   20 CONTINUE
*
*       Scale the masses to <M> = 1 for now and set consistent total mass.
   30 ZMBAR1 = ZMASS/FLOAT(N)
      DO 35 I = 1,N
          BODY(I) = BODY(I)/ZMBAR1
   35 CONTINUE
      ZMASS = FLOAT(N)
*
*       Set up initial coordinates & velocities (uniform or Plummer model).
   40 IF (KZ(22).LE.1) THEN
          CALL SETUP
      END IF
*
   50 RETURN
*
      END
