      SUBROUTINE DATA
*
*
*       Initial conditions.
*       -------------------
*
      INCLUDE 'common4.h'
      REAL*4  RAN2
      CHARACTER*1 CHAR(80)
      LOGICAL LREADN,LREADP,LREADF
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
*       Set provisional total mass (rescaled in routine SCALE).
      ZMASS = FLOAT(N)
*
*       Read mass function parameters, # primordials, Z-abundance, epoch
*       and interval for calling routine HRPLOT.
      READ (5,*)  ALPHA, BODY1, BODYN, NBIN0, ZMET, EPOCH0, DTPLOT
*
*       Check option for reading initial conditions from input file.
      IF (KZ(22).GE.2) THEN
          ZMASS = 0.0D0
          DO 5 I = 1,N
              READ (10,*)  BODY(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3)
              ZMASS = ZMASS + BODY(I)
    5     CONTINUE
*       End reading NBODY input data format.
*
*        Read TREE input format
      IF (KZ(22).EQ.3) THEN
          READ(10,*) N
          READ(10,*) DUMDY
          READ(10,*) DUMDY
          PRINT*,' N=',N
          DO 51 I = 1,N
          READ(10,*)BODY(I)
   51     CONTINUE
          PRINT*,' masses read ',BODY(1),BODY(N)
          DO 52 I = 1,N
   52     READ(10,*)(X(K,I),K=1,3)
          DO 53 I = 1,N
   53     READ(10,*)(XDOT(K,I),K=1,3)
          NTOT = N
          PRINT*,' rank ',rank,N,' body data read from unit 10 '
          CALL FLUSH(6)
      END IF
*      End read TREE input format
*        Read STARLAB input format
      IF (KZ(22).EQ.4) THEN
          I = 0
          IS = 0
          LREADF = .FALSE.
          LREADP = .FALSE.
  61      CONTINUE
          READ(10,'(2A1)')(CHAR(K),K=1,2)
*
          LREADN=(.NOT.LREADF).AND.CHAR(1).EQ.'('.AND.CHAR(2).EQ.'P'
          IF(LREADN)THEN
          LREADF=.TRUE.
          READ(10,111)N
          PRINT*,' Read N=',N
          NTOT = N
          END IF
*
          LREADP=CHAR(1).EQ.'('.AND.CHAR(2).EQ.'D'
          IF(LREADP.AND.IS.EQ.0)THEN
          IS = 1
          ELSE
          IF(LREADP)THEN
          I = I + 1
          READ(10,*)CHAR(1),CHAR(2),BODY(I)
          READ(10,*)CHAR(1),CHAR(2),(X(K,I),K=1,3)
          READ(10,*)CHAR(1),CHAR(2),(XDOT(K,I),K=1,3)
          END IF
          END IF
          IF(I.LT.N)GO TO 61
          PRINT*,N,' Particles read from Starlab File'
 111  FORMAT(5X,I5)
      END IF
*
          IF (KZ(22).GE.2) GO TO 50
      END IF
*
*       Include the case of equal masses (ALPHA = 1 or BODY1 = BODYN).
      IF (ALPHA.EQ.1.0.OR.BODY1.EQ.BODYN) THEN
          DO 10 I = 1,N
              BODY(I) = 1.0
   10     CONTINUE
          GO TO 40
      END IF
*
      CALL IMF(BODY1,BODYN)
*
*       Scale the masses to <M> = 1 for now and set consistent total mass.
   30 ZMBAR = ZMASS/FLOAT(N)
      DO 35 I = 1,N
          BODY(I) = BODY(I)/ZMBAR
   35 CONTINUE
      ZMTOT = ZMASS
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
