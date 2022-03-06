# 1 "data.F"
# 1 "<built-in>"
# 1 "<command line>"
# 1 "data.F"
      SUBROUTINE DATA
*
*
* Initial conditions.
* -------------------
*
      INCLUDE 'common6.h'
      REAL*8 A(8)
      CHARACTER*1 CHAR(80)
      LOGICAL LREADN,LREADP,LREADF
*
# 21 "data.F"
*
* Set provisional total mass (rescaled in routine SCALE).
      ZMASS = FLOAT(N)
*
* Check option for reading initial conditions from input file.
      if(rank.eq.0)then
*
      IF (KZ(22).EQ.2.OR.KZ(22).EQ.6) THEN
          DO 5 I = 1,N
              READ (10,*) BODY(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3)
    5 CONTINUE
* Read tidal radius if cutoff required
      IF (KZ(23).GE.3) READ (10,*) RTIDE
          PRINT*,' rank ',rank,N,' body data read from unit 10 ',
     * ' RTIDE =',RTIDE
      END IF
* End reading NBODY input data format.
*
* Read TREE input format
      IF (KZ(22).EQ.3.OR.KZ(22).EQ.7) THEN
          READ(10,*) N
          READ(10,*) DUMDY
          READ(10,*) DUMDY
          PRINT*,' N=',N
          DO 51 I = 1,N
          READ(10,*)BODY(I)
   51 CONTINUE
          PRINT*,' masses read ',BODY(1),BODY(N)
          DO 52 I = 1,N
   52 READ(10,*)(X(K,I),K=1,3)
          DO 53 I = 1,N
   53 READ(10,*)(XDOT(K,I),K=1,3)
          NTOT = N
          PRINT*,' rank ',rank,N,' body data read from unit 10 '
          CALL FLUSH(6)
      END IF
* End read TREE input format
* Read STARLAB input format
      IF (KZ(22).EQ.4.OR.KZ(22).EQ.8) THEN
          I = 0
          IS = 0
          LREADF = .FALSE.
          LREADP = .FALSE.
  61 CONTINUE
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
 111 FORMAT(5X,I5)
      END IF
*
      end if
# 101 "data.F"
*
* Read mass function parameters, # primordials, Z-abundance & epoch.
* And plot interval for HR diagram.
       if(rank.eq.0)then
       READ (5,*) ALPHA, BODY1, BODYN, NBIN0, ZMET, EPOCH0,DTPLOT
       end if
# 116 "data.F"
*
      IF(KZ(22).GE.2)GO TO 50
*
* Include the case of equal masses (ALPHA = 1 or BODY1 = BODYN).
      IF (ALPHA.EQ.1.0.OR.BODY1.EQ.BODYN) THEN
          DO 10 I = 1,N
              BODY(I) = 1.0
   10 CONTINUE
          IF (KZ(20).NE.6) GO TO 40
      END IF
*
* Choose between two realistic IMFs and standard Salpeter function.
      IF (KZ(20).EQ.1) THEN
          CALL IMF(BODY1,BODYN)
          GO TO 30
      ELSE IF (KZ(20).GE.2) THEN
          CALL IMF2(BODY1,BODYN)
          GO TO 30
      END IF
* Metallicity should be set to minimum value.
      IF(ZMET.EQ.0.D0) ZMET = 1.D-4
*
* Generate a power-law mass function with exponent ALPHA.
      ALPHA1 = ALPHA - 1.0
      FM1 = 1.0/BODY1**ALPHA1
      FMN = (FM1 - 1.0/BODYN**ALPHA1)/(FLOAT(N) - 1.0)
      ZMASS = 0.0D0
      CONST = 1.0/ALPHA1
*
* Assign individual masses sequentially.
      DO 20 I = 1,N
          FMI = FM1 - FLOAT(I - 1)*FMN
          BODY(I) = 1.0/FMI**CONST
          ZMASS = ZMASS + BODY(I)
   20 CONTINUE
*
      if(rank.eq.0)then
      WRITE (6,15) ALPHA, BODY1, BODYN, ZMASS, NBIN0, ZMET, EPOCH0
   15 FORMAT (/,12X,'STANDARD IMF    ALPHA =',F5.2,
     & '  BODY1 =',F5.1,'  BODYN =',F5.2,' ZMASS =',1P,E12.5,0P,
     & ' NBIN0=',I5,' ZMET =',F5.2,' EPOCH0 =',F5.2)
      end if
*
* Scale the masses to <M> = 1 for now and set consistent total mass.
   30 ZMBAR1 = ZMASS/FLOAT(N)
      DO 35 I = 1,N
          BODY(I) = BODY(I)/ZMBAR1
   35 CONTINUE
      ZMASS = FLOAT(N)
*
* Set up initial coordinates & velocities (uniform or Plummer model).
   40 IF (KZ(22).LE.1) THEN
          CALL SETUP
      END IF
*
   50 CONTINUE
*
      RETURN
*
      END
