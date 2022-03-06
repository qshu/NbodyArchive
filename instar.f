      SUBROUTINE INSTAR
*
*
*       Initialization of stars.
*       ------------------------
*
      INCLUDE 'common6.h'
      REAL  LUMS(6),TSCLS(6),ZM,TM,TN,RM,TE,AGE
      DATA  PC,AU  /3.0857E+18,2.0627E+05/
*
*
*       Initialize mass loss variables & counters.
      TPHYS = 0.0D0
      ZMRG = 0.0D0
      ZMHE = 0.0D0
      ZMRS = 0.0D0
      ZMWD = 0.0D0
      ZMSN = 0.0D0
      ZMDOT = 0.0D0
      EMDOT = 0.0D0
      NMDOT = 0
      NRG = 0
      NHE = 0
      NRS = 0
      NWD = 0
      NSN = 0
      NTZ = 0
      NBS = 0
      AGE = 0.0
      TMDOT = 1.0E+10
      IF (KZ(27).EQ.0) THEN
          TSTAR = TSCALE
      END IF
*
      DO 10 I = 1,N
*
*       Obtain stellar parameters at current epoch.
          ZM = BODY(I)*ZMBAR
          CALL STAR(ZM,TM,TN,TSCLS,LUMS)
          CALL HRDIAG(ZM,AGE,TM,TN,LUMS,TSCLS,RM,TE,KW)
*
*       Convert from solar radii to scaled units (assume Sun = 0.005 AU).
          RADIUS(I) = 0.005*RM*(RBAR/AU)
*
*       Take initial evolution time 0.1 - 0.5 the main sequence time scale.
          FM = 0.5 - 0.4*(BODY(I) - BODY(N))/(BODY1 - BODY(N))
          TEV(I) = FM*TM/TSTAR
*
*       Determine the time for next stellar evolution check.
          IF (TEV(I).LT.TMDOT) THEN
              TMDOT = TEV(I)
          END IF
*
*       Initialize the stellar classification type (KW = 1 - 8).
          KSTAR(I) = KW
*
*       Save the initial mass of the 500 heaviest stars in sequential order.
          IF (I.LE.500) THEN
              BODY0(I) = BODY(I)
          END IF
   10 CONTINUE
*
      RETURN
*
      END
