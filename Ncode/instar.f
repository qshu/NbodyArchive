      SUBROUTINE INSTAR
*
*
*       Initialization of stars.
*       ------------------------
*
      INCLUDE 'common6.h'
      REAL*8  LUMS(10),TSCLS(20),GB(10)
      REAL*8  M0,M1,LUM,MC
      EXTERNAL  RZAMSF
*
*
*       Initialize mass loss variables & counters.
      TPHYS = 0.0D0
      TPLOT = 0.0
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
      NBH = 0
      NTZ = 0
      NBS = 0
      NKICK = 0
      NBKICK = 0
      AGE = 0.D0
      TMDOT = 1.0E+10
      IF (KZ(27).EQ.0) THEN
          TSTAR = TSCALE
      END IF
*
*       Set the Hydrogen & Helium abundances.
      yhel = 0.24d0 + 2.d0*zmet
      xhyd = 1.d0 - yhel - zmet
*
*       Set evolution parameters which depend solely on metallicity.
      CALL zcnsts(zmet,zpars)
*
      DO 10 I = 1,N
*
*       Obtain stellar parameters at current epoch.
          M1 = BODY(I)*ZMBAR
          M0 = M1
          IF (M1.GT.0.1) THEN
              KW = 1
              AGE = 0.0
              CALL STAR(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
              CALL HRDIAG(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                    RM,LUM,KW,MC,RCC)
*       Set initial look-up time 0.05 of main sequence time (but < 10**7 yrs).
              TEV(I) = MIN(0.05*TM,10.0D0)/TSTAR
          ELSE
              KW = 0
              RM = RZAMSF(M1)
              TEV(I) = 1.0D+10/TSTAR
          END IF
*
*       Convert from solar radii to scaled units (assume Sun = 0.005 AU).
          RADIUS(I) = 0.005*RM/(AU*RBAR)
*
*       Determine the time for next stellar evolution check.
          IF (TEV(I).LT.TMDOT) THEN
              TMDOT = TEV(I)
          END IF
*
*       Initialize the stellar classification type (KW = 1 - 8).
          KSTAR(I) = KW
*
*       Save the initial mass of all the stars in sequential order.
          BODY(I) = M1/ZMBAR
          BODY0(I) = M0/ZMBAR
   10 CONTINUE
*
*       Define first quantized step < 1000 yrs (minimum interval for MDOT).
      DT = 1.0E-03/TSCALE
      CALL STEPK(DT,DTN)
      STEPX = DTN
*
      RETURN
*
      END
