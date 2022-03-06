      SUBROUTINE INTIDE
*
*
*       Input & scaling for tidal capture.
*       ----------------------------------
*
      INCLUDE 'common6.h'
      REAL  LUMS(6),TSCLS(6),ZM,TM,TN,RM,TE,AGE
*
*
*       Read parameters for tidal capture simulation.
      READ (5,*)  RSTAR, RTSTAR, RSYNC, EPOCH
*
*       Define GM in cgs units and length scale in pc.
      GM = 6.67E-08*1.989E+33
      PC = 3.0857E+18
*
*       Form time scale in seconds and velocity scale in km/sec.
      TSTAR = SQRT(PC/GM)*PC
      VSTAR = 1.0E-05*SQRT(GM/PC)
*
*       Convert time scale from units of seconds to million years.
      TSTAR = TSTAR/(3.15E+07*1.0E+06)
*
*       Ensure ZMBAR & RBAR > 0 (=0: assume <M>/Sun = 1, RBAR = 1 pc).
      IF (ZMBAR.LE.0.0D0) ZMBAR = FLOAT(N)/ZMASS
      IF (RBAR.LE.0.0D0) RBAR = 1.0
*
*       Scale to working units of RBAR in pc & ZMBAR in solar masses.
      TSTAR = TSTAR*SQRT(RBAR**3/(ZMASS*ZMBAR))
      VSTAR = VSTAR*SQRT(ZMASS*ZMBAR/RBAR)
*
*       Convert stellar radii from A.U. to internal length units.
      AU = PC/1.495979E+13
      RSTAR = RSTAR/(AU*RBAR)
      RTSTAR = RTSTAR/(AU*RBAR)
      RSYNC = RSYNC/(AU*RBAR)
*
*       Express evolution epoch & time scale in internal time units.
      AGE = EPOCH
      EPOCH = EPOCH/TSTAR
      TSCALE = TSTAR
*
*       Initialize event counter.
      DO 5 K = 1,10
          NTYPE(K) = 0
    5 CONTINUE
*
      WRITE (6,10)  RSTAR, RTSTAR, RSYNC, EPOCH
   10 FORMAT (/,12X,'TIDAL CAPTURE:   RSTAR =',1PE8.1,'  RTSTAR =',E8.1,
     &                                '  RSYNC =',E8.1,'  EPOCH =',E8.1)
*
      WRITE (6,20)  RBAR, ZMBAR, VSTAR, TSTAR
   20 FORMAT (/,12X,'PHYSICAL SCALING:   R* =',F6.3,'  M* =',F7.1,
     &                                    '  V* =',F6.3,'  T* =',F6.3,/)
*
*       Physical scaling: X, M, V, T from RBAR*X, ZMBAR*M, VSTAR*V, TSTAR*T.
*
      READ (5,*)  IMS, IEV, RMS, REV, IMF
      WRITE (6,25)  IMS, IEV, RMS, REV, IMF
   25 FORMAT (/,12X,'STELLAR RADII:   IMS =',I5,'  IEV =',I4,
     &           '  RMS/RSTAR =',F5.2,'  REV/RSTAR =',F6.3,'  IMF =',I2)
*
*       Assign individual radii for main-sequence and evolved stars.
      DO 30 I = 1,N
          IF (IMF.GT.0) THEN
              ZM = BODY(I)*ZMBAR
*
*       Obtain radius at current epoch from stellar evolution models.
              CALL STAR(ZM,TM,TN,TSCLS,LUMS)
              CALL HRDIAG(ZM,AGE,TM,TN,LUMS,TSCLS,RM,TE,KW)
*
*       Convert from solar radii to scaled units (assume Sun = 0.005 AU).
              RADIUS(I) = 0.005*RM*(RBAR/AU)
*       Scale by non-zero fudge factor (experimental).
              IF (RMS.GT.0.0) RADIUS(I) = RMS*RADIUS(I)
              GO TO 30
          END IF
*
*       Retain the old scheme just in case.
          IF (I.LE.IMS) THEN
              RADIUS(I) = RMS*RSTAR
          ELSE
              RADIUS(I) = REV*RSTAR
          END IF
   30 CONTINUE
*
      WRITE (6,40)  BODY(1), BODY(N), RADIUS(1), RADIUS(N)
   40 FORMAT (/,12X,'SCALED RADII:   M(1) =',F8.4,'  M(N) =',F8.4,
     &                                '  R(1) =',1PE8.1,'  R(N) =',E8.1)
*
      RETURN
*
      END
