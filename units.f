      SUBROUTINE UNITS
*
*
*       Initialization of units & scaling factors.
*       ------------------------------------------
*
      INCLUDE 'common6.h'
*
      
*       Define GM & PC in cgs units and AU in pc.
      GM = 6.67D-08*1.989D+33
      PC = 3.0857D+18
      AU = 2.0627E+05
*
*       Form scaling factors for binary periods A*SQRT(A/M) to yrs and days.
      IF(ZMBAR.GT.0.D0)THEN
      YRS = (RBAR*AU)**1.5/SQRT(ZMBAR)
      DAYS = 365.24*YRS
      ELSE
      YRS = 0.D0
      DAYS = 0.D0
      END IF
*
*       Specify conversion factors for lengths to solar radii & AU.
      SU = PC/(AU*6.96D+10)*RBAR*AU
      RAU = RBAR*AU
*
*       Copy solar mass scaling to new variable (M = BODY*<M>).
      SMU = ZMBAR
*
*       Form time scale in seconds and velocity scale in km/sec.
      TSTAR = SQRT(PC/GM)*PC
      VSTAR = 1.0D-05*SQRT(GM/PC)
*
*       Convert time scale from units of seconds to million years.
      TSTAR = TSTAR/(3.15D+07*1.0D+06)
*
*       Ensure ZMBAR & RBAR > 0 (=0: assume <M>/Sun = 1, RBAR = 1 pc).
      IF (ZMBAR.LE.0.0D0) ZMBAR = FLOAT(N)/ZMASS
      IF (RBAR.LE.0.0D0) RBAR = 1.0
*
*       Scale to working units of RBAR in pc & ZMBAR in solar masses.
      TSTAR = TSTAR*SQRT(RBAR**3/(ZMASS*ZMBAR))
      VSTAR = VSTAR*SQRT(ZMASS*ZMBAR/RBAR)
      CLIGHT = 2.99792458D5/VSTAR  
*
*       Copy TSTAR to secondary time-scale factor.
      TSCALE = TSTAR
*
*       Physical scaling: X, M, V, T from RBAR*X, ZMBAR*M, VSTAR*V, TSTAR*T.
      if(rank.eq.0)then
      WRITE (6,10) RBAR, ZMBAR, VSTAR, TSTAR, BODYM*ZMBAR, SU, RAU,YRS,
     &             CLIGHT
   10 FORMAT (/,12X,'PHYSICAL SCALING:    R* =',F4.1,'  M* =',E15.8,
     &              '  V* =',E15.8,'  T* =',1PE9.2,'  <M> =',E15.8,
     &              '  SU =',1P,E8.1,' AU=',1PE9.2,' YRS=',E9.2,/,
     &              ' CLIGHT =',E9.2)
      end if
*
      RETURN
*
      END
