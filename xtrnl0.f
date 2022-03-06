      SUBROUTINE XTRNL0
*
*
*       External force initialization.
*       ------------------------------
*
      INCLUDE 'common4.h'
      REAL*8  R2,RHOS
      COMMON/WORK1/  R2(NMAX),RHOS(NMAX)
*
*
*       Form time scale in seconds and velocity scale in km/sec.
      TSCALE = SQRT(PC/GM)*PC
      VSTAR = 1.0E-05*SQRT(GM/PC)
*
*       Convert time scale from units of seconds to million years.
      TSCALE = TSCALE/(3.15D+07*1.0D+06)
*
*       Copy total mass from basic power-law IMF (ZMBAR used for scaling).
      IF (KZ(20).EQ.0) THEN
          ZMBAR = ZMTOT
      END IF
*
*       Specify Oort's constants (units of km/sec/kpc).
      A = 14.4
      B = -12.0
*       Adopt local density from Gilmore & Kuijken (solar mass/pc**3).
      RHO = 0.11
*       Convert rotation constants to units of cm/sec/pc.
      A = 100.0*A
      B = 100.0*B
*
*       Specify the tidal term in star cluster units (solar mass & pc).
      TIDAL(1) = 4.0*A*(A - B)*(PC/GM)
*
*       Initialize the Y-component to zero.
      TIDAL(2) = 0.0
*
      IF(KZ(14).NE.-1)THEN
*       Specify the vertical force gradient.
         TIDAL(3) = -(2.0*TWOPI*RHO + 2.0*(A - B)*(A + B)*(PC/GM))
*       Adopt twice the angular velocity for Coriolis terms.
         TIDAL(4) = 2.0*(A - B)*SQRT(PC/GM)
      ENDIF
*
      FAC = 1.0E-10/(PC/GM)
      WRITE (6,1)  ZMBAR*ZMASS, FAC*TIDAL(1), FAC*TIDAL(3), PC/GM
    1 FORMAT (/,12X,'TOTAL MASS =',F8.1,'  TIDAL(1&3) =',1P,2E10.2,
     &              '  PC/GM =',E10.2)
*
*       Scale to working units of RBAR in pc & ZMBAR in solar masses.
      RHOBAR = ZMASS*ZMBAR/RBAR**3
      DO 5 K = 1,3
          TIDAL(K) = TIDAL(K)/RHOBAR
    5 CONTINUE
      TIDAL(4) = TIDAL(4)/SQRT(RHOBAR)
*
*       Scale to working units of RBAR in pc & ZMBAR in solar masses.
      TSCALE = TSCALE/SQRT(RHOBAR)
      TSTAR = TSCALE
      VSTAR = VSTAR*SQRT(ZMASS*ZMBAR/RBAR)
*
*       Define tidal radius in scaled units.
      RTIDE = (ZMASS/TIDAL(1))**0.3333
*
      WRITE (6,40)  (TIDAL(K),K=1,4), TSCALE, RTIDE
   40 FORMAT (/,12X,'TIDAL PARAMETERS:  ',1P,4E10.2,'  TSCALE =',E9.2,
     &                              ' (10**6 YRS)','  RTIDE =',0PF6.2,/)
*
      IF (KZ(12).GT.0) THEN
*       Specify disk shock interval (Myr, truncated) and set next shock time.
          DTSHOCK = TWOPI*UT/TIDAL(4)
          DTSHOCK = 1.0 + INT(DTSHOCK)
          TSHOCK = DTSHOCK
          NSHOCK = 0
          WRITE (6,50)  DTSHOCK
   50     FORMAT (/,12X,'DISK SHOCK:    DTSHOCK =',F7.2)
      END IF
*
      RETURN
*
      END
