      SUBROUTINE XTRNL0
*
*
*       External force initialization.
*       ------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Check option for cluster in circular galactic orbit.
      IF (KZ(14).NE.1) GO TO 20
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
*       Define GM in cgs units and length scale in pc.
      GM = 6.67D-08*1.989D+33
      PC = 3.0857D+18
*
*       Specify the tidal term in star cluster units (solar mass & pc).
      TIDAL(1) = 4.0*A*(A - B)*(PC/GM)
*
*       Initialize the Y-component to zero.
      TIDAL(2) = 0.0
*
*       Specify the vertical force gradient.
      TIDAL(3) = -(2.0*TWOPI*RHO + 2.0*(A - B)*(A + B)*(PC/GM))
*
      FAC = 1.0E-10/(PC/GM)
      WRITE (6,5)  ZMBAR*ZMASS, FAC*TIDAL(1), FAC*TIDAL(3), PC/GM
    5 FORMAT (/,12X,'TOTAL MASS =',F7.1,'  TIDAL(1&3) =',1P,2E10.2,
     &              '  PC/GM =',E10.2)
*
*       Adopt twice the angular velocity for Coriolis terms.
      TIDAL(4) = 2.0*(A - B)*SQRT(PC/GM)
*
*       Define time scale in seconds and velocity scale in km/sec.
      TSCALE = SQRT(PC/GM)*PC
      VSTAR = 1.0E-05*SQRT(GM/PC)
*
*       Convert time scale from units of seconds to million years.
      TSCALE = TSCALE/(3.15D+07*1.0D+06)
*
*       Scale to working units of RBAR in pc & ZMBAR in solar masses.
      DO 10 K = 1,3
          TIDAL(K) = TIDAL(K)*RBAR**3/ZMBAR
   10 CONTINUE
      TIDAL(4) = TIDAL(4)*SQRT(RBAR**3/ZMBAR)
      TSCALE = TSCALE*SQRT(RBAR**3/(ZMASS*ZMBAR))
      VSTAR = VSTAR*SQRT(ZMASS*ZMBAR/RBAR)
*
*       Define tidal radius in scaled units.
      RTIDE = (ZMASS/TIDAL(1))**0.3333
*
      WRITE (6,15)  (TIDAL(K),K=1,4), TSCALE, RTIDE
   15 FORMAT (/,12X,'TIDAL PARAMETERS:  ',1P4E10.2,'  TSCALE =',E9.2,
     &                              ' (10**6 YRS)','  RTIDE =',0PF6.2,/)
*
   20 RETURN
*
      END
