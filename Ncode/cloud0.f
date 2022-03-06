      SUBROUTINE CLOUD0
*
*
*       Cloud initialization.
*       ---------------------
*
      INCLUDE 'common6.h'
      COMMON/CLOUDS/  XCL(3,MCL),XDOTCL(3,MCL),BODYCL(MCL),RCL2(MCL),
     &                CLM(MCL),CLMDOT(MCL),CLDOT,VCL,SIGMA,RB2,PCL2,
     &                TCL,STEPCL,NCL,NEWCL
*
*
*       Initialize cloud variables.
      NCL = 0
      TCL = 0.0D0
      STEPCL = 1.0E+06
      NEWCL = 0
      PCL2 = 0.0
*
*       Read the cloud parameters.
      READ (5,*)  NCL, RB2, VCL, SIGMA, (CLM(J),J=1,NCL),
     &            (RCL2(J),J=1,NCL)
      WRITE (6,100)  NCL, RB2, VCL, SIGMA
  100 FORMAT (/,12X,'CLOUDS =',I3,'  BOUNDARY RADIUS (PC) =',F6.1,
     &    '  MEAN CLOUD VELOCITY (KM/SEC) =',F5.1,'  DISPERSION =',F5.1)
      RBAR1 = RBAR
      IF (RBAR.EQ.0.0) RBAR1 = 1.0
*       Set cloud parameters in scaled units.
      RB2 = RB2/RBAR1
      A1 = 0.047*SQRT(ZMASS*ZMBAR/RBAR1)
*       Rms velocity of cluster members in km/sec.
      A2 = A1/SQRT(0.5D0*ZMASS)
*       Velocity unit.
      VCL = VCL/A2
*       Cloud velocity in scaled units.
      SIGMA = SIGMA/A2
*       Specify conservative cloud integration step using crossing time.
      STEPCL = 0.002*TCR*RB2/VCL
*
*       Adopt a quantized value.
      CALL STEPK(STEPCL,DTN)
      STEPCL = DTN
*
*       Scale radii & masses to model units.
      DO 101 J = 1,NCL
          RCL2(J) = RCL2(J)/RBAR1
          CLM(J) = CLM(J)/ZMBAR
  101 CONTINUE
*
      WRITE (6,102)  RB2, VCL, SIGMA, STEPCL, (CLM(J),J=1,NCL),
     &               (RCL2(J),J=1,NCL)
  102 FORMAT (/,12X,'SCALED CLOUD PARAMETERS',3F7.1,F9.6,5F7.1,5F6.1)
      CLDOT = 0.1*RB2/VCL
*       Time scale for 'sun-rise' is 0.05 of the cloud crossing time.
      CLDOT = 1.0/CLDOT
*
*       Define the square of cloud half-mass radii & growth times.
      DO 104 J = 1,NCL
          RCL2(J) = RCL2(J)**2
          CLMDOT(J) = CLM(J)*CLDOT
  104 CONTINUE
*
*       Set square boundary radius & impact parameter.
      RB2 = RB2**2
      PCL2 = RB2
*       Define density centre for routine CLOUD.
      DO 105 K = 1,3
          RDENS(K) = 0.0
  105 CONTINUE
*
*       Initialize new clouds on the boundary.
      DO 110 ICL = 1,NCL
          CALL CLOUD(ICL)
  110 CONTINUE
*
  120 RETURN
*
      END
