      PROGRAM NEWMASS
*
*
*       Generation of initial coordinates & velocities.
*       -----------------------------------------------
*       for special mass distributions R.Sp. Sta. Cruz Aug 2000
*
      INCLUDE 'common6.h'
      REAL*8  A(8)
      REAL*4  RAN2
      REAL*8 C(3),CMR(3),CMRDOT(3),RHOC1(NMAX)
      INTEGER IDUM,NRAND
*
*
      XFAC = 1.D0
*       Choose between uniform density and Plummer model.
      TWOPI = 8.0D0*ATAN(1.D0)
      XMINT = 1.D0
      IPASS = 0
*
      PRINT*,' Enter N,NRAND,Q,mu,IMOD:'
      READ*,N,NRAND,QQ,XMU,IMOD
*
      IDUM = -NRAND
*
*       xmu is m2/m1, qq is M2/Mtot
*
      IF(XMU.LT.1.D0)THEN
      PRINT*,' mu should always be greater than 1'
      STOP
      END IF
*
      XM1TOT = 1.D0 - QQ
      XM2TOT = QQ
      XM1 = (1.D0 - (1.D0-XMU)*XM1TOT)/(XMU*N)
      XM2 = XMU*XM1
*       Theoretical values of particle numbers
      XN1TOT = XM1TOT/XM1
      XN2TOT = XM2TOT/XM2
      XMAV = 1.D0/N
*
      PRINT*,' Comp 1 M,N,m=',XM1TOT,XN1TOT,XM1
      PRINT*,' Comp 2 M,N,m=',XM2TOT,XN2TOT,XM2
      PRINT*,' Total  M,N,m=',XM1TOT+XM2TOT,XN1TOT+XN2TOT,
     *      (XM1TOT+XM2TOT)/(XN1TOT+XN2TOT),XMAV
      PRINT*,' XMAV/XM1,XMAV/XM2=',XMAV/XM1,XMAV/XM2
*
*       Now prepare the real particles by random numbers
*
      N1 = INT(XN1TOT)
      N2 = INT(XN2TOT)
      ZMASS = 0.D0
      DO 1 I = 1,N
      IF(I.LE.N2)THEN
      BODY(I) = XM2
      IND(I) = 2
      ZMASS = ZMASS + BODY(I)
      ELSE
      BODY(I) = XM1
      IND(I) = 1
      ZMASS = ZMASS + BODY(I)
      END IF
 1    CONTINUE
*
      PRINT*,' N,ZMASS,N1,N2=',N,ZMASS,N-N2,N2
*       Final Normalization of mass
      DO 2 I = 1,N
 2    BODY(I) = BODY(I)/ZMASS
*
      ZMASS = 1.D0
*       Save random number sequence to start at same point second time.
*
       WRITE(91,*)IDUM,IY,IDUM2,
     *           (IV(KK),KK=1,NTAB)
 1000 CONTINUE
*       Retry to get correct total mass second time
      REWIND 91
      READ(91,*,ERR=9995,END=9995)IDUM,IY,IDUM2,
     *           (IV(KK),KK=1,NTAB)
      IPASS = IPASS + 1
*       Initialize centre of mass terms.
      DO 25 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   25 CONTINUE
*
      RMAX = 20.D0
      RHALF = (0.5D0**(-0.6666667D0) - 1.D0)**(-0.5D0)
*       Compute maximum value for comparison with random number
      FXMAX = 7.D0
*
      PRINT*,' Before first RAN'
      CALL FLUSH(6)
      PRINT*,' First RAN = ',RAN2()
*       Generate initial conditions using rejection technique.
      DO 140 I = 1,N
*
 130  RR = RMAX*RAN2()
      RR2 = RR*RR
*       Determine average mass variation function
      XRR = RR/RHALF
      XR2 = XFAC*XRR*XRR
      XRI2 = XFAC/XRR/XRR
      ALPHA = XM1/XMAV
      BETA = XM2/XMAV
*       Ansatz for f(r) where 1/<m> = 1/<m>_0 f(r)
*       Divide by XMINT to get correct total mass in second pass.
      IF(IMOD.EQ.0)THEN
      XFR = 1.D0
      ELSE IF(IMOD.EQ.1)THEN
      XFR = (XR2+1.D0)/(ALPHA*XR2+BETA)
      ELSE IF(IMOD.EQ.2)THEN
      XFR = (XRI2+1.D0)/(ALPHA*XRI2+BETA)
      ELSE
      PRINT*,' IMOD=',IMOD,' not defined '
      END IF
*
      XFUNC(I) = XFR
      XCORR = 1.D0/(XMU-1.D0)
*       Which particle sort coming next?
      IF(I.LE.N2)THEN
*       Density correction factor must be 1/qi if XFR=1
      RHOC(I) = -XMU*XCORR*(ALPHA*XFR-1.D0)/XM2TOT
      RHOC1(I) = XCORR*(BETA*XFR-1.D0)/XM1TOT
      ELSE
      RHOC(I) = XCORR*(BETA*XFR-1.D0)/XM1TOT
      RHOC1(I) = -XMU*XCORR*(ALPHA*XFR-1.D0)/XM2TOT
      END IF
*
      RHO(I) = 3.D0/(1.D0+RR2)**2.5D0
      RHOR = RHOC(I)*RHO(I)*RR2
      IF(RHOR.GT.FXMAX)PRINT*,' Warning I,R,RHOR,FXMAX=',I,RR,RHOR,FXMAX
*
      XRAN = RAN2()
      IF(XRAN*FXMAX.LT.RHOR)THEN
*       Take the radius if not rejected by Monte Carlo
      RI = RR
*
          A(2) = RAN2()
          A(3) = RAN2()
          X(3,I) = (1.0 - 2.0*A(2))*RI
          X(1,I) = SQRT(RI**2 - X(3,I)**2)*COS(TWOPI*A(3))
          X(2,I) = SQRT(RI**2 - X(3,I)**2)*SIN(TWOPI*A(3))
  132     A(4) = RAN2()
          A(5) = RAN2()
          A(6) = A(4)**2*(1.0 - A(4)**2)**3.5
          IF (0.1*A(5).GT.A(6)) GO TO 132
*
          A(8) = A(4)*SQRT(2.0)/(1.0 + RI**2)**0.25
          A(6) = RAN2()
          A(7) = RAN2()
          XDOT(3,I) = (1.0 - 2.0*A(6))*A(8)
          XDOT(1,I) = SQRT(A(8)**2 - XDOT(3,I)**2)*COS(TWOPI*A(7))
          XDOT(2,I) = SQRT(A(8)**2 - XDOT(3,I)**2)*SIN(TWOPI*A(7))
*       
*       Accumulate c.m. terms.
          DO 135 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
  135     CONTINUE
*
          ELSE
*
          IREJ = IREJ + 1
          GOTO 130
*
          END IF
*
  140 CONTINUE
*
*       Scale coordinates & velocities to analytical expectation values.
      SX = 1.5*TWOPI/16.0
      SV = SQRT(ZMASS/SX)
      DO 150 I = 1,N
          DO 145 K = 1,3
*             X(K,I) = X(K,I) - CMR(K)/ZMASS
*             XDOT(K,I) = XDOT(K,I) - CMRDOT(K)/ZMASS
              X(K,I) = SX*X(K,I)
              XDOT(K,I) = SV*XDOT(K,I)
  145     CONTINUE
      RHO(I) = RHO(I)/SX**3/(TWOPI*2.D0)
  150 CONTINUE
*
*
      PRINT*,' Going to call ENERGY IPASS = ',IPASS
      CALL FLUSH(6)
*
      CALL ENERGY(ZKIN,POT)
*
      PRINT*,' ENERGY: ZKIN,POT,ETOT,Q=',ZKIN,POT,ZKIN-POT,ZKIN/POT
*
*       Scale total energy to standard units (E = -0.25 for QVIR < 1).
*     QVIR = 0.5D0
*     E0 = -0.25
*     E0 = ZKIN - POT
*     ETOT = (QVIR - 1.0)*POT
*     SX = E0/ETOT
*     SX = QVIR*POT/ZKIN
*
*     DO 70 I = 1,N
*         DO 68 K = 1,3
*             X(K,I) = X(K,I)/SX
*             XDOT(K,I) = XDOT(K,I)*SQRT(SX)
*  68     CONTINUE
*     RHO(I) = RHO(I)*SX**3
*  70 CONTINUE

*     PRINT*,' Going to call ENERGY after scaling'
*     CALL FLUSH(6)
*
*     CALL ENERGY(ZKIN,POT)
*
*     PRINT*,' ENERGY: ZKIN,POT,ETOT,Q=',ZKIN,POT,ZKIN-POT,ZKIN/POT
*
      CALL LAGR
*
      IF(DABS(XMINT-1.D0).GT.1.D-2)THEN
      IF(XMINT.GT.1.D0)THEN
      XFAC = (XFAC + XFAC/(1.D0+1.D2*DABS(XMINT-1.D0)))/2.D0
      ELSE
      XFAC = (XFAC + XFAC*(1.D0+1.D2*DABS(XMINT-1.D0)))/2.D0
      END IF
      PRINT*,' IPASS,XMINT,XMINT-1,XFAC=',IPASS,XMINT,XMINT-1.D0,XFAC
      GOTO 1000
      END IF
*
      DO 996 I=1,N
      WRITE(10,1110)BODY(I),(X(K,I),K=1,3),(XDOT(K,I),K=1,3)
 996  CONTINUE
      DO 997 I=N2+1,N
      WRITE(11,1111)BODY(I),(X(K,I),K=1,3),(XDOT(K,I),K=1,3),PHIDBL(I),
     *RHO(I),RHOC(I),RHOC1(I),XMR1(I),XMR2(I),XMRPL(I),RHOPL(I),XFUNC(I)
 997  CONTINUE
      DO 998 I=1,N2
      WRITE(12,1111)BODY(I),(X(K,I),K=1,3),(XDOT(K,I),K=1,3),PHIDBL(I),
     *RHO(I),RHOC(I),RHOC1(I),XMR1(I),XMR2(I),XMRPL(I),RHOPL(I),XFUNC(I)
 998  CONTINUE
*
 1110 FORMAT(1X,1P,7(1X,D15.7))
 1111 FORMAT(1X,1P,16(1X,D15.7))
*
      STOP
*
 9995 PRINT*,' Error in reading random number data'
      STOP
      END
