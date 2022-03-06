      SUBROUTINE TCIRC2(ES0,I1,I2,ICIRC,TC)
*
*
*       Circularization time for chain binary.
*       --------------------------------------
*
*       Theory of Rosemary Mardling, Ap. J. XX, YYY, 1995.
*       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*
      INCLUDE 'common4.h'
      PARAMETER  (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &            NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIJ(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      REAL*8  WW(3),QQ(3),W(2),Q(2),AT0(2),M21,WG(2),QG(2),
     &        WSCALE(2),QSCALE(2),A(2),B(2)
      DATA  WW  /2.119,3.113,8.175/
      DATA  QQ  /0.4909,0.4219,0.2372/
      DATA  A  /6.306505,-7.297806/
      DATA  B  /32.17211,13.01598/
      SAVE  IDIAG
      DATA  IDIAG/0/
*
*
*       Define index J1 as biggest radius to be used with AT0(1).
      IF (RADIUS(I1).GE.RADIUS(I2)) THEN
          J1 = I1
          J2 = I2
      ELSE
          J1 = I2
          J2 = I1
      END IF
*
*       Define oscillation period (dimensionless time) and damping constants.
      DO 5 K = 1,2
          IF (K.EQ.1) THEN
              IK = J1
          ELSE
              IK = J2
          END IF
*       Specify polytropic index for each star (n = 3, 2 or 3/2).
          IF (ISTAR(IK).EQ.3.OR.ISTAR(IK).EQ.5) THEN
              CALL GIANT2(K,IK,WG,QG,WSCALE,QSCALE,XN,QL)
              W(K) = WG(1)
              Q(K) = QG(1)
          ELSE
              QL = 1.0D+04
              IP = 3
              IF (ISTAR(IK).GE.3) IP = 2
              IF (ISTAR(IK).EQ.4.OR.ISTAR(IK).EQ.6) IP = 3
              IF (ISTAR(IK).EQ.0) IP = 1
              W(K) = WW(IP)
              Q(K) = QQ(IP)
          END IF
          TL = TWOPI*SIZE(IK)*SQRT(SIZE(IK)/M(IK)/W(K))
          AT0(K) = 1.0/(QL*TL)
    5 CONTINUE
*
*       Form mass, radius & pericentre ratio.
      IF (SIZE(I1).GE.SIZE(I2)) THEN
          M21 = M(I2)/M(I1)
          R21 = SIZE(I2)/SIZE(I1)
	  RP1 = QPERI/SIZE(I1)
      ELSE
	  M21 = M(I1)/M(I2)
          R21 = SIZE(I1)/SIZE(I2)
	  RP1 = QPERI/SIZE(I2)
      END IF
*
*	Evaluate damping coefficient.
      RR = RP1*(1.0 + ES0)
      CONST = 2.0*(AT0(1)*(Q(1)/W(1))**2*(1.0 + M21)*M21 + 
     &             AT0(2)*(Q(2)/W(2))**2*((1.0 + M21)/M21**2)*R21**8)/
     &                                                         RR**8
*
*       Adopt WD scaling for any NS to avoid numerical problem.
      IF (ISTAR(I1).EQ.9.OR.ISTAR(I2).EQ.9) THEN
          CONST = 1.0D-04*CONST
      END IF
*
*	Form rational function approximation to Hut solution.
      FF = (( A(2)*ES0 + A(1))*ES0 + 1.0 )/
     &     (( B(2)*ES0 + B(1))*ES0 + 1.0 )
*     FF = MIN(FF,0.999D0)
*
*       Evaluate circularization time (in units of 10**6 yrs).
      TC = TSTAR*(1.0 - FF)/CONST
*
*       Activate KSTIDE indicator if TC < 2x10**9 yrs or hyperbolic orbit.
      IF (TC.LT.2000.0.OR.ES0.GT.1.0) THEN
          SEMI = QPERI/(1.0 - ES0)
          IDIAG = IDIAG + 1
          IF (IDIAG.LE.10.OR.ES0.GT.0.9) THEN
*             IF(NDIAG.LT.20000)THEN
                 WRITE (6,20)  I1, ES0, RP1, M21, TC, SEMI
   20            FORMAT (' TCIRC2:    I1 E RP/S M21 TC A ',
     &                                I4,F8.4,F8.1,F6.2,1P,2E10.2)
*                NDIAG = NDIAG + 1
*             ENDIF
          END IF
          ICIRC = 1
      END IF
*
      RETURN
*
      END
