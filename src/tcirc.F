      SUBROUTINE TCIRC(RP,ES0,I1,I2,ICIRC,TC)
!VAMPINST include
*
*
*       Pre-main sequence circularization.
*       ----------------------------------
*
*       Theory of Piet Hut, A & A 99, 126 (1981).
*       Developed by Rosemary Mardling (31/1/97).
*       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*
      INCLUDE 'common6.h'
      REAL*8  WW(3),QQ(3),W(2),Q(2),AT0(2),M21,
     &        A(2),B(2),C(6)
      DATA  WW  /2.119,3.113,8.175/
      DATA  QQ  /0.4909,0.4219,0.2372/
      DATA  A  /6.306505,-7.297806/
      DATA  B  /32.17211,13.01598/
      DATA  C  /5.101417,24.71539,-9.627739,1.733964,
     &                            -2.314374,-4.127795/
      DATA  ECCM /0.002/
*
*
*       Specify index J1 as biggest radius to be used with AT0(1).
!VAMPINST subroutine_start
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
          QL = 1.0D+04
          IP = 3
          IF (KSTAR(IK).GE.3) IP = 2
          IF (KSTAR(IK).EQ.4) IP = 3
          IF (KSTAR(IK).EQ.0) IP = 1
          W(K) = WW(IP)
          Q(K) = QQ(IP)
          TL = TWOPI*RADIUS(IK)*SQRT(RADIUS(IK)/BODY(IK)/W(K))
          AT0(K) = 1.0/(QL*TL)
    5 CONTINUE
*
*       Form mass, radius & pericentre ratio.
      IF (RADIUS(I1).GE.RADIUS(I2)) THEN
          M21 = BODY(I2)/BODY(I1)
          R21 = RADIUS(I2)/RADIUS(I1)
	  RP1 = RP/RADIUS(I1)
      ELSE
	  M21 = BODY(I1)/BODY(I2)
          R21 = RADIUS(I1)/RADIUS(I2)
	  RP1 = RP/RADIUS(I2)
      END IF
*
*	Evaluate damping coefficient.
      RR = RP1*(1.0 + ES0)
      CONST = 2.0*(AT0(1)*(Q(1)/W(1))**2*(1.0 + M21)*M21 + 
     &             AT0(2)*(Q(2)/W(2))**2*((1.0 + M21)/M21**2)*R21**8)/
     &                                                         RR**8
*
*	Form rational function approximation to Hut solution.
      FF = (( A(2)*ES0 + A(1))*ES0 + 1.0 )/
     &     (( B(2)*ES0 + B(1))*ES0 + 1.0 )
*
      TIME0 = 0.0
      IF (TC.LT.0.0) TIME0 = TC
*
*       Obtain the new eccentricity.
      Z = (TIME - TIME0)*CONST + FF
      ECC = (-1.0 + C(1)*Z - SQRT(C(2)*Z**2 + C(3)*Z + C(4)))
     &		                             /(C(5) + C(6)*Z)
*
      ECC = MAX(ECC,ECCM)
      RP = RP*(1.0 + ES0)/(1.0 + ECC)
      ES0 = ECC
*
!VAMPINST return
          RETURN
      RETURN
*
!VAMPINST subroutine_end
      END
