      SUBROUTINE ECIRC(RP,ES0,I1,I2,TG,TC,ECC1,EDOT)
*
*
*       Eccentricity for given circularization time.
*       --------------------------------------------
*
*       Theory of Rosemary Mardling, Ap. J. XX, YYY, 1995.
*       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*
      INCLUDE 'common6.h'
      REAL*8  WW(3),QQ(3),W(2),Q(2),AT0(2),M21,WG(2),QG(2),
     &        WSCALE(2),QSCALE(2),A(2),B(2),C(6)
      DATA  WW  /2.119d0,3.113d0,8.175d0/
      DATA  QQ  /0.4909d0,0.4219d0,0.2372d0/
      DATA  A  /6.306505d0,-7.297806d0/
      DATA  B  /32.17211d0,13.01598d0/
      DATA  C  /5.101417d0,24.71539d0,-9.627739d0,1.733964d0,
     &                            -2.314374d0,-4.127795d0/
*
*
*       Specify index J1 as biggest radius to be used with AT0(1).
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
          IF (KSTAR(IK).EQ.3.OR.KSTAR(IK).EQ.5.OR.
     &        KSTAR(IK).EQ.6.OR.KSTAR(IK).EQ.9) THEN
              IPAIR = KVEC(I1)
              CALL GIANT(IPAIR,IK,WG,QG,WSCALE,QSCALE,ZN,QL)
              W(K) = WG(1)
              Q(K) = QG(1)
          ELSE
              QL = 1.0D+04
              IP = 3
              IF (KSTAR(IK).GE.3) IP = 2
              IF (KSTAR(IK).EQ.4.OR.KSTAR(IK).EQ.7) IP = 3
              IF (KSTAR(IK).EQ.8) IP = 3
              IF (KSTAR(IK).EQ.0) IP = 1
              W(K) = WW(IP)
              Q(K) = QQ(IP)
          END IF
          TL = TWOPI*RADIUS(IK)*SQRT(RADIUS(IK)/BODY(IK)/W(K))
          AT0(K) = 1.d0/(QL*TL)
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
*     Evaluate damping coefficient.
      RR = RP1*(1.d0 + ES0)
      CONST = 2.d0*(AT0(1)*(Q(1)/W(1))**2*(1.d0 + M21)*M21 + 
     &             AT0(2)*(Q(2)/W(2))**2*((1.d0 + M21)/M21**2)*R21**8)/
     &                                                         RR**8
*
*       Adopt WD scaling for any NS to avoid numerical problem.
      IF (KSTAR(I1).EQ.13.OR.KSTAR(I2).EQ.13) THEN
          CONST = 1.0D-04*CONST
      END IF
*
*     Form rational function approximation to Hut solution.
      FF = (( A(2)*ES0 + A(1))*ES0 + 1.d0 )/
     &     (( B(2)*ES0 + B(1))*ES0 + 1.d0 )
*     FF = MIN(FF,0.999D0)
*
*       Determine eccentricity corresponding to t_{circ} = t_{grow}.
      Z = TG*CONST/TSTAR + FF
      ECC1 = (-1.d0 + C(1)*Z - SQRT(C(2)*Z**2 + C(3)*Z + C(4)))
     &     /(C(5) + C(6)*Z)
*
*       Evaluate circularization time (in units of 10**6 yrs).
      TC = TSTAR*(1.d0 - FF)/CONST
*
*       Obtain (de/dt) due to tidal circularization.
      FE = 1.d0 + 3.75d0*ES0**2 + 1.875d0*ES0**4 + (5.d0/64.d0)*ES0**6
      FE = (9.d0*TWOPI/10.d0)*ES0*(1.d0 - ES0**2)**1.5d0*FE
      EDOT = -CONST*FE
*
      RETURN
*
      END
