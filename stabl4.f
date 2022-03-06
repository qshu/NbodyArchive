      SUBROUTINE STABL4(ITERM)
*
*
*       Stability test of four-body system.
*       -----------------------------------
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      LOGICAL  SWITCH,GTYPE,GTYPE0,AUX
      COMMON/CREG/  M(4),X(12),XD(12),P(12),Q(12),TIME4,ENERGY,EPSR2,
     &              XR(9),W(9),R(6),TA(6),MIJ(6),CM(10),RMAX4,TMAX,
     &              DS,TSTEP,EPS,NSTEP4,NAME4(4),KZ15,KZ27,NREG,NFN
      COMMON/TPR/   SWITCH,GTYPE,GTYPE0
      COMMON/KSAVE/  K1,K2
      REAL*8  RC(3),VC(3),RC0(3),VC0(3)
*
*
*       Transform to physical variables (retain value of SWITCH).
      AUX = SWITCH
      SWITCH = .FALSE.
      CALL ENDREG
      SWITCH = AUX
*
*       Specify indices of two least dominant bodies (denoted K3 & K4).
      K3 = 0
      DO 1 L = 1,4
          IF (L.EQ.K1.OR.L.EQ.K2) GO TO 1
          IF (K3.EQ.0) THEN
              K3 = L
          ELSE
              K4 = L
          END IF
    1 CONTINUE
*
*       Initialize scalars for orbital elements.
      RB0 = 0.0D0
      RB = 0.0D0
      RB1 = 0.0D0
      RDOT = 0.0D0
      RDOT1 = 0.0D0
      VREL2 = 0.0D0
      VREL20 = 0.0D0
      VREL21 = 0.0D0
*
*       Define binary masses of smallest & widest pair (K1 & K2 and K3 & K4).
      MB0 = M(K1) + M(K2)
      MB = M(K3) + M(K4)
*
*       Form separations & velocities of MB0, MB and their relative orbit.
      DO 10 K = 1,3
          J1 = 3*(K1-1) + K
          J2 = 3*(K2-1) + K
          J3 = 3*(K3-1) + K
          J4 = 3*(K4-1) + K
          RC0(K) = (M(K1)*X(J1) + M(K2)*X(J2))/MB0
          RC(K) = (M(K3)*X(J3) + M(K4)*X(J4))/MB
          VC0(K) = (M(K1)*XD(J1) + M(K2)*XD(J2))/MB0
          VC(K) = (M(K3)*XD(J3) + M(K4)*XD(J4))/MB
          RB0 = RB0 + (X(J1) - X(J2))**2
          RB = RB + (X(J3) - X(J4))**2
          RB1 = RB1 + (RC(K) - RC0(K))**2
          RDOT = RDOT + (X(J3) - X(J4))*(XD(J3) - XD(J4))
          RDOT1 = RDOT1 + (RC(K) - RC0(K))*(VC(K) - VC0(K))
          VREL2 = VREL2 + (XD(J3) - XD(J4))**2
          VREL20 = VREL20 + (XD(J1) - XD(J2))**2
          VREL21 = VREL21 + (VC(K) - VC0(K))**2
   10 CONTINUE
*
*       Determine semi-major axis of inner binary.
      RB0 = SQRT(RB0)
      SEMI0 = 2.0/RB0 - VREL20/MB0
      SEMI0 = 1.0/SEMI0
      E0 = 1.0 - RB0/SEMI0
*
*       Form semi-major axis & eccentricity of wide pair.
      RB = SQRT(RB)
      SEMI = 2.0D0/RB - VREL2/MB
      SEMI = 1.0/SEMI
      E = SQRT((1.0D0 - RB/SEMI)**2 + RDOT**2/(SEMI*MB))
*
*       Evaluate orbital elements of relative c.m. motion.
      RB1 = SQRT(RB1)
      SEMI1 = 2.0D0/RB1 - VREL21/CM(7)
      SEMI1 = 1.0/SEMI1
      E1 = SQRT((1.0D0 - RB1/SEMI1)**2 + RDOT1**2/(SEMI1*CM(7)))
*
*       Obtain standard stability ratio (outer pericentre / inner apocentre).
      RATIO = SEMI1*(1.0D0 - E1)/(SEMI0*(1.0D0 + E0))
*
*       Form hierarchical stability ratio (Kiseleva & Eggleton 1995).
      Q0 = MB/MB0
      Q1 = MAX(M(K2)/M(K1),M(K1)/M(K2))
      Q3 = Q0**0.33333
      Q13 = Q1**0.33333
      AR = 1.0 + 3.7/Q3 - 2.2/(1.0 + Q3) + 1.4/Q13*(Q3 - 1.0)/(Q3 + 1.0)
*
      PCRIT = AR*SEMI0*(1.0D0 + E0)
      PMIN = SEMI1*(1.0D0 - E1)
*
*       Set negative termination index if system is stable or SEMI < 0.
      IF (PCRIT.LT.PMIN) THEN
          ITERM = -1
          WRITE (6,20)  SEMI, SEMI1, E, E1, RB1, RATIO, PCRIT, PMIN,
     &                  SEMI0, RB0
   20     FORMAT ('  STABQ:  A A1 E E1 R1 RATIO PCR PM A0 RB0 ',
     &                            1P,2E10.2,0P,2F6.2,F9.5,F6.2,1P,4E9.1)
      ELSE
          ITERM = 0
      END IF
*
      RETURN
*
      END
