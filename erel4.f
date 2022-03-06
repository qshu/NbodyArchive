      SUBROUTINE EREL4(IM,EB,SEMI)
*
*
*       Dominant two-body energy in four-body system.
*       ---------------------------------------------
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      COMMON/CREG/  M(4),X(12),XD(12),P(12),Q(12),TIME4,ENERGY,EPSR2,
     &              XR(9),W(9),R(6),TA(6),MIJ(6),CM(10),RMAX4,TMAX,
     &              DS,TSTEP,EPS,NSTEP4,NAME4(4),KZ15,KZ27,NREG,NFN
      COMMON/CLOSE/  RIJ(4,4),RCOLL,QPERI,SIZE(4),ECOLL3
      COMMON/CCOLL/  QK(12),PK(12),ICALL,ICOLL,NDISS4
      COMMON/SAVEP/  PI(12)
      COMMON/KSAVE/  K1,K2
*
*
*       Obtain physical momenta of current configuration defined by QK & PK.
      CALL TRANS4
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
*       Form kinetic energy terms of dominant c.m. (K1 + K2) and K3 & K4.
      P0 = 0.0D0
      P3 = 0.0D0
      P4 = 0.0D0
      J1 = 3*(K1 - 1)
      J3 = 3*(K3 - 1)
      J4 = 3*(K4 - 1)
*
      DO 5 K = 1,3
          P0 = P0 + (PI(J3+K) + PI(J4+K))**2
          P3 = P3 + PI(J3+K)**2
          P4 = P4 + PI(J4+K)**2
    5 CONTINUE
*
*       Evaluate potential energy due to non-singular terms.
      POT = 0.0D0
      DO 10 I = 1,6
          IF (I.EQ.IM) GO TO 10
          POT = POT + MIJ(I)/R(I)
   10 CONTINUE
*
*       Obtain binding energy from total energy and perturbing function.
      MB = M(K1) + M(K2)
      VP = 0.5D0*(P0/MB + P3/M(K3) + P4/M(K4)) - POT
      EB = ENERGY - VP
*
*       Set semi-major axis.
      SEMI = -M(K1)*M(K2)/(2.0D0*EB)
*
      RETURN
*
      END
