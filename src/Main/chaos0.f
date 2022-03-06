      SUBROUTINE CHAOS0(QPERI,ECC,EB0,ZJ0,M1,M2,S1,S2,W,ECRIT,AR,BR,
     &                                                           IDIS)
*
*
*       Initial chaos boundary parameters.
*       ----------------------------------
*
*       Theory of Rosemary Mardling, Ap. J. XX, YYY, 1995.
*       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*
      IMPLICIT REAL*8 (A-H,M,O-Z)
      REAL*8  W(4)
      DATA  C0,C1,C2,C3,C4,C5
     &   /-1.763D-03,0.1277d0,-3.743D-03,3.525d-03,-6.179d-04,2.744d-05/
      DATA  CB /1.836d0/
      DATA  A0,A1,A2,A3 /0.079d0,-0.022d0,1.275d0,-0.670d0/
*
*
*       Form chaos boundary parameter and decide which star used for scaling.
      P = ((M2/M1)*W(2)/W(1))**0.4d0*(S1/S2)
      IF (P.GT.1.d0) THEN
          S = S1
          M = M1
          M21 = M2/M1
          W2 = W(1)  
      ELSE
          S = S2
          M = M2
          M21 = M1/M2
          W2 = W(2)  
      END IF
*
*       Specify energy & angular momentum mass factors.
      CE = 0.5d0*M1*M2
      ZMU = M1*M2/(M1 + M2)
      CJ = ZMU*SQRT(M1 + M2)
*
*       Define RP.
      RP = 0.5d0*QPERI*(1.d0 + ECC)
*
*       Write polynomial for E(c).
      Y = LOG((2.d0*RP/(CB*S))**5*(W2/(1.d0 + M21))**2)
      EC = ((((C5*Y + C4)*Y + C3)*Y + C2)*Y + C1)*Y + C0
*
*       Compare EC to ECDIS (disruption value for EC).
      IOUT = IDIS
      IDIS = 0
      ECDIS = ((A3*ECC + A2)*ECC + A1)*ECC + A0
      IF (EC.LT.ECDIS) THEN
          IDIS = 1
      END IF
*
*	Check whether spiral is indicated (unless collision or ECC > 1).
      IF ((ECC.LT.EC).AND.(IDIS.EQ.0).AND.(ECC.LT.1.d0)) THEN
         IDIS = -1
         RETURN
      END IF
*
*       Obtain R(c) & Ecrit.
      RC = 2.d0*RP/(1.d0 + EC)
      ECRIT = CE*(EC - 1.d0)/RC
*
*       Calculate RP(b) & E(b).
      ALF2 = (2.d0/SQRT(W2))*SQRT(S**3/M)
      CFAC = (2.d0*ECRIT - EB0)/CE
      YFAC = (ZJ0 + 2.d0*ALF2*(ECRIT - EB0))/CJ
      TMP = 1.d0 + CFAC*YFAC**2
      TMP = MAX(TMP,0.0D0)
      EB = SQRT(TMP)
      RB = (EB - 1.d0)*CE/(2.d0*ECRIT - EB0)
*
*       Form the A & B parameters for determining minimum eccentricity.
      AR = (RB - RC)/(EB - EC)
      BR = (EB*RC - EC*RB)/(EB - EC)
*       Note: AR & BR in units of length and ECRIT in mass/radius.
*
*       ---------------------------------------------------------------
*       Include occasional diagnostic output (skip CHRECT/CHAOS).
      IF (IOUT.GE.0) THEN
          WRITE (6,1)  M1, M2, M21, QPERI/S, EC, ECDIS, EB
    1     FORMAT (' INIT:    M1 M2 M21 QP/S EC ECD EB ',
     &                       1P,2E10.2,0P,2F5.1,3F6.2)
      END IF
*       ---------------------------------------------------------------
      RETURN
*
      END
