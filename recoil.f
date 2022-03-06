      SUBROUTINE RECOIL(IEND)
*
*
*       Binary analysis and final output.
*       ---------------------------------
*
      INCLUDE 'COMMON1.CH'
      INCLUDE 'COMMON2.CH'
      REAL*8  M,MIJ,MB,MB1,R2(NMX,NMX)
      INTEGER  IJ(NMX)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
      COMMON/CALLS/  TPR,ICALL,NFN,NREG
      SAVE
*
*
*       Sort particle separations (I1 & I2 form closest pair).
      CALL R2SORT(IJ,R2)
      I1 = IJ(1)
      I2 = IJ(2)
      I3 = IJ(3)
      I4 = IJ(4)
*       Set redundant index for N = 3.
      IF (N.EQ.3) I4 = I1
      IF (N.EQ.2) THEN
          I3 = I2
          I4 = I1
      END IF
*
*       Form output diagnostics.
      VREL2 = 0.0D0
      VREL21 = 0.0D0
      RDOT = 0.0D0
      DO 10 K = 1,3
          J1 = 3*(I1-1) + K
          J2 = 3*(I2-1) + K
          J3 = 3*(I3-1) + K
          J4 = 3*(I4-1) + K
          VREL2 = VREL2 + (V(J1) - V(J2))**2
          VREL21 = VREL21 + (V(J3) - V(J4))**2
          RDOT = RDOT + (X(J1) - X(J2))*(V(J1) - V(J2))
   10 CONTINUE
*
*       Evaluate orbital elements.
      RB = SQRT(R2(I1,I2))
      RB1 = SQRT(R2(I3,I4))
      R13 = SQRT(R2(I1,I3))
      R24 = SQRT(R2(I2,I4))
      MB = M(I1) + M(I2)
      SEMI = 2.0D0/RB - VREL2/MB
      SEMI = 1.0/SEMI
      ECC = SQRT((1.0D0 - RB/SEMI)**2 + RDOT**2/(SEMI*MB))
      EB = -M(I1)*M(I2)/(2.0D0*SEMI)
      ET = -M(I1)*M(I2)/(2.0D0*SEMI*CM(8))
*
      IF (N.GT.3) THEN
          MB1 = M(I3) + M(I4)
          SEMI1 = 2.0D0/RB1 - VREL21/MB1
          SEMI1 = 1.0/SEMI1
          EB1 = -M(I3)*M(I4)/(2.0D0*SEMI1)
          IF (EB1.GT.0.0) EB1 = 0.0
      ELSE
          R24 = 1.0E+10
          IF (TIMEC.LE.0.0D0) EB1 = 0.0
      END IF
*
*       Ensure that perturbed boundary exceeds system size.
      RMAX = MAX(SQRT(R2(I1,I3)),SQRT(R2(I2,I4)))
      RMAXC = MAX(RMAXC,1.2*RMAX)
*
*       Save initial energies (note: I1 & I2 may not be most energetic).
      IF (IEND.EQ.0) THEN
          EB0 = EB
          EB10 = EB1
          E0 = ENERGY
          DE = 0.0
          GO TO 50
      ELSE IF (IEND.EQ.1) THEN
          DE = DE + ABS(ENERGY - E0)
          E0 = ENERGY
*       Initialize second binary on subsequent absorption.
          IF (EB10.EQ.0.0D0) EB10 = EB1
          GO TO 50
      END IF
*
*       Determine nearest and most distant binary perturber (N = 4).
      IF (R13.LT.R24) THEN
          IM = I3
          RM = R13
          IMAX = I4
          RMAX = R24
      ELSE
          IM = I4
          RM = R24
          IMAX = I3
          RMAX = R13
      END IF
*
*       Estimate perturbations on smallest binary.
      IF (N.GT.3) THEN
          GB = 2.0*MB1/MB*(RB/RM)**3
          IF (EB1.LT.0.0) GB = GB*M(IM)/MB1
          G4 = 2.0*M(IMAX)/(MB + M(IM))*(RM/RMAX)**3
      ELSE
          GB = M(I3)/MB*(RB/R13)**3
          G4 = 0.0
          EB1 = 0.0
      END IF
*
*       Form net binary energy change (copied to CHCOLL in routine CHTERM).
      CM(9) = EB + EB1 - (EB0 + EB10)
*
*       Define relative energy production and total energy change.
      DB = CM(9)/(EB0 + EB10)
      DE = DE + ABS(ENERGY - E0)
*
*       Print final binary if relative energy increase exceeds 0.1.
      IF (DB.GT.0.1) THEN
*       Scale binding energy of relative motion by initial value.
          E1 = (ENERGY - EB - EB1)/EB0
          EB = EB/ENERGY
          EB1 = EB1/ENERGY
          WRITE (6,20)  MB, SEMI, ECC, EB, GB, G4, EB1, E1, ET
   20     FORMAT (/,' CHAIN BINARY','  MB =',F7.4,'  A =',1PE8.1,
     &              '  E =',0PF5.2,'  EB =',F5.2,'  GB =',1PE8.1,
     &              '  G4 =',E8.1,'  EB1 =',0PF5.2,'  E1 =',F5.2,
     &              '  ET =',F6.3)
      END IF
*
*       Check optional print diagnostics of chain regularization.
      IF (KZ30.GT.1) THEN
          TCR = MASS**2.5/ABS(2.0D0*ENERGY)**1.5
          TC = TIMEC/TCR
          EC = ENERGY/CM(8)
          WRITE (6,30)  I1, I2, I3, I4, RB, R13, R24, RGRAV, TC, NSTEP1,
     &                  NREG, DB, EC
   30     FORMAT (/,' END CHAIN  ',4I3,'  RB =',1PE8.1,'  R13 =',E8.1,
     &              '  R24 =',E8.1,'  RG =',E8.1,'  TC =',0PF9.1,'  #',
     &                 I5,I4,'  DB =',F5.2,'  EC =',F6.3)
          WRITE (6,40)  NSTEP1, NFN, NPERT, ENERGY, DE,
     &                  (1.0/RINV(K),K=1,N-1)
   40     FORMAT (' END CHAIN   # NFN NP E DE R ',
     &                          I6,I8,I4,F10.5,1P,6E9.1)
      END IF
*
   50 RETURN
*
      END
