      SUBROUTINE XTRNLV(I1,I2)
*
*
*       External potential energy.
*       --------------------------
*
      INCLUDE 'common6.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2
      REAL*8 XI(3)
      SAVE FIRST
      LOGICAL FIRST
      DATA FIRST /.TRUE./
*
*
      ET = 0.0D0
*       Skip external potential during scaling (parameters not defined).
      IF (FIRST.AND.TIME.EQ.0.0D0) THEN
          FIRST = .FALSE.
          GO TO 30
      ELSE IF (KZ(14).GT.2) THEN
*       Ignore the general case for now.
          GO TO 30
      END IF
*
*       See whether to include a linearized galactic tidal force.
      IF (KZ(14).LE.2) THEN
          DO 10 I = I1,I2
              ET = ET - 0.5D0*BODY(I)*(TIDAL(1)*X(1,I)**2 +
     &                                 TIDAL(3)*X(3,I)**2)
   10     CONTINUE
      ELSE IF (KZ(14).EQ.3) THEN
          I = I1
          RG2 = 0.0
          RI2 = 0.0
          DO 20 K = 1,3
              XI(K) = X(K,I)
              RG2 = RG2 + RG(K)**2
              RI2 = RI2 + (RG(K) + XI(K))**2
   20     CONTINUE
*
*       Include galaxy point mass term for body #I in differential form.
          IF (GMG.GT.0.0D0) THEN
              ET = ET + GMG*(1.0/SQRT(RI2) - 1.0/SQRT(RG2))
          END IF
*
*       Add optional Miyamoto disk potential.
          IF (DISK.GT.0.0D0) THEN
              R2 = (RG(1) + XI(1))**2 + (RG(2) + XI(2))**2
              BZ = SQRT(B**2 + (RG(3) + XI(3))**2)
              AZ = SQRT(R2 + (A + BZ)**2)
              R20 = RG(1)**2 + RG(2)**2
              BZ0 = SQRT(B**2 + RG(3)**2)
              AZ0 = SQRT(R20 + (A + BZ0)**2)
              ET = ET + DISK*(1.0/AZ - 1.0/AZ0)
          END IF
*
*       Check addition of differential logarithmic potential.
          IF (V02.GT.0.0D0) THEN
              ET = ET + 0.5*V02*(LOG(RI2) - LOG(RG2))
          END IF
*       Form the differential potential energy due to tides.
          ET = BODY(I)*ET
      END IF
*
*       Place sum in ETIDE and single particle contribution in HT.
   30 IF (I2.GT.I1) THEN
          ETIDE = ET
      ELSE
          HT = ET
      END IF
*
      RETURN
*
      END
