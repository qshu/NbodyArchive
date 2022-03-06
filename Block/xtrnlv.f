      SUBROUTINE XTRNLV(I1,I2)
*
*
*       External potential and virial energy.
*       -------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2,GMB,AR,GAM,ZDUM(7)
      REAL*8 XI(3),XIDOT(3),FIRR(3),FREG(3),FD(3),FDR(3)
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
*       Treat all cases uniformly (but note Coliolis term in FIRR).
      ELSE IF (KZ(14).LE.4.AND.I2.NE.I1) THEN
          VIR = 0.0
          DO 5 I = I1,I2
              DO 2 K = 1,3
*       Initialize scalars for force contributions.
                  XI(K) = X(K,I)
                  XIDOT(K) = XDOT(K,I)
                  FIRR(K) = 0.0
                  FREG(K) = 0.0
                  FD(K) = 0.0
                  FDR(K) = 0.0
    2         CONTINUE
*       Evaluate the virial energy using definition sum {m*r*F}.
              CALL XTRNLF(XI,XIDOT,FIRR,FREG,FD,FDR,1)
*       Note that FIRR(K) only contributes for #14 = 1 or 2.
              DO 4 K = 1,3
                  VIR = VIR + BODY(I)*XI(K)*(FIRR(K) + FREG(K))
    4         CONTINUE
    5     CONTINUE
*       Quit for pure 3D tidal field (case of full N summation).
          IF (KZ(14).EQ.3) GO TO 40
*       Note ETIDE is accumulated after each regular step (REGINT/GPUCOR).
      END IF
*
*       See whether to include a linearized galactic tidal force.
      IF (KZ(14).LE.2) THEN
          DO 10 I = I1,I2
              CALL XTRNLV2(I,ET1)
              ET = ET + ET1
   10     CONTINUE
      ELSE IF (KZ(14).EQ.3) THEN
          CALL XTRNLV2(I1,ET)
      END IF
*
*       Place sum in ETIDE and single particle contribution in HT.
   30 IF (I2.GT.I1) THEN
          ETIDE = ET
      ELSE
          HT = ET
      END IF
*
   40 RETURN
*
      END
