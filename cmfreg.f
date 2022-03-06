      SUBROUTINE CMFREG(I,RS2,RCRIT2,VRFAC,NNB,XI,XID,FIRR,FREG,FD,FDR)
*
*
*       Regular & irregular c.m. force & first derivative.
*       --------------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  XI(3),XID(3),FIRR(3),FREG(3),DV(3),FD(3),FDR(3)
*       
*       mass dependency of neighbours:
      REAL*8 rs2fac
*
*
*       Adopt accurate force for perturbed c.m. particle.
      IF (I.GT.N) THEN
          IPAIR = I - N
          I2 = 2*IPAIR
          I1 = I2 - 1
          RPERT2 = CMSEP2*R(IPAIR)**2
          BODYIN = 1.0/BODY(I)
          IF (GAMMA(IPAIR).GE.GMIN) GO TO 10
      END IF
*
*       Copy all KS pairs to JLIST and find MAX(R) for joint treatment.
      NNB1 = NPAIRS
      RMAX1 = 0.0
      DO 1 LJ = 1,NNB1
          JLIST(LJ) = N + LJ
          RMAX1 = MAX(RMAX1,R(LJ))
    1 CONTINUE
*
*       Set minimum square distance for c.m. approximation.
      RCM2 = MAX(RCRIT2,CMSEP2*RMAX1**2)
*       Use accurate force algorithm for J > N.
      GO TO 30
*
*       Use fast force loop for particles satisfying c.m. approximation.
   10 RCM2 = MAX(RCRIT2,RPERT2)
      NNB1 = 0
      DO 20 J = IFIRST,NTOT
          A1 = X(1,J) - XI(1)
          A2 = X(2,J) - XI(2)
          A3 = X(3,J) - XI(3)
          DV(1) = XDOT(1,J) - XID(1)
          DV(2) = XDOT(2,J) - XID(2)
          DV(3) = XDOT(3,J) - XID(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          CALL rs2modf(rs2fac,BODY(I),BODY(J))
*
*       Form a list of c.m. particles to be resolved.
          IF (RIJ2.LT.(rs2fac*RCM2)) THEN
              NNB1 = NNB1 + 1
              JLIST(NNB1) = J
              GO TO 20
          END IF
*
          DR2I = 1.0/RIJ2
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
          FREG(1) = FREG(1) + A1*DR3I
          FREG(2) = FREG(2) + A2*DR3I
          FREG(3) = FREG(3) + A3*DR3I
          FDR(1) = FDR(1) + (DV(1) - A1*DRDV)*DR3I
          FDR(2) = FDR(2) + (DV(2) - A2*DRDV)*DR3I
          FDR(3) = FDR(3) + (DV(3) - A3*DRDV)*DR3I
   20 CONTINUE
*
*       Begin dual purpose force loop (all RIJ2 < RCM2 or J > N).
   30 DO 60 LJ = 1,NNB1
          JDUM = JLIST(LJ)
          A1 = X(1,JDUM) - XI(1)
          A2 = X(2,JDUM) - XI(2)
          A3 = X(3,JDUM) - XI(3)
          DV(1) = XDOT(1,JDUM) - XID(1)
          DV(2) = XDOT(2,JDUM) - XID(2)
          DV(3) = XDOT(3,JDUM) - XID(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          CALL rs2modf(rs2fac,BODY(I),BODY(JDUM))
*       First see whether the distance exceeds c.m. approximation limit.
          IF (RIJ2.GT.(rs2fac*RCM2)) GO TO 56
*
          J = JDUM
*       Test whether particle #J satisfies neighbour criteria.
          IF (RIJ2.GT.(rs2fac*RCRIT2)) GO TO 54
*
          IF (JDUM.EQ.I) GO TO 60
          IF (RIJ2.GT.(rs2fac*RS2).AND.STEP(J).GT.SMIN) THEN
*       Omission of small step particle may give large correction terms.
              A7 = A1*DV(1) + A2*DV(2) + A3*DV(3)
*       Accept member if maximum penetration factor exceeds 8 per cent.
              IF (A7.GT.VRFAC) GO TO 54
          END IF
*
*       Obtain force due to current neighbours.
          NNB = NNB + 1
          ILIST(NNB) = J
          ICM = 1
*       Irregular force indicator in case I > N or J > N are resolved.
          IF (J.GT.N) THEN
*       See whether c.m. approximation applies (skip perturbed case).
              J1 = 2*(J - N) - 1
              IF (RIJ2.LT.rs2fac*CMSEP2*R(J-N)**2.AND.LIST(1,J1).GT.0)
     &             GO TO 50
          END IF
*
          IF (I.GT.N) THEN
*       See whether c.m. force needs summation over each component.
              IF (RIJ2.LT.(rs2fac*RPERT2).AND.GAMMA(IPAIR).GE.GMIN)
     &            GO TO 40
          END IF
*
          DR2I = 1.0/RIJ2
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
          FIRR(1) = FIRR(1) + A1*DR3I
          FIRR(2) = FIRR(2) + A2*DR3I
          FIRR(3) = FIRR(3) + A3*DR3I
          FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
          FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
          FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
          GO TO 60
*
*       Obtain relevant force on c.m (ICM = 0 denotes resolved pair).
   40     K = J
   42     L = I1
*       Individual components I1 & I2 are resolved in routine INTGRT.
   45     A1 = X(1,K) - X(1,L)
          A2 = X(2,K) - X(2,L)
          A3 = X(3,K) - X(3,L)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          CALL rs2modf(rs2fac,BODY(K),BODY(L))
          DV(1) = XDOT(1,K) - XDOT(1,L)
          DV(2) = XDOT(2,K) - XDOT(2,L)
          DV(3) = XDOT(3,K) - XDOT(3,L)
*
          DR2I = 1.0/RIJ2
          DR3I = BODY(K)*BODY(L)*DR2I*SQRT(DR2I)*BODYIN
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
          IF (ICM.NE.0) THEN
              FIRR(1) = FIRR(1) + A1*DR3I
              FIRR(2) = FIRR(2) + A2*DR3I
              FIRR(3) = FIRR(3) + A3*DR3I
              FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
              FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
              FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
          ELSE
              FREG(1) = FREG(1) + A1*DR3I
              FREG(2) = FREG(2) + A2*DR3I
              FREG(3) = FREG(3) + A3*DR3I
              FDR(1) = FDR(1) + (DV(1) - A1*DRDV)*DR3I
              FDR(2) = FDR(2) + (DV(2) - A2*DRDV)*DR3I
              FDR(3) = FDR(3) + (DV(3) - A3*DRDV)*DR3I
          END IF
*
          L = L + 1
          IF (L.EQ.I2) GO TO 45
          K = K + 1
          IF (K.EQ.J + J) GO TO 42
*
          GO TO 60
*
   50     J = J - N
          J2 = J + J
*       Sum over the components (unperturbed case is OK).
          K = J2 - 1
*       Treat all interactions between components of two perturbed pairs.
          IF (I.GT.N) THEN
              IF (RIJ2.LT.(rs2fac*RPERT2).AND.GAMMA(I-N).GT.GMIN)
     &            GO TO 42
          END IF
*
   52     A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          DV(1) = XDOT(1,K) - XID(1)
          DV(2) = XDOT(2,K) - XID(2)
          DV(3) = XDOT(3,K) - XID(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          CALL rs2modf(rs2fac,BODY(I),BODY(K))
*
          DR2I = 1.0/RIJ2
          DR3I = BODY(K)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
          IF (ICM.NE.0) THEN
              FIRR(1) = FIRR(1) + A1*DR3I
              FIRR(2) = FIRR(2) + A2*DR3I
              FIRR(3) = FIRR(3) + A3*DR3I
              FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
              FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
              FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
          ELSE
              FREG(1) = FREG(1) + A1*DR3I
              FREG(2) = FREG(2) + A2*DR3I
              FREG(3) = FREG(3) + A3*DR3I
              FDR(1) = FDR(1) + (DV(1) - A1*DRDV)*DR3I
              FDR(2) = FDR(2) + (DV(2) - A2*DRDV)*DR3I
              FDR(3) = FDR(3) + (DV(3) - A3*DRDV)*DR3I
          END IF
*
          K = K + 1
          IF (K.EQ.J2) GO TO 52
          GO TO 60
*
*       Resolve components if I > N or J > N and no c.m. approximation.
   54     IF (J.GT.N) THEN
              IF (RIJ2.LT.rsfac*CMSEP2*R(J-N)**2) THEN
*       Set zero indicator to denote current pair resolved for regular force.
                  ICM = 0
                  GO TO 50
              END IF
          END IF
*
          IF (I.GT.N) THEN
              IF (RIJ2.LT.(rs2fac*RPERT2)) THEN
                  ICM = 0
                  GO TO 40
              END IF
          END IF
*
*       Obtain the regular force due to single body or c.m. particle.
   56     DR2I = 1.0/RIJ2
          DR3I = BODY(JDUM)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
          FREG(1) = FREG(1) + A1*DR3I
          FREG(2) = FREG(2) + A2*DR3I
          FREG(3) = FREG(3) + A3*DR3I
          FDR(1) = FDR(1) + (DV(1) - A1*DRDV)*DR3I
          FDR(2) = FDR(2) + (DV(2) - A2*DRDV)*DR3I
          FDR(3) = FDR(3) + (DV(3) - A3*DRDV)*DR3I
   60 CONTINUE
*
*       Check force correction due to regularized chain (G < GMIN in REGINT).
      IF (I.GT.N.AND.NCH.GT.0) THEN
          IF (GAMMA(I-N).GT.GMIN) THEN
              CALL KCPERT(I,I1,FIRR,FD)
          END IF
      END IF 
*
      RETURN
*
      END
