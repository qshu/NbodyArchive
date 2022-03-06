      SUBROUTINE HIVEL(IH)
*
*
*       High-velocity particle search.
*       ------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Copy membership and add 1 for analogy with HARP.
      NHI = LISTV(1) + 1
*
*       Specify square velocity limit in terms of current state.
      VMAX2 = 16.0*ECLOSE
*
*       Check for removal of distant high-velocity particle.
      IF (IH.LT.0) THEN
    1     LH = 0
          DO 2 L = 2,NHI
              I = LISTV(L)
              VI2 = X0DOT(1,I)**2 + X0DOT(2,I)**2 + X0DOT(3,I)**2
              RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                                       (X(3,I) - RDENS(3))**2
*       Save index of fast particle outside 3*<R> or VI < VMAX/4.
              IF (RI2.GT.9.0*RSCALE**2.OR.VI2.LT.VMAX2/16.0) THEN
                  LH = L
                  LI = I
                  RL2 = RI2
                  VL2 = VI2
              END IF
    2     CONTINUE
*       Reduce membership and remove any distant member from the list.
          IF (LH.GT.0) THEN
              if(rank.eq.0)
     &        WRITE (29,3)  LI, NAME(LI), SQRT(RL2), SQRT(VL2)
    3         FORMAT (' HIVEL REMOVE    I NAM R V ',2I6,2F6.1)
              NHI = NHI - 1
              LISTV(1) = LISTV(1) - 1
              DO 4 L = LH,NHI
                  LISTV(L) = LISTV(L+1)
    4         CONTINUE
              GO TO 1
          END IF
      END IF
*
*       Set index to known neutron star or terminated KS/chain.
      IF (IH.GT.0) THEN
          I1 = IH
          I2 = IH
      ELSE
*       Include two first single particles or up to five chain members.
          I1 = IFIRST
          I2 = IFIRST + 1
          IF (IPHASE.EQ.8) THEN
              I2 = IFIRST + 4
*       Search all particles after escaper removal (defined by IPHASE = -1).
          ELSE IF (IH.EQ.0.AND.IPHASE.EQ.-1) THEN
              I2 = NTOT
              NHI = 1
              LISTV(1) = 0
          ELSE IF (IPHASE.EQ.1) THEN
*       See whether the first few locations may have been exchanged.
              DO 5 L = 2,NHI
                  IF (LISTV(L).LT.IFIRST + 3) THEN
                      I2 = NTOT
                  END IF
    5         CONTINUE
          END IF
      END IF
*
*       Add any new high-velocity particles (skip F**2 > N & STEP < DTMIN).
      DO 10 I = I1,I2
          FI2 = F(1,I)**2 + F(2,I)**2 + F(3,I)**2
*       Adopt geometric mean of close encounter and core mass force.
          IF (FI2.GT.FLOAT(N).OR.STEP(I).LT.DTMIN.OR.
     &        STEPR(I).LT.20.0*DTMIN) GO TO 10
          VI2 = X0DOT(1,I)**2 + X0DOT(2,I)**2 + X0DOT(3,I)**2
          RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                                   (X(3,I) - RDENS(3))**2
          IF (VI2.GT.VMAX2.AND.RI2.LT.9.0*RSCALE**2) THEN
              DO 8 L = 2,NHI
                  IF (I.EQ.LISTV(L)) GO TO 10
    8         CONTINUE
*       Check maximum membership and possible ghost particle.
              IF (NHI.GE.MLV.OR.STEP(I).GT.1.0) GO TO 10
              NHI = NHI + 1
              LISTV(1) = LISTV(1) + 1
              NFAST = NFAST + 1
              LISTV(NHI) = I
              if(rank.eq.0)
     &        WRITE (29,9)  TIME+TOFF, NHI, I, NAME(I), KSTAR(I),
     &                      SQRT(VI2), SQRT(RI2), STEP(I)
    9         FORMAT (' HIVEL ADD    T NHI I NM K* VI R DT ',
     &                               F10.4,I4,2I6,I4,2F6.2,1P,E10.2)
          END IF
   10 CONTINUE
*
*       Increase counter if loop is over all particles.
*     IF (I2.EQ.NTOT) THEN
*         NHIVEL = NHIVEL + 1
*     END IF
*
      RETURN
*
      END

