      SUBROUTINE ENERGY
*
*
*       Total energy.
*       -------------
*
      INCLUDE 'common6.h'
*
      REAL*8 A1B,A2B,A3B,A4B
*       Sum the total energy of regularized pairs.
      EBIN = 0.0D0
      DO 10 IPAIR = 1,NPAIRS
*       Skip pairs with zero mass of c.m. particle (merged binary ghost).
          IF (BODY(N+IPAIR).GT.0.0D0) THEN
*       Predict coordinates, velocities & binding energy.
              CALL RESOLV(IPAIR,1)
              EBIN = EBIN + BODY(2*IPAIR-1)*BODY(2*IPAIR)*HT/
     &                                                     BODY(N+IPAIR)
          END IF
   10 CONTINUE
*
*       Calculate the potential energy.
      ZKIN = 0.D00
      POT = 0.0
*
      DO 20 I = 1,NTOT
      JMIN = I + 1
      IF (I.LE.2*NPAIRS) THEN
*       Binding energy of regularized pairs is included explicitly above.
          IPAIR = KVEC(I)
          JMIN = 2*IPAIR + 1
      END IF
*
      IPAIR = 0
      IF (I.GT.N)  THEN
*       Binding energy at center of mass position without binary members
          IPAIR = I - N
      END IF
*
      POTJ = 0.D00
      POTI = 0.D00
*       POTI contains potential at particles position to be stored later (R.Sp.)
*
      DO 30 J = 1,N
      IF (J.EQ.I .OR. J.EQ.2*IPAIR-1 .OR. J.EQ.2*IPAIR .OR.
     *    BODY(J).EQ.0.0D0 .OR. BODY(I).EQ.0.0D0)  GO TO 30
          A1 = X(1,I) - X(1,J)
          A2 = X(2,I) - X(2,J)
          A3 = X(3,I) - X(3,J)
      A4 = BODY(J)/DSQRT (A1*A1 + A2*A2 + A3*A3) 
      POTI = POTI - A4 
*  also J.LT.N?
      IF(J.GE.JMIN) POTJ = POTJ + A4 
   30 CONTINUE
*     
*       Store potential in shared vector first (R.Sp.)
      PHIDBL(I) = POTI 
      POT = POT + BODY(I)*POTJ
   20 CONTINUE
*
*       Sum the kinetic energy (include c.m. bodies but not components).
      EBLCKHL=0.0
      DO 40 I = IFIRST,NTOT
      ZKIN = ZKIN + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 +
     &                                          XDOT(3,I)**2)
*       Accumulation of Black Hole Energy
*      IF(TIME.GT.0.0) THEN
      EBLCKHL=EBLCKHL-CMBLHOLE*BODY(I)/DSQRT(X(1,I)**2+
     &        X(2,I)**2+X(3,I)**2)
*      END IF
   40 CONTINUE
      ZKIN = 0.5D0*ZKIN
*     PRINT*,'BLACK HOLE BINDING ENERGY=',EBLCKHL,ZKIN,POT
*
*       Obtain the tidal potential if external field is present.
      ETIDE = 0.0D0
      IF (KZ(14).GT.0) THEN
          CALL XTRNLV(1,N)
      END IF
*
*       Check differential potential energy due to chain subsystem.
      IF (NCH.GT.0) THEN
          CALL CHPOT(DP)
          POT = POT + DP
      END IF
*
*       Total energy = ZKIN - POT + ETIDE + EBIN + ESUB + EMERGE + ECOLL.
*
      RETURN
*
      END





