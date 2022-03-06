      SUBROUTINE UNPERT(IPAIR,KPERT)
*
*
*       Unperturbed two-body motion.
*       ----------------------------
*
      INCLUDE 'common6.h'
*
*
*       Set component & c.m. index and semi-major axis.
      I1 = 2*IPAIR - 1
      I = N + IPAIR
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
*
*       Add elapsed unperturbed Kepler periods and update the time T0(I1).
      TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
      K = NINT((TIME - T0(I1))/TK)
      T0(I1) = TIME
*
*       Reset unperturbed counter if > 2*10**9 and update the frequency.
      IF (NCOUNT(13).GT.2000000000) THEN
          NCOUNT(13) = 0
          NCOUNT(10) = NCOUNT(10) + 1
      END IF
      NCOUNT(13) = NCOUNT(13) + K
      IF (NKSPER.GT.2000000000) THEN
          NKSPER = 0
          NPRECT = NPRECT + 1
      END IF
      NKSPER = NKSPER + K
*
*       Evaluate interval for unperturbed motion (GAMMA < GMIN).
      CALL TPERT(IPAIR,GMIN,DT)
*
*       Restore KS indicator and set restart index if interval < period.
      IF (DT.LT.TK) THEN
          KSLOW(IPAIR) = 1
          KPERT = 1
          GO TO 10
      END IF
*
*       Specify a conservative number of unperturbed orbits.
      K = 1 + INT(0.5D0*DT/TK)
      STEP(I1) = FLOAT(K)*TK
*
*       Check merger condition before continuing unperturbed motion.
C Correct:
c 10   IF (KZ(15).NE.0.AND.STEP(I).LT.DTMIN) THEN
C Old:
      IF (KZ(15).NE.0.AND.STEP(I).LT.DTMIN) THEN
*       Skip possible second call during termination (KSTERM calls KSINT).
          IF (IPHASE.EQ.0) THEN
              CALL IMPACT(I)
          END IF
      END IF
*
*       Set KPERT = -1 to prevent perturbed restart in KSINT.
      KPERT = -1
*
C Correct:
c      RETURN
C Old:
 10   RETURN
*
      END



