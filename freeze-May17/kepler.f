      SUBROUTINE KEPLER(I,DTI)
*
*
*       Step reduction of hierarchy.
*       ----------------------------
*
*
      INCLUDE 'common6.h'
      COMMON/ISAVE/  LI0,LI,NS,NSLIST(LMAX)
*
      IPAIR = I - N
      I1 = 2*IPAIR - 1
*
*       Only examine hierarchical configurations (NP < 5 OR STEP < DTMIN).
      IF (LIST(1,I1).GT.4.AND.DTI.GT.DTMIN) GO TO 30
*
*       Set new close encounter index, list membership & Kepler period.
      ICLOSE = I
      NP1 = LIST(1,I1) + 1
      SEMI = -0.5D0*BODY(I)/H(I-N)
      TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
*
*       Consider the most dominant perturbers having comparable step to c.m.
      DO 20 L = 2,NP1
          J = LIST(L,I1)
          IF (STEP(J).GT.4.0*STEP(I)) GO TO 20
          RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                  (X(3,I) - X(3,J))**2
          RIJ = SQRT(RIJ2)
          DT2 = 0.14*RIJ*SQRT(ETAI*RIJ/(BODY(I) + BODY(J)))
          DT = 0.25D0*SQRT(ETAI)*RIJ*TK/SEMI
*       Compare predicted c.m. step with conservative Kepler expressions.
          DT = MIN(DT,DT2)
          IF (DTI.LT.DT) GO TO 20
          DTI = DT
*
*       Check whether to reduce step of dominant perturber.
          IF (STEP(J).LT.DT) GO TO 20
          IF (KZ(36).LE.0) THEN
              A3 = T0(J) + STEP(J)
*       Adopt maximum reduction factor 0.5 ensuring T0 + STEP > TIME.
              STEP(J) = STEP(J) - 0.5*(A3 - TIME)
          ELSE IF (T0(J).EQ.TIME) THEN
              STEP(J) = 0.5D0*STEP(J)
              TIMENW(J) = T0(J) + STEP(J)
              GO TO 20
          END IF
          JCLOSE = J
*
*       Consider body #J for NSLIST if new & old step straddle TLIST.
          A3 = T0(J) + STEP(J)
          IF (A3.GT.TLIST) GO TO 20
*
*       Note that body #J may be a member if treated recently.
          DO 5 LL = 1,NS
              IF (NSLIST(LL).EQ.J) GO TO 20
    5     CONTINUE
*
*       Include #J in NSLIST or reduce sequential position in NLIST.
          IF (A3.GT.TLIST) THEN
              NS = NS + 1
              NSLIST(NS) = J
          ELSE
              NL = NLIST(1) + 1
*       Skip NLIST modification if body #I is last member.
              IF (LI.GE.NL) GO TO 20
*
*       Find current NLIST index (LJ = 0 is OK for next loop).
              LJ = 0
              DO 10 LL = LI+1,NL
                  IF (NLIST(LL).EQ.J) THEN
                      LJ = LL
                      GO TO 12
                  END IF
   10         CONTINUE
*
*       Determine new NLIST index > LI in case #J needs to be moved.
   12         LN = 0
              DO 14 LL = LI+1,LJ
                  JJ = NLIST(LL)
                  IF (A3.LT.T0(JJ) + STEP(JJ)) THEN
                      LN = LL
                      GO TO 16
                  END IF
   14         CONTINUE
*
*       Move NLIST members down by one and replace #J in location LN.
   16         IF (LN.GT.0) THEN
                  LL = LJ
   18             NLIST(LL) = NLIST(LL-1)
                  LL = LL - 1
                  IF (LL.GT.LN) GO TO 18
                  NLIST(LN) = J
              END IF
          END IF
   20 CONTINUE
*
   30 RETURN
*
      END
