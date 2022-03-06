      SUBROUTINE NBPOT2(NB,NP,POTS,KLIST,KPERT)
*
*
*       Potential energy of subsystem.
*       ------------------------------
*
      INCLUDE 'common6.h'
      INTEGER KLIST(LMAX),KPERT(LMAX)
*
*
*     WRITE (6,2)  (KLIST(K),K=1,NB)
*   2 FORMAT (' KLIST   ',10I5)
*     WRITE (6,4)  (KPERT(K),K=1,NP)
*   4 FORMAT (' KPERT   '10I5,4(/,5X,10I5))
*       Obtain potential energy of subsystem KLIST(NB) & KPERT(NP).
      POTS = 0.0D0
      DO 10 L = 1,NB
          I = KLIST(L)
          DO 5 K = 1,NP
              J = KPERT(K)
              IF (J.EQ.I) GO TO 5
              IF (J.GT.N) THEN
*       Skip outer binary c.m. or components during merger correction.
                  JP = J - N
                  IF (I.LT.IFIRST) THEN
                      IF (JP.EQ.KVEC(I)) GO TO 5
                  END IF
*       Resolve c.m. of other perturbed KS pairs.
                  IF (LIST(1,2*JP-1).GT.0) THEN
                      J = 2*JP - 1
                  END IF
              END IF
    1         RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                      (X(3,I) - X(3,J))**2
              POTS = POTS + BODY(I)*BODY(J)/SQRT(RIJ2)
*
*       Check for second component of perturbed KS pair.
              IF (KPERT(K).GT.N) THEN
                  IF (J.EQ.2*JP - 1) THEN
                      J = J + 1
                      GO TO 1
                  END IF
              END IF
    5     CONTINUE
   10 CONTINUE
*
*       Include any external potentials.
      IF (KZ(14).GT.0) THEN
          DO 20 L = 1,NB
              I = KLIST(L)
              CALL XTRNLV2(I,ET)
              POTS = POTS - ET
*       Note positive sign convention for potential energy.
   20     CONTINUE
      END IF
*
      RETURN
*
      END
