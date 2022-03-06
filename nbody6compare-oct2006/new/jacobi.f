      SUBROUTINE JACOBI(NESC)
*
*
*       Jacobi escape criterion.
*       ------------------------
*
      INCLUDE 'common6.h'
      REAL*8  PHI(50),PSI(50)
*
*
*       Specify escape energy (tidal field or isolated system).
      IF (KZ(14).GT.0.AND.KZ(14).LT.2) THEN
          ECRIT = -1.5*(TIDAL(1)*ZMASS**2)**0.333
      ELSE
          ECRIT = 0.0
      END IF
*
*       Count all escapers.
      NESC = 0
      IF (N.GT.10) GO TO 51
      DO 50 I = IFIRST,NTOT
          VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
          POTI = 0.0
          DO 40 J = IFIRST,NTOT
              IF (J.EQ.I) GO TO 40
              RIJ2 = 0.0
              DO 30 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
   30         CONTINUE
              POTI = POTI + BODY(J)/SQRT(RIJ2)
   40     CONTINUE
          EI = 0.5*VI2 - POTI
          IF (EI.GT.ECRIT) NESC = NESC + 1
   50 CONTINUE
   51 CONTINUE
*
      DO 60 I = IFIRST,NTOT
          RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                                   (X(3,I) - RDENS(3))**2
          VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
          IF (RI2.GT.9.0*RSCALE**2) THEN
              VC2 = ZMASS/SQRT(RI2) + ECRIT
              IF (VI2.LT.VC2) GO TO 60
          ELSE IF (RI2.GT.RSCALE**2) THEN
              VC2 = ZMASS/SQRT(RI2) + ECRIT
              IF (VI2.LT.VC2) GO TO 60
          ELSE
              VC2 = ZMASS/SQRT(RI2) + ECRIT
              IF (VI2.LT.VC2) GO TO 60
          END IF
          POTI = 0.0
          DO 55 J = IFIRST,NTOT
              IF (J.EQ.I) GO TO 55
              RIJ2 = 0.0
              DO 52 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
   52         CONTINUE
              POTI = POTI + BODY(J)/SQRT(RIJ2)
   55     CONTINUE
          EI = 0.5*VI2 - POTI
          IF (EI.GT.ECRIT) NESC = NESC + 1
   60 CONTINUE
*
      IF (N.GT.10) RETURN
      X1 = 0.05
      EPS2 = (0.01*RSCALE)**2
      DO 80 L = 1,50
      PHI(L) = 0.0
      PSI(L) = 0.0
      DO 70 I = IFIRST,NTOT
          RIJ2 = (X(1,I) - RDENS(1) - X1)**2 + (X(2,I) - RDENS(2))**2 +
     &                                         (X(3,I) - RDENS(3))**2
          PHI(L) = PHI(L) + BODY(I)/SQRT(RIJ2 + EPS2)
   70 CONTINUE
      PSI(L) = PSI(L) + ZMASS/SQRT(X1**2 + RC**2)
      WRITE (17,75) X1, PHI(L), PSI(L)
   75 FORMAT (' ',1P,3E12.4)
      CALL FLUSH(17)
      X1 = SQRT(2.0)*X1 
      IF (X1.GT.2.0*RTIDE) RETURN
   80 CONTINUE
*
      RETURN
*
      END
