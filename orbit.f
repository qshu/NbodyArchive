      SUBROUTINE ORBIT(I,J,SEMI,ECC,GI)
*
*
*       Two-body elements.
*       -----------------
*
      INCLUDE 'common6.h'
*
*       Predict, for parallel (RSp.)
          CALL JPRED(I)
*
*       Find the dominant neighbour if body #J not specified.
      GI = 0.0
      IF (J.LE.0) THEN
          NNB = LIST(1,I)
          FMAX = 1.0D-10
          JM = LIST(2,I)
          DO 2 L = 1,NNB
              JJ = LIST(L+1,I)
*        Predict coordinate, for parallel (RSp.)
              CALL JPRED(JJ)
              JLIST(L) = JJ
              RIJ2 = 0.0
              DO 1 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,JJ))**2
    1         CONTINUE
              FIJ = (BODY(I) + BODY(JJ))/RIJ2
*       Exclude any c.m. bodies from dominant motion.
              IF (FIJ.GT.FMAX.AND.JJ.LE.N) THEN
                  FMAX = FIJ
                  JM = JJ
              END IF
    2     CONTINUE
*       Avoid rare case of halo orbit with zero neighbour number.
          IF (JM.EQ.I) JM = I - 1
*
*       Obtain the relative perturbation for decision-making.
          CALL FPERT(I,JM,NNB,PERT)
          GI = PERT/FMAX
          J = JM
      END IF
*
*       Determine the semi-major axis and eccentricity.
      RIJ2 = 0.0
      VIJ2 = 0.0
      RDOT = 0.0
      DO 5 K = 1,3
          RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
          VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,J))**2
          RDOT = RDOT + (X(K,I) - X(K,J))*(XDOT(K,I) - XDOT(K,J))
    5 CONTINUE
      RIJ = SQRT(RIJ2)
      A1 = 2.0/RIJ - VIJ2/(BODY(I) + BODY(J))
      SEMI = 1.0/A1
      ECC2 = (1.0 - RIJ/SEMI)**2 + RDOT**2/(SEMI*(BODY(I) + BODY(J)))
      ECC = SQRT(ECC2)
*
      RETURN
*
      END
