      SUBROUTINE BINDAT
*
*
*       Binary data bank.
*       -----------------
*
      INCLUDE 'common6.h'
      REAL  EB(KMAX),ECC(KMAX),RCM(KMAX),ECM(KMAX),AS(30)
*
*
*       Form binding energy and central distance for each KS pair.
      ZMBIN = 0.0
      DO 10 JPAIR = 1,NPAIRS
          J2 = 2*JPAIR
          J1 = J2 - 1
          ICM = N + JPAIR
          ZMBIN = ZMBIN + BODY(ICM)
*       Avoid division by zero for merged ghost binary.
          IF (BODY(J1).GT.0.0) THEN
              EB(JPAIR) = BODY(J1)*BODY(J2)*H(JPAIR)/
     &                                             (BODY(J1) + BODY(J2))
              SEMI = -0.5*BODY(ICM)/H(JPAIR)
              ECC2 = (1.0 - R(JPAIR)/SEMI)**2 +
     &                                  TDOT2(JPAIR)**2/(BODY(ICM)*SEMI)
              ECC(JPAIR) = SQRT(ECC2)
          ELSE
              EB(JPAIR) = 0.0
              ECC(JPAIR) = 0.0
          END IF
          RCM(JPAIR) = SQRT((X(1,ICM) - RDENS(1))**2 +
     &                      (X(2,ICM) - RDENS(2))**2 +
     &                      (X(3,ICM) - RDENS(3))**2)
*       Obtain binding energy (per unit mass) of c.m. motion.
          VJ2 = XDOT(1,ICM)**2 + XDOT(2,ICM)**2 + XDOT(3,ICM)**2
          POTJ = 0.0
          DO 5 J = IFIRST,NTOT
              IF (J.EQ.ICM) GO TO 5
              RIJ2 = (X(1,ICM) - X(1,J))**2 + (X(2,ICM) - X(2,J))**2 +
     &                                        (X(3,ICM) - X(3,J))**2
              POTJ = POTJ + BODY(J)/SQRT(RIJ2)
    5     CONTINUE
          ECM(JPAIR) = 0.5*VJ2 - POTJ
*       Check for external tidal field (note that HT includes mass).
          IF (KZ(14).NE.0) THEN
              CALL XTRNLV(ICM,ICM)
              ECM(JPAIR) = ECM(JPAIR) + HT/(BODY(ICM) + 1.0E-20)
          END IF
   10 CONTINUE
*
*       Copy relevant binary diagnostics to single precision.
      AS(1) = TIME
      AS(2) = RSCALE
      AS(3) = RTIDE
      AS(4) = RC
      AS(5) = 0.0
      AS(6) = 0.0
      AS(7) = 0.0
      DO 20 K = 1,10
          AS(K+7) = E(K)
   20 CONTINUE
      AS(18) = SBCOLL
      AS(19) = BBCOLL
      AS(20) = ZKIN
      AS(21) = POT
      AS(22) = EBIN0
      AS(23) = EBIN
      AS(24) = ESUB
      AS(25) = EMERGE
      AS(26) = BE(3)
      AS(27) = ZMASS
      AS(28) = ZMBIN
      AS(29) = CHCOLL
      AS(30) = 0.0
*
*       Write formatted data bank on unit 9.
      OPEN (UNIT=9,FILE='OUT9')
      WRITE (9,30)  NPAIRS, MODEL, NRUN, N, NC, NMERGE, (AS(K),K=1,7)
   30 FORMAT (3I4,I6,2I4,2X,F7.1,2F7.2,F7.3,3F9.4)
      WRITE (9,35)  (AS(K),K=8,17)
   35 FORMAT (10F11.6)
      WRITE (9,40)  (AS(K),K=18,30)
   40 FORMAT (13F10.5)
      DO 50 JPAIR = 1,NPAIRS
          J1 = 2*JPAIR - 1
          J2 = 2*JPAIR
          WRITE (9,45)  EB(JPAIR), ECC(JPAIR), ECM(JPAIR), RCM(JPAIR),
     &                  NAME(J1), NAME(J2)
   45     FORMAT (F8.5,F7.3,F8.3,F7.3,2I6)
   50 CONTINUE
      CLOSE (UNIT=9)
*
      RETURN
*
      END
