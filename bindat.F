      SUBROUTINE BINDAT
*
*
*       Binary data bank.
*       -----------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX)
      REAL*4  EB(KMAX),ECC(KMAX),RCM(KMAX),ECM(KMAX),PB(KMAX),AS(30)
*
*
*       Form binding energy and central distance for each KS pair.
      ZMBIN = 0.0
      DO 10 JPAIR = 1,NPAIRS
          J2 = 2*JPAIR
          J1 = J2 - 1
          ICM = N + JPAIR
          ZMBIN = ZMBIN + BODY(ICM)
          BODYCM = BODY(ICM)
*       Determine merger & ghost index for negative c.m. name (skip ghost).
          IF (NAME(ICM).LT.0.AND.BODY(ICM).GT.0.0) THEN
              CALL FINDJ(J1,J,IMERGE)
*       Employ actual masses and two-body distance for energy & eccentricity.
              BODYCM = CM(1,IMERGE) + CM(2,IMERGE)
              EB(JPAIR) = CM(1,IMERGE)*CM(2,IMERGE)*HM(IMERGE)/BODYCM
              SEMI = -0.5*BODYCM/HM(IMERGE)
              RJ = SQRT(XREL(1,IMERGE)**2 + XREL(2,IMERGE)**2 +
     &                                      XREL(3,IMERGE)**2)
*       Assume that merged binary is near apo or peri (hence ignore TDOT2).
              ECC2 = (1.0 - RJ/SEMI)**2
*       Include separate diagnostics for the hierarchy (inner comps J1 & J).
              SEMI1 = -0.5*BODY(ICM)/H(JPAIR)
              ECC1 = (1.0 - R(JPAIR)/SEMI1)**2 +
     &                                 TDOT2(JPAIR)**2/(BODY(ICM)*SEMI1)
              E0 = SQRT(ECC2)
              E1 = SQRT(ECC1)
              Q = BODY(J2)/BODYCM
              XFAC = (1.0 + Q)*(1.0 + E1)/SQRT(ABS(1.0 - E1))
              FE = 2.8
              PCR = FE*XFAC**0.4*SEMI
              PM = SEMI1*(1.0 - E1)/PCR
              IF (J.LT.0) J = J1
              RM = SEMI*(1.0 - E0)/MAX(RADIUS(J1),RADIUS(J),1.0D-20)
              RM = MIN(RM,99.9D0)
              P0 = DAYS*SEMI*SQRT(ABS(SEMI)/BODYCM)
              P1 = DAYS*SEMI1*SQRT(ABS(SEMI1)/BODY(ICM))
              if(rank.eq.0)then
              WRITE (84,2) TTOT, NAME(J1), NAME(J), KSTAR(J1), KSTAR(J),
     &                     KSTARM(IMERGE), E0, E1, PM, RM, P0, P1, SEMI1
    2         FORMAT (' BINDAT:    T NM K* E0 E1 PM/PC PM0/R* P0 P1 A1',
     &                             F8.1,2I5,3I4,2F7.3,2F6.1,1P,3E9.1)
              CALL FLUSH(84)
              end if
          ELSE IF (BODY(J1).GT.0.0D0) THEN
*       Form binding energy and eccentricity for standard case.
              EB(JPAIR) = BODY(J1)*BODY(J2)*H(JPAIR)/
     &                                             (BODY(J1) + BODY(J2))
              SEMI = -0.5*BODY(ICM)/H(JPAIR)
              ECC2 = (1.0 - R(JPAIR)/SEMI)**2 +
     &                                  TDOT2(JPAIR)**2/(BODY(ICM)*SEMI)
          ELSE
              IM = 0
*       Search merger table to identify corresponding index of c.m. name.
              DO 5 K = 1,NMERGE
                  IF (NAMEG(K).EQ.NAME(ICM)) THEN
                      IM = K
                  END IF
    5         CONTINUE
              BODYJ1 = CM(3,IM)
              BODYJ2 = CM(4,IM)
              BODYCM = BODYJ1 + BODYJ2
              BODYCM = MAX(BODYCM,1.0D-10)
              EB(JPAIR) = BODYJ1*BODYJ2*H(JPAIR)/BODYCM
              SEMI = -0.5*BODYCM/H(JPAIR)
              ECC2 = (1.0 - SEMI/R(JPAIR))**2
          END IF
          ECC(JPAIR) = SQRT(ECC2)
          EB(JPAIR) = MAX(EB(JPAIR),-9.99999)
          PB(JPAIR) = DAYS*SEMI*SQRT(ABS(SEMI)/BODYCM)
          PB(JPAIR) = MIN(PB(JPAIR),99999.9)
*       Obtain binding energy (per unit mass) of c.m. motion.
          VJ2 = XDOT(1,ICM)**2 + XDOT(2,ICM)**2 + XDOT(3,ICM)**2
          IF (BODY(ICM).EQ.0.0D0) VJ2 = 0.0
          POTJ = 0.0
          DO 9 J = IFIRST,NTOT
              IF (J.EQ.ICM.or.nomass(j).eq.1) GO TO 9
              RIJ2 = (X(1,ICM) - X(1,J))**2 + (X(2,ICM) - X(2,J))**2 +
     &                                        (X(3,ICM) - X(3,J))**2 
              POTJ = POTJ + BODY(J)/SQRT(RIJ2)
    9     CONTINUE
          ECM(JPAIR) = 0.5*VJ2 - POTJ
*       Check for external tidal field (note that HT includes mass).
          IF (KZ(14).GT.0) THEN
              CALL XTRNLV(ICM,ICM)
              ECM(JPAIR) = ECM(JPAIR) + HT/(BODY(ICM) + 1.0E-20)
          END IF
          RCM(JPAIR) = SQRT((X(1,ICM) - RDENS(1))**2 +
     &                      (X(2,ICM) - RDENS(2))**2 +
     &                      (X(3,ICM) - RDENS(3))**2)
          RCM(JPAIR) = MIN(RCM(JPAIR),99.9)
   10 CONTINUE
*
*       Copy relevant binary diagnostics to single precision.
      AS(1) = TTOT
      AS(2) = RSCALE
      AS(3) = RTIDE
      AS(4) = RC
      AS(5) = TPHYS
      AS(6) = -1.5*(TIDAL(1)*ZMASS**2)**0.3333
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
      AS(30) = ECOLL
*
*       Write formatted data bank on unit 9.
*
#ifdef PARALLEL
      if(rank.eq.0)then
#endif
      WRITE (9,30)  NPAIRS, MODEL, NRUN, N, NC, NMERGE, (AS(K),K=1,7)
   30 FORMAT (3I4,I6,2I4,2X,F7.1,2F7.2,F7.3,3F9.4)
      WRITE (9,35)  (AS(K),K=8,17)
   35 FORMAT (10F11.6)
      WRITE (9,40)  (AS(K),K=18,30)
   40 FORMAT (13F10.5)
*
      DO 50 JPAIR = 1,NPAIRS
          J1 = 2*JPAIR - 1
          J2 = 2*JPAIR
          KCM = KSTAR(N+JPAIR)
          IF (NAME(N+JPAIR).LT.0) THEN
              KCM = -10
          END IF
          WRITE (9,45)  EB(JPAIR), ECC(JPAIR), ECM(JPAIR), RCM(JPAIR),
     &                  BODY(J1)*ZMBAR, BODY(J2)*ZMBAR, PB(JPAIR),
     &                  NAME(J1), NAME(J2), KSTAR(J1), KSTAR(J2), KCM
   45     FORMAT (F8.5,F7.3,F7.2,F6.2,2F5.1,F8.1,2I6,3I4)
   50 CONTINUE
      CALL FLUSH(9)
#ifdef PARALLEL
      end if
#endif
*
      RETURN
*
      END
*
      SUBROUTINE PLNDAT
*
*
*       Planet Binary data bank.
*       ------------------------
*
      INCLUDE 'common6.h'
      REAL*4  EB(KMAX),ECC(KMAX),RCM(KMAX),ECM(KMAX),PB(KMAX),AS(30)
      LOGICAL  FIRST
      SAVE  FIRST
      DATA  FIRST /.TRUE./
      REAL*4 SEMI0(KMAX),ECC0(KMAX),BODCM(KMAX),SEMIA(KMAX)
*
*
*       Form binding energy and central distance for each KS pair.
      ZMBIN = 0.0
      DO 10 JPAIR = 1,NPAIRS
          J2 = 2*JPAIR
          J1 = J2 - 1
          ICM = N + JPAIR
          ZMBIN = ZMBIN + BODY(ICM)
          BODYCM = BODY(J1) + BODY(J2)
*       Form binding energy and eccentricity for standard case.
              SEMI = -0.5*BODY(ICM)/H(JPAIR)
              SEMIA(JPAIR) = SEMI
              ECC2 = (1.0 - R(JPAIR)/SEMI)**2 +
     &                                  TDOT2(JPAIR)**2/(BODY(ICM)*SEMI)
          ECC(JPAIR) = SQRT(ECC2)
          PB(JPAIR) = DAYS*SEMI*SQRT(ABS(SEMI)/BODYCM)
          PB(JPAIR) = MIN(PB(JPAIR),99999.9)
*
          IF(FIRST)THEN
          SEMI0(JPAIR) = SEMI
          ECC0(JPAIR) = ECC(JPAIR)
          BODCM(JPAIR) = BODYCM
          END IF
*
*       Obtain binding energy (per unit mass) of c.m. motion.
          VJ2 = XDOT(1,ICM)**2 + XDOT(2,ICM)**2 + XDOT(3,ICM)**2
          POTJ = 0.0
          DO 9 J = IFIRST,NTOT
              IF (J.EQ.ICM.or.nomass(j).eq.1) GO TO 9
              RIJ2 = (X(1,ICM) - X(1,J))**2 + (X(2,ICM) - X(2,J))**2 +
     &                                        (X(3,ICM) - X(3,J))**2 
              POTJ = POTJ + BODY(J)/SQRT(RIJ2)
    9     CONTINUE
          ECM(JPAIR) = 0.5*VJ2 - POTJ
*       Check for external tidal field (note that HT includes mass).
          IF (KZ(14).GT.0) THEN
              CALL XTRNLV(ICM,ICM)
              ECM(JPAIR) = ECM(JPAIR) + HT/(BODY(ICM) + 1.0E-20)
          END IF
          RCM(JPAIR) = SQRT((X(1,ICM) - RDENS(1))**2 +
     &                      (X(2,ICM) - RDENS(2))**2 +
     &                      (X(3,ICM) - RDENS(3))**2)
          RCM(JPAIR) = MIN(RCM(JPAIR),99.9)
   10 CONTINUE
*
*       Copy relevant binary diagnostics to single precision.
      AS(1) = TTOT
      AS(2) = RSCALE
      AS(3) = RTIDE
      AS(4) = RC
      AS(5) = TPHYS
      AS(6) = -1.5*(TIDAL(1)*ZMASS**2)**0.3333
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
      AS(30) = ECOLL
*
*       Write formatted data bank on unit 39.
#ifdef PARALLEL
      if(rank.eq.0)then
#endif
*
      WRITE (39,30)  NPAIRS, MODEL, NRUN, N, NC, NMERGE, (AS(K),K=1,7)
   30 FORMAT (3I4,I6,2I4,2X,F7.1,2F7.2,F7.3,3F9.4)
      WRITE (39,35)  (AS(K),K=8,17)
   35 FORMAT (10F11.6)
      WRITE (39,40)  (AS(K),K=18,30)
   40 FORMAT (13F10.5)
*
      DO 50 JPAIR = 1,NPAIRS
          J1 = 2*JPAIR - 1
          J2 = 2*JPAIR
          KCM = KSTAR(N+JPAIR)
          IF (NAME(N+JPAIR).LT.0) THEN
              KCM = -10
          END IF
       WRITE(39,45)SEMIA(JPAIR)*rau,ECC(JPAIR), ECM(JPAIR), RCM(JPAIR),
     &                  BODY(J1)*ZMBAR, BODY(J2)*ZMBAR, PB(JPAIR),
     &                  NAME(J1), NAME(J2), KSTAR(J1), KSTAR(J2), KCM,
     &                  SEMI0(JPAIR)*rau,ECC0(JPAIR),BODCM(JPAIR)
   45     FORMAT (1X,1P,7D13.5,2I6,3I4,3D13.5)
   50 CONTINUE
      CALL FLUSH(39)
#ifdef PARALLEL
      end if
#endif
*
      RETURN
*
      END
*
      SUBROUTINE ELMNTS(KK)
*
*
*       Planet orbital data computation and save
*       ----------------------------------------
*
      INCLUDE 'common6.h'
*         Variables for Planetary Data
      REAL*8  SEMIPL(NMAX),ECCPLN(NMAX),XIN(NMAX),OMEG(NMAX),
     & EBPLN(NMAX),ECMPLN(NMAX),SEMIOL(NMAX),ECCOLD(NMAX),
     & XINOLD(NMAX),OMEGOL(NMAX),XPERIH(NMAX),XPEROL(NMAX),
     & XRELPL(3),VRELPL(3),XLI(3),EI(3),XNI(3),XEZ(3),XNE(3)
      INTEGER JKL(NMAX)
      DATA XEZ/0.0,0.0,1.0/
*
          IPL = NAME(KK)
          IPAIR = KVEC(KK)
          ICOMP = 2*IPAIR - 1
          JCOMP = 2*IPAIR
          I = N + IPAIR
*      Distance to closest neighbour
           NNB = LIST(1,I)
           RIJMIN = 1.D30
          DO 555 L = 2,NNB
           J = LIST(L,I)
           RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &            (X(3,I) - X(3,J))**2
           IF(RIJ2.LT.RIJMIN**2)THEN
           RIJMIN = DSQRT(RIJ2)
           JKLMIN = J
           END IF
 555       CONTINUE
*
          RI = R(IPAIR)
          HI = H(IPAIR)
          RI2 = 0.0
          VI2 = 0.0
          RVI = 0.0
          DO 333 K = 1,3
              XRELPL(K) = X(K,ICOMP) - X(K,JCOMP)
              VRELPL(K) = XDOT(K,ICOMP) - XDOT(K,JCOMP)
              RI2 = RI2 + XRELPL(K)**2
              VI2 = VI2 + VRELPL(K)**2
              RVI = RVI + XRELPL(K)*VRELPL(K)
  333     CONTINUE
*
          RPLI = DSQRT(X(1,I)**2 + X(2,I)**2 + X(3,I)**2)
          VPLI = DSQRT(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2)
*
*       Construct Runge-Lenz vector (Heggie & Rasio, 1995, IAU174, Eq.(5)).
          EI2 = 0.0
          XLI2 = 0.0
          CHECK = 0.0
          DO 6 K = 1,3
              KP1 = 1 + MOD(K,3)
              KP2 = 1 + MOD(K+1,3)
              EI(K) = (VI2*XRELPL(K) - RVI*VRELPL(K))/BODY(I) -
     &                                        XRELPL(K)/SQRT(RI2)
              XLI(K) = XRELPL(KP1)*VRELPL(KP2)-XRELPL(KP2)*VRELPL(KP1)
              EI2 = EI2 + EI(K)**2
              XLI2 = XLI2 + XLI(K)**2
              CHECK = CHECK + EI(K)*XLI(K)
    6     CONTINUE
*       Define vector parallel to line of nodes (node vector)
          DO 7 K = 1,3
              KP1 = 1 + MOD(K,3)
              KP2 = 1 + MOD(K+1,3)
              XNI(K) = XEZ(KP1)*XLI(KP2)-XEZ(KP2)*XLI(KP1)
    7     CONTINUE
*       Determine angle of node by using (n.e) (nxe) n node vector
*
          DO 71 K = 1,3
              KP1 = 1 + MOD(K,3)
              KP2 = 1 + MOD(K+1,3)
              XNE(K) = XNI(KP1)*EI(KP2) - XNI(KP2)*EI(KP1)
   71     CONTINUE
          XAX = 0.D0
          XAY = 0.D0
          DO 8 K = 1,3
              XAY = XAY + XNE(K)*XLI(K)
              XAX = XAX + XNI(K)*EI(K)
    8     CONTINUE
          XAY = XAY/DSQRT(XLI2)
*
          IF (XLI(3).GT.0.D0) THEN
              XINCL = ASIN(DSQRT((XLI(1)**2+XLI(2)**2)/XLI2))
              ELSE
              XINCL = TWOPI/2.D0-ASIN(DSQRT((XLI(1)**2+XLI(2)**2)/XLI2))
          END IF
*       Only defined if orbit not in plane
          IF(DABS(DSIN(XINCL)).GT.1.D-6)THEN
*
          XNX = XNI(1)/DSQRT(XLI2)/DABS(DSIN(XINCL))
          XNY = XNI(2)/DSQRT(XLI2)/DABS(DSIN(XINCL))
          XAX = XAX/DSQRT(EI2*XLI2)/DABS(DSIN(XINCL))
          XAY = XAY/DSQRT(EI2*XLI2)/DABS(DSIN(XINCL))
*       angle of line of node XOMEG, perihelion angle XPERI
          IF(XNX.GT.0.D0)THEN
              IF(XNY.GT.0.D0)THEN
                 XOMEG = ACOS(XNX)
              ELSE
                 XOMEG = TWOPI - ACOS(XNX)
              END IF
          ELSE
              IF(XNY.GT.0.D0)THEN
                 XOMEG = TWOPI/2.D0 - ACOS(-XNX)
              ELSE
                 XOMEG = TWOPI/2.D0 + ACOS(-XNX)
              END IF
          END IF
*
*
          IF(XAX.GT.0.D0)THEN
             IF(XAY.GT.0.D0)THEN
                 XPERI = ACOS(XAX)
              ELSE
                 XPERI = TWOPI - ACOS(XAX)
              END IF
          ELSE
              IF(XAY.GT.0.D0)THEN
                 XPERI = TWOPI/2.D0 - ACOS(-XAX)
              ELSE
                 XPERI = TWOPI/2.D0 + ACOS(-XAX)
              END IF
          END IF

          ELSE
*
          XOMEG = -1.D0
          XPERI = -1.D0
*
          END IF
*
          SEMI = -0.5*BODY(I)/HI
          ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
          SEMIPL(IPL) = SEMI
          ECCPLN(IPL) = DSQRT(ECC2)
          XIN(IPL) = XINCL
          OMEG(IPL) = XOMEG
          XPERIH(IPL) = XPERI
          EBPLN(IPL) = HI
          ECMPLN(IPL) = BODY(I)*(XDOT(1,I)**2
     &            + XDOT(2,I)**2 + XDOT(3,I)**2)*5.D-1
*
          IF (TIME.GT.0.D0) THEN
              DAU = (SEMIPL(IPL) - SEMIOL(IPL))*RAU
              DEL = (ECCPLN(IPL) - ECCOLD(IPL))
              DINC = (XIN(IPL) - XINOLD(IPL))
              DOMEG = (OMEG(IPL) - OMEGOL(IPL))
              DPERI = (XPERIH(IPL) - XPEROL(IPL))
*
              WRITE(77,777)IPL,TIME,BODY(I),SEMIPL(IPL)*RAU,
     *                     EBPLN(IPL),ECCPLN(IPL),ECMPLN(IPL),DAU,DEL,
     *                     XIN(IPL),OMEG(IPL),XPERIH(IPL),
     *                     DINC,DOMEG,DPERI,RPLI,VPLI,RIJMIN,JKLMIN
 777              FORMAT(1X,I8,1P,17(D12.5,1X),I6)
          END IF
*
          SEMIOL(IPL) = SEMIPL(IPL)
          ECCOLD(IPL) = ECCPLN(IPL)
          XINOLD(IPL) = XIN(IPL)
          OMEGOL(IPL) = OMEG(IPL)
          XPEROL(IPL) = XPERIH(IPL)
*
      RETURN
*
      END



