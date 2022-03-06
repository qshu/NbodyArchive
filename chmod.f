      SUBROUTINE CHMOD(ISUB,KCASE)
*
*
*       Modification of chain member(s).
*       --------------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &            NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK,XCM(3),VCM(3)
      INTEGER  ISORT(NMX)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX,1)
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/INCOND/  X4(3,NMX),XDOT4(3,NMX)
*
*
*       First identify the dominant perturber.
      PMAX = 0.0
      NNB = LISTC(1,1)
      DO 1 L = 2,NNB+1
          J = LISTC(L,1)
          RIJ2 = (X(1,J) - X(1,ICH))**2 + (X(2,J) - X(2,ICH))**2 +
     &                                    (X(3,J) - X(3,ICH))**2
          PIJ = BODY(J)/(RIJ2*SQRT(RIJ2))
          IF (PIJ.GT.PMAX) THEN
              PMAX = PIJ
              RJMIN2 = RIJ2
              JCLOSE = J
          END IF
    1 CONTINUE
*
*       Form the scalar product R*V for sign of radial velocity.
      RDOT = (X(1,JCLOSE) - X(1,ICH))*(XDOT(1,JCLOSE) - XDOT(1,ICH)) +
     &       (X(2,JCLOSE) - X(2,ICH))*(XDOT(2,JCLOSE) - XDOT(2,ICH)) +
     &       (X(3,JCLOSE) - X(3,ICH))*(XDOT(3,JCLOSE) - XDOT(3,ICH))
*
*       Check acceptance criterion (RSUM + RIJ < RMIN & RDOT < 0).
      IF (RSUM + SQRT(RJMIN2).LT.RMIN.AND.RDOT.LT.0.0) THEN
          RSUM = RSUM + SQRT(RJMIN2)
          IF (KZ(30).GT.2) THEN
              WRITE (6,2)  JCLOSE, NAME(JCLOSE), SQRT(RJMIN2), RSUM
    2         FORMAT (' CHMOD:   JCLOSE NMJ RJMIN RSUM ',2I5,1P,2E9.1)
          END IF
          CALL ABSORB(ISUB)
*       Activate indicator for new chain treatment.
          KCASE = 1
          GO TO 50
      END IF
*
*       Place index of the smallest INVERSE distance in ISORT(1).
      CALL HPSORT(NN-1,RINV,ISORT)
*
*       Determine index of escaper candidate (single star or binary).
      KCASE = 0
      JESC = 0
*       Distinguish two cases of each type (beginning & end of chain).
      IF (ISORT(1).EQ.1) THEN
          IESC = INAME(1)
          KCASE = 1
      ELSE IF (ISORT(1).EQ.NN-1) THEN
          IESC = INAME(NN)
          KCASE = 1
*       Check for possible binary escaper (NN = 3 implies single escaper).
      ELSE IF (ISORT(1).EQ.2) THEN
          IESC = INAME(1)
          JESC = INAME(2)
          IBIN = 1
          KCASE = 2
*       Ignore wrong identification if NN = 4 (termination on binary escape).
      ELSE IF (ISORT(1).EQ.NN-2) THEN
          IESC = INAME(NN-1)
          JESC = INAME(NN)
          IBIN = NN - 1
          KCASE = 2
      END IF
*
*       Include safety test for abnormal configuration.
      IF (KCASE.EQ.0) GO TO 60
*
      IF (KZ(30).GT.2) THEN
          WRITE (6,4)  IESC, JESC,NSTEP1,ISORT(1),(1.0/RINV(K),K=1,NN-1)
    4     FORMAT (' CHMOD:   IESC JESC # ISORT1 R ',2I3,I5,I3,1P,5E9.1)
      END IF
*
*       Copy chain variables to standard form.
      LK = 0
      DO 10 L = 1,NCH
          DO 5 K = 1,3
              LK = LK + 1
              X4(K,L) = XCH(LK)
              XDOT4(K,L) = VCH(LK)
    5     CONTINUE
   10 CONTINUE
*
*       First check case of escaping binary (JESC > 0 & RB < RGRAV).
      IF (JESC.GT.0) THEN
          RB = 1.0/RINV(IBIN)
          IF (RB.GT.RGRAV) GO TO 30
      ELSE
          GO TO 30
      END IF
*
*       Form coordinates & velocities of escaping binary (local c.m. frame).
      BCM = BODYC(IESC) + BODYC(JESC)
      RI2 = 0.0
      RDOT = 0.0
      DO 25 K = 1,3
          XCM(K) = (BODYC(IESC)*X4(K,IESC) + BODYC(JESC)*X4(K,JESC))/BCM
          VCM(K) = (BODYC(IESC)*XDOT4(K,IESC) +
     &              BODYC(JESC)*XDOT4(K,JESC))/BCM
          RI2 = RI2 + XCM(K)**2
          RDOT = RDOT + XCM(K)*VCM(K)
   25 CONTINUE
*
*       Convert to relative distance & radial velocity w.r. to inner part.
      FAC = BODY(ICH)/(BODY(ICH) - BCM)
      RI = SQRT(RI2)
      RDOT = FAC*RDOT/RI
      RI = FAC*RI
      RESC = MAX(3.0*RGRAV,0.5*RSUM)
*
*       Employ parabolic escape criterion (or terminate if RI > RMIN).
      IF (RI.GT.RESC.AND.RDOT.GT.0.0) THEN
          IF (RDOT**2.LT.2.0*BODY(ICH)/RI) THEN
              KCASE = 0
              IF (RI.GT.RMIN) KCASE = -1
              GO TO 60
          ELSE
              WRITE (6,28)  IESC, JESC, RI, RDOT**2, 2.0*BODY(ICH)/RI
   28         FORMAT (' CHAIN BINARY ESCAPE:   IESC JESC RI RDOT2',
     &                                       ' 2*M/R ',2I3,1P,3E9.1)
*       Enforce termination if NCH <= 4 (chain initialization not possible).
              IF (NCH.LE.4) THEN
                  KCASE = -1
                  GO TO 50
              END IF
              GO TO 40
          END IF
      ELSE
          KCASE = 0
          GO TO 60
      END IF
*
*       Form relative distance and radial velocity for single particle.
   30 RI = SQRT(X4(1,IESC)**2 + X4(2,IESC)**2 + X4(3,IESC)**2)
      RDOT = X4(1,IESC)*XDOT4(1,IESC) + X4(2,IESC)*XDOT4(2,IESC) +
     &                                  X4(3,IESC)*XDOT4(3,IESC)
      FAC = BODY(ICH)/(BODY(ICH) - BODYC(IESC))
      RDOT = FAC*RDOT/RI
      RI = FAC*RI
*
*       Check approximate escape criterion (or terminate if RI > RMIN).
      IF (RI.GT.2.0*RGRAV.AND.RDOT.GT.0.0) THEN
          IF (RDOT**2.LT.2.0*BODY(ICH)/RI) THEN
              KCASE = 0
              IF (RI.GT.RMIN) KCASE = -1
              GO TO 60
          END IF
          IF (KZ(30).GT.1) THEN
              WRITE (6,35)  IESC, RI, RDOT**2, 2.0*BODY(ICH)/RI
   35         FORMAT (' CHAIN SINGLE ESCAPE:   IESC RI RDOT2 2*M/R ',
     &                                         I3,1P,3E9.1)
          END IF
*      Ensure single body is removed in case of wide binary.
          IF (JESC.GT.0) KCASE = 1
      ELSE
          KCASE = 0
          GO TO 60
      END IF
*
*       Reduce chain membership (NCH > 3) or specify termination.
   40 IF (NCH.GT.3) THEN
*       Subtract largest chain distance from system size (ignore binary).
          RSUM = RSUM - 1.0/RINV(ISORT(1))
          CALL REDUCE(IESC,JESC,ISUB)
      ELSE
          KCASE = -1
      END IF
*
*       Set phase indicator < 0 to ensure new NLIST in routine INTGRT.
   50 IPHASE = -1
*
   60 RETURN
*
      END
