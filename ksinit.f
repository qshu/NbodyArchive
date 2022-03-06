      SUBROUTINE KSINIT
*
*
*       Initialization of KS regularization.
*       ------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  Q(3),RDOT(3),A1(3,4)
*
*
*       Set new global indices of the components and current pair index.
      ICOMP = 2*NPAIRS - 1
      JCOMP = ICOMP + 1
      IPAIR = NPAIRS
*
*       Add body #N in case the only neighbour was removed in KSREG.
      IF (LIST(1,ICOMP).EQ.0) THEN
          LIST(2,ICOMP) = N
          LIST(2,JCOMP) = N
          LIST(1,ICOMP) = 1
          LIST(1,JCOMP) = 1
      END IF
*
*       Specify mass, neighbour radius & name for new c.m.
      BODY(NTOT) = BODY(ICOMP) + BODY(JCOMP)
      RS(NTOT) = RS(ICOMP)
      NAME(NTOT) = NZERO + NAME(ICOMP)
*
*       Define c.m. coordinates & velocities and set XDOT for components.
      DO 10 K = 1,3
          X(K,NTOT) = (BODY(ICOMP)*X(K,ICOMP) + BODY(JCOMP)*X(K,JCOMP))/
     &                                                        BODY(NTOT)
          X0DOT(K,NTOT) = (BODY(ICOMP)*X0DOT(K,ICOMP) + BODY(JCOMP)*
     &                                        X0DOT(K,JCOMP))/BODY(NTOT)
          XDOT(K,NTOT) = X0DOT(K,NTOT)
          XDOT(K,ICOMP) = X0DOT(K,ICOMP)
          XDOT(K,JCOMP) = X0DOT(K,JCOMP)
   10 CONTINUE
*
*       Obtain force polynomial for c.m. with components ICOMP & JCOMP.
      NNB = LIST(1,ICOMP)
*
*       Predict current coordinates & velocities for the neighbours.
      CALL XVPRED(ICOMP,NNB)
*
*       Obtain new polynomials & steps.
      CALL FPOLY1(ICOMP,JCOMP,1)
      CALL FPOLY2(NTOT,NTOT,1)
*
*       See whether to include new c.m. in NLIST.
      IF (T0(NTOT) + STEP(NTOT).LT.TLIST) THEN
          CALL NLMOD(NTOT,1)
      END IF
*
*       Skip KS initialization at merger termination (H, U & UDOT in RESET).
      IF (IPHASE.EQ.7) THEN
          CALL KSLIST(IPAIR)
          EB = -BODYM*ECLOSE
          GO TO 50
      END IF
*
*       Define relative coordinates and velocities in physical units.
      DO 20 K = 1,3
          Q(K) = X(K,ICOMP) - X(K,JCOMP)
          RDOT(K) = X0DOT(K,ICOMP) - X0DOT(K,JCOMP)
   20 CONTINUE
*
*       Introduce regularized variables using definition of 1985 paper.
      R(IPAIR) = SQRT(Q(1)**2 + Q(2)**2 + Q(3)**2)
*
*       Initialize the regularized coordinates according to sign of Q(1).
      IF (Q(1).LE.0.0D0) THEN
          U(3,IPAIR) = 0.0D0
          U(2,IPAIR) = SQRT(0.5D0*(R(IPAIR) - Q(1)))
          U(1,IPAIR) = 0.5D0*Q(2)/U(2,IPAIR)
          U(4,IPAIR) = 0.5D0*Q(3)/U(2,IPAIR)
      ELSE
          U(4,IPAIR) = 0.0D0
          U(1,IPAIR) = SQRT(0.5D0*(R(IPAIR) + Q(1)))
          U(2,IPAIR) = 0.5D0*Q(2)/U(1,IPAIR)
          U(3,IPAIR) = 0.5D0*Q(3)/U(1,IPAIR)
      END IF
*
*       Set current transformation matrix.
      A1(1,1) =  U(1,IPAIR)
      A1(1,2) = -U(2,IPAIR)
      A1(1,3) = -U(3,IPAIR)
      A1(1,4) =  U(4,IPAIR)
      A1(2,1) =  U(2,IPAIR)
      A1(2,2) =  U(1,IPAIR)
      A1(2,3) = -U(4,IPAIR)
      A1(2,4) = -U(3,IPAIR)
      A1(3,1) =  U(3,IPAIR)
      A1(3,2) =  U(4,IPAIR)
      A1(3,3) =  U(1,IPAIR)
      A1(3,4) =  U(2,IPAIR)
*
*       Form regularized velocity and set initial KS coordinates.
      DO 30 K = 1,4
          UDOT(K,IPAIR) = 0.50D0*(A1(1,K)*RDOT(1) + A1(2,K)*RDOT(2) +
     &                                              A1(3,K)*RDOT(3))
*       Note that A1(J,K) is the transpose of A1(K,J).
          U0(K,IPAIR) = U(K,IPAIR)
   30 CONTINUE
*
*       Evaluate initial binding energy per unit mass.
      H(IPAIR) = (2.0D0*(UDOT(1,IPAIR)**2 + UDOT(2,IPAIR)**2 +
     &                   UDOT(3,IPAIR)**2 + UDOT(4,IPAIR)**2) -
     &                                              BODY(NTOT))/R(IPAIR)
*
*       Select perturbing particles for relative motion at apocentre.
      EB = H(IPAIR)*BODY(ICOMP)*BODY(JCOMP)/BODY(NTOT)
*       Use twice semi-major axis for hard binary or 2*RMIN otherwise.
      IF (EB.LT.-0.5*BODYM*ECLOSE) THEN
          RAP = -BODY(NTOT)/H(IPAIR)
      ELSE
          RAP = 2.0*RMIN
      END IF
*
*       Set cubed distance for significant perturbation.
      RCRIT3 = 2.0*RAP**3/(BODY(NTOT)*GMIN)
      ITER = 0
      RJMIN3 = RS(NTOT)**3
      JMIN = IFIRST
   35 NNB = 0
*
*       Form perturber list using tidal approximation (GAMMA < GMIN).
      NNB1 = LIST(1,NTOT) + 1
      DO 40 L = 2,NNB1
          J = LIST(L,NTOT)
          RIJ2 = (X(1,NTOT) - X(1,J))**2 + (X(2,NTOT) - X(2,J))**2 +
     &                                     (X(3,NTOT) - X(3,J))**2
          RIJ3 = RIJ2*SQRT(RIJ2)
          IF (RIJ3.LT.BODY(J)*RCRIT3) THEN
              NNB = NNB + 1
              LIST(NNB+1,ICOMP) = J
              IF (RIJ3.LT.RJMIN3) THEN
                  RJMIN3 = RIJ3
                  JMIN = J
              END IF
          END IF
   40 CONTINUE
*
*       Check whether to increase the perturber sphere.
      IF (NNB.EQ.0) THEN
          RCRIT3 = 3.0*RCRIT3
          GO TO 35
      END IF
*
*       Compare perturber range & termination criterion for soft binary.
      IF (EB.GT.-0.5*BODYM*ECLOSE.AND.ITER.EQ.0) THEN
          GIMAX = (RMIN*ABS(H(IPAIR))/BODY(NTOT))**3
          GIMAX = MAX(GIMAX,GMAX)
          IF (EB.GT.0.0) GIMAX = GMAX
          GIJ = 2.0*BODY(JMIN)*R(IPAIR)**3/(BODY(NTOT)*RJMIN3)
          RIMAX = R(IPAIR)*(GIMAX/GIJ)**0.333
*       Increase perturber sphere if estimated termination exceeds 2*RMIN.
          IF (RIMAX.GT.2.0*RMIN) THEN
              RCRIT3 = RCRIT3*(RIMAX/(2.0*RMIN))**3
              ITER = 1
              GO TO 35
          END IF
      END IF
*
*       Specify the membership (adopt zero for second component).
      LIST(1,ICOMP) = NNB
   50 LIST(1,JCOMP) = 0
*       Set large step for second component to avoid NLIST detection.
      STEP(JCOMP) = 1.0E+06
*
*       Obtain polynomials for perturbed KS motion (standard case & merger).
      CALL KSPOLY(IPAIR,1)
*
*       Set 2*SEMI as termination scale for hard binary if 2*SEMI < RS/20.
      IF (EB.LT.-0.5*BODYM*ECLOSE.AND.RAP.LT.0.05*RS(NTOT)) THEN
          R0(IPAIR) = MAX(RMIN,-BODY(NTOT)/H(IPAIR))
      ELSE
          R0(IPAIR) = R(IPAIR)
      END IF
*
*       Increase regularization counters (#9 & NKSHYP for hyperbolic orbits).
      NCOUNT(8) = NCOUNT(8) + 1
      NKSREG = NKSREG + 1
      IF (H(IPAIR).GT.0.0) NCOUNT(9) = NCOUNT(9) + 1
      IF (H(IPAIR).GT.0.0) NKSHYP = NKSHYP + 1
*
      IF (KZ(10).NE.0) THEN
          RI = SQRT((X(1,NTOT) - RDENS(1))**2 +
     &              (X(2,NTOT) - RDENS(2))**2 +
     &              (X(3,NTOT) - RDENS(3))**2)
      if(rank.eq.0)then
          WRITE (6,60)  TIME, NAME(ICOMP), NAME(JCOMP), DTAU(IPAIR),
     &                  R(IPAIR), RI, H(IPAIR), IPAIR, GAMMA(IPAIR),
     &                  STEP(NTOT), LIST(1,ICOMP)
   60     FORMAT (/,' NEW KSREG    TIME =',F7.2,2I6,F12.3,1PE10.1,
     &                                   0PF7.2,F9.2,I5,F8.3,1PE10.1,I5)
      end if
      END IF
*
*       Modify the termination criterion according to value of NPAIRS.
      IF (NPAIRS.GT.KMAX - 3) GMAX = 0.8*GMAX
      IF (NPAIRS.LT.KMAX - 5.AND.GMAX.LT.0.001) GMAX = 1.2*GMAX
      IF (NPAIRS.EQ.KMAX) WRITE (6,70)  NPAIRS, TIME
   70 FORMAT (5X,'WARNING!   MAXIMUM KS PAIRS   NPAIRS TIME',I5,F8.2)
*
*       Initialize prediction variables to avoid skipping KS components.
      DO 75 KCOMP = 1,2
          JDUM = 2*NPAIRS - 2 + KCOMP
          DO 72 K = 1,3
              X0(K,JDUM) = X(K,JDUM)
              X0DOT(K,JDUM) = 0.0D0
              F(K,JDUM) = 0.0D0
              FDOT(K,JDUM) = 0.0D0
              D2(K,JDUM) = 0.0D0
              D3(K,JDUM) = 0.0D0
              D2R(K,JDUM) = 0.0D0
              D3R(K,JDUM) = 0.0D0
   72     CONTINUE
   75 CONTINUE
*
*       See whether either component has been regularized recently.
      NNB = LISTD(1) + 1
      K = 0
*       Check case of initial binary and loop over disrupted pairs.
      IF (IABS(NAME(ICOMP) - NAME(JCOMP)).EQ.1) THEN
          IF (NAME(ICOMP).LE.2*NBIN0) K = -1
      END IF
      DO 80 L = 2,NNB
          IF (NAME(ICOMP).EQ.LISTD(L).OR.NAME(JCOMP).EQ.LISTD(L)) K = -1
   80 CONTINUE
*
*       Ensure that mergers are treated as new binaries.
      IF (IPHASE.EQ.6) K = 0
*       Set flags to distinguish primordial binaries & standard KS motion.
      LIST(2,JCOMP) = K
      KSLOW(IPAIR) = 1
*
*       Check diagnostic output of new hard binary.
      IF (KZ(8).GT.0.AND.K.EQ.0) THEN
          IF (EB.GT.-0.5*BODYM*ECLOSE) GO TO 100
          SEMI = -0.5*BODY(NTOT)/H(IPAIR)
          RI = SQRT((X(1,NTOT) - RDENS(1))**2 +
     &              (X(2,NTOT) - RDENS(2))**2 +
     &              (X(3,NTOT) - RDENS(3))**2)
          IF (IPHASE.EQ.6) K = -1
          WRITE (8,90)  TIME, NAME(ICOMP), NAME(JCOMP), K, BODY(ICOMP),
     &                  BODY(JCOMP), EB, SEMI, R(IPAIR), GAMMA(IPAIR),
     &                  RI
   90     FORMAT (' NEW BINARY   T =',F7.1,'  NAME = ',2I5,I3,
     &                        '  M =',2F8.4,'  EB =',F9.4,'  A =',F7.4,
     &                          '  R =',F7.4,'  G =',F6.3,'  RI =',F5.2)
      END IF
*
  100 RETURN
*
      END
