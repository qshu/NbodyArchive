      SUBROUTINE KSTERM
*
*
*       Termination of KS regularization.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX)
      REAL*8  SAVE(9)
*
*
*       Copy pair index from COMMON.
      IPAIR = KSPAIR
*       Define index of first component & corresponding c.m.
      I1 = 2*IPAIR - 1
      ICM = N + IPAIR
*
*       Prepare termination at block time (KS, triple, quad, chain or merge).
      IF (TIME.LT.TBLOCK.AND.IPHASE.LT.10) THEN
          TIME0 = MIN(T0(ICM) + STEP(ICM),TBLOCK)
    1     DT0 = TIME0 - T0(I1)
*       Integrate up to current block time in case interval is too large.
          IF (DT0.GT.STEP(I1)) THEN
              TIME = T0(I1) + STEP(I1)
              CALL KSINT(I1)
              DTU = DTAU(IPAIR)
              STEP(I1) = ((ONE6*TDOT3(IPAIR)*DTU + 0.5*TDOT2(IPAIR))*DTU
     &                                                   + R(IPAIR))*DTU
              GO TO 1
          END IF
*
*       Determine the last regularized step by Newton-Raphson iteration.
          DTU = DT0/R(IPAIR)
          ITER = 0
    2     Y0 = DT0 - ((ONE6*TDOT3(IPAIR)*DTU + 0.5*TDOT2(IPAIR))*DTU +
     &                                                     R(IPAIR))*DTU
          YPR = -((0.5*TDOT3(IPAIR)*DTU + TDOT2(IPAIR))*DTU + R(IPAIR))
          DTU = DTU - Y0/YPR
          DT1 = ((ONE6*TDOT3(IPAIR)*DTU + 0.5*TDOT2(IPAIR))*DTU +
     &                                                     R(IPAIR))*DTU
          ITER = ITER + 1
          IF (ABS(DT0 - DT1).GT.1.0E-06*STEP(I1).AND.ITER.LT.20) GO TO 2
*
*       Advance the KS solution to next block time and terminate at TIME0.
          DTAU(IPAIR) = DTU
          STEP(I1) = DT1
          TIME = T0(I1) + DT1
          CALL KSINT(I1)
          TIME = TIME0
*
*       Predict X & XDOT for ICM & JCOMP (note TIME = TBLOCK if second call).
          CALL XVPRED(ICM,0)
          IF (JCOMP.GE.IFIRST) THEN
              CALL XVPRED(JCOMP,0)
          END IF
*
*       Retain final KS variables for explicit restart at merge termination.
          IF (IPHASE.EQ.6) THEN
              HM(NMERGE) = H(IPAIR)
              DO 3 K = 1,4
                  UM(K,NMERGE) = U(K,IPAIR)
                  UMDOT(K,NMERGE) = UDOT(K,IPAIR)
    3         CONTINUE
          END IF
      END IF
*
*       Form square regularized velocity for the explicit binding energy.
      UPR2 = 0.0
      DO 5 K = 1,4
          UPR2 = UPR2 + UDOT(K,IPAIR)**2
    5 CONTINUE
*
*       Form KS scaling factors from energy and angular momentum relation.
      A1 = 0.25D0*BODY(ICM)/UPR2
*       Solve for C1 from H = (2*U'*U'*C1**2 - M)/(U*U*C2**2) with C2 = 1/C1.
      A2 = A1**2 + 0.5D0*H(IPAIR)*R(IPAIR)/UPR2
*
*       Check for undefined case (circular orbit or eccentric anomaly = 90).
      IF (A2.GT.0.0D0) THEN
          IF (A1.LT.1.0) THEN
*       Choose square root sign from eccentric anomaly (e*cos(E) = 1 - R/a).
              C1 = SQRT(A1 + SQRT(A2))
          ELSE
              C1 = SQRT(A1 - SQRT(A2))
          END IF
      ELSE
          C1 = 1.0
      END IF
*       Specify KS coordinate scaling from angular momentum conservation.
      C2 = 1.0/C1
*
*       Transform KS variables to yield the prescribed elements.
      R(IPAIR) = 0.0D0
      DO 6 K = 1,4
          U(K,IPAIR) = C2*U(K,IPAIR)
          UDOT(K,IPAIR) = C1*UDOT(K,IPAIR)
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
    6 CONTINUE
*
*       Check optional diagnostic output for disrupted new hard binary.
      IF (KZ(8).EQ.0) GO TO 10
      IF (LIST(2,I1+1).NE.0.OR.H(IPAIR).GT.0.0) GO TO 10
      IF (GAMMA(IPAIR).GT.0.5.AND.JCOMP.GT.0.OR.IPHASE.EQ.7) THEN
          IF (JCOMP.EQ.0.OR.IPHASE.EQ.7) JCOMP = I1
          K = 0
          IF (JCOMP.GT.N) THEN
              J2 = 2*(JCOMP - N)
              K = LIST(2,J2)
          END IF
          I2 = I1 + 1
          SEMI = -0.5*BODY(ICM)/H(IPAIR)
          EB = -0.5*BODY(I1)*BODY(I2)/SEMI
          RI = SQRT((X(1,ICM) - RDENS(1))**2 +
     &              (X(2,ICM) - RDENS(2))**2 +
     &              (X(3,ICM) - RDENS(3))**2)
          WRITE (8,8)  TIME, NAME(I1), NAME(I2), K, NAME(JCOMP),
     &                 BODY(JCOMP), EB, SEMI, R(IPAIR), GAMMA(IPAIR), RI
    8     FORMAT (' END BINARY   T =',F7.1,'  NAME = ',2I5,I3,I6,
     &                      '  M(J) =',F8.4,'  EB =',F10.5,'  A =',F8.5,
     &                          '  R =',F8.5,'  G =',F5.2,'  RI =',F5.2)
      END IF
*
   10 IF (KZ(11).GT.0) THEN
          RI = SQRT((X(1,ICM) - RDENS(1))**2 +
     &              (X(2,ICM) - RDENS(2))**2 +
     &              (X(3,ICM) - RDENS(3))**2)
      if(rank.eq.0)then
          WRITE (6,15)  TIME, BODY(I1), BODY(I1+1), DTAU(IPAIR),
     &                  R(IPAIR), RI, H(IPAIR), IPAIR, GAMMA(IPAIR),
     &                  STEP(I1), LIST(1,I1), LIST(1,ICM)
   15     FORMAT (/,' END KSREG    TIME =',F7.2,2F8.4,F8.3,1PE10.1,
     &                                  0PF7.2,F9.2,I5,F8.3,1PE10.1,2I5)
      end if
      END IF
*
*       Obtain global coordinates & velocities (skip mass-loss case).
      IF (IPHASE.NE.-1) THEN
          CALL RESOLV(IPAIR,2)
      END IF
*
*       Modify c.m. neighbour radius by density contrast and set new values.
      NNB = LIST(1,ICM) + 1
*       Check predicted neighbour number and form volume ratio.
      NBP = MIN(ALPHA*SQRT(FLOAT(NNB)*RS(ICM))/(RS(ICM)**2),ZNBMAX)
      NBP = MAX(NBP,INT(ZNBMIN))
      A0 = FLOAT(NBP)/FLOAT(NNB)
*       Copy all neighbours in the case of merger.
      IF (IPHASE.EQ.6) THEN
          A0 = 1.0
          NBP = NNB - 1
      END IF
      IF (RS(ICM).GT.-100.0*BODY(ICM)/H(IPAIR)) A0 = 1.0
*       Accept old c.m. values for small length scale ratio or H > 0.
      RS(I1) = RS(ICM)*A0**0.3333
      RS(I1+1) = RS(I1)
      RS2 = RS(I1)**2
*
*       Select neighbours for components inside the modified c.m. sphere.
   20 NNB1 = 1
      DO 25 L = 2,NNB
          J = LIST(L,ICM)
          RIJ2 = (X(1,ICM) - X(1,J))**2 + (X(2,ICM) - X(2,J))**2 +
     &                                    (X(3,ICM) - X(3,J))**2
*       Ensure that at least the predicted neighbour number is reached.
          IF (RIJ2.LT.RS2.OR.L + NBP.GT.NNB1 + NNB) THEN
              NNB1 = NNB1 + 1
              ILIST(NNB1) = J
          END IF
   25 CONTINUE
*
*       Check that there is space for adding dominant component later.
      IF (NNB1.GE.NNBMAX) THEN
          RS2 = 0.9*RS2
          GO TO 20
      END IF
*
*       Reduce pair index, total number & single particle index.
      NPAIRS = NPAIRS - 1
      NTOT = N + NPAIRS
      IFIRST = 2*NPAIRS + 1
*
*       Save name of components & flag for modifying LISTD in UPDATE.
      JLIST(1) = NAME(I1)
      JLIST(2) = NAME(I1+1)
      JLIST(3) = LIST(2,I1+1)
*
*       Skip adjustment of tables if last or only pair being treated.
      IF (IPAIR.EQ.NPAIRS + 1) GO TO 60
*
*       Move the second component before the first.
      DO 50 KCOMP = 2,1,-1
          I = 2*IPAIR - 2 + KCOMP
*
          DO 30 K = 1,3
              SAVE(K) = X(K,I)
              SAVE(K+3) = X0DOT(K,I)
   30     CONTINUE
*       Current velocity has been set in routine RESOLV.
          SAVE(7) = BODY(I)
          SAVE(8) = RS(I)
          SAVE(9) = RADIUS(I)
          NAMEI = NAME(I)
          LAST = 2*NPAIRS - 1 + KCOMP
*
*       Move up global variables of other components.
          DO 40 J = I,LAST
              DO 35 K = 1,3
                  X(K,J) = X(K,J+1)
*       Copy latest X & X0DOT (= 0) of single components for predictor.
                  X0(K,J) = X(K,J)
                  X0DOT(K,J) = X0DOT(K,J+1)
   35         CONTINUE
              BODY(J) = BODY(J+1)
              RS(J) = RS(J+1)
              RADIUS(J) = RADIUS(J+1)
              NAME(J) = NAME(J+1)
              STEP(J) = STEP(J+1)
              T0(J) = T0(J+1)
              K = LIST(1,J+1) + 1
              IF (K.EQ.1) K = 2
*       Transfer unmodified neighbour lists (include flag in 2nd comp).
              DO 38 L = 1,K
                  LIST(L,J) = LIST(L,J+1)
   38         CONTINUE
   40     CONTINUE
*
*       Set new component index and copy basic variables.
          I = LAST + 1
          DO 45 K = 1,3
              X(K,I) = SAVE(K)
              X0DOT(K,I) = SAVE(K+3)
              XDOT(K,I) = SAVE(K+3)
   45     CONTINUE
          BODY(I) = SAVE(7)
          RS(I) = SAVE(8)
          RADIUS(I) = SAVE(9)
          NAME(I) = NAMEI
   50 CONTINUE
*
*       Update all regularized variables.
      CALL REMOVE(IPAIR,2)
*
*       Remove old c.m. from all COMMON tables (no F & FDOT correction).
      CALL REMOVE(ICM,3)
*
*       Set new global index of first & second component.
   60 ICOMP = 2*NPAIRS + 1
      JCOMP = ICOMP + 1
*
*       Save c.m. neighbour list for routine FPOLY1/2 (may be renamed below).
      ILIST(1) = NNB1 - 1
      DO 65 L = 1,NNB1
          LIST(L,ICOMP) = ILIST(L)
   65 CONTINUE
*
*       Modify all relevant COMMON list arrays.
      CALL UPDATE(IPAIR)
*
*       Check replacing of single KS component by corresponding c.m.
   70 IF (LIST(2,ICOMP).LT.ICOMP) THEN
          J = 0.5001*(LIST(2,ICOMP) + 1) + N
          DO 80 L = 2,NNB1
              IF (L.LT.NNB1.AND.LIST(L+1,ICOMP).LT.J) THEN
                  LIST(L,ICOMP) = LIST(L+1,ICOMP)
              ELSE
                  LIST(L,ICOMP) = J
              END IF
   80     CONTINUE
*       Check again until first neighbour > ICOMP.
          GO TO 70
      END IF
*
*       Make space for dominant component and copy members to JCOMP list.
      DO 90 L = NNB1,2,-1
          LIST(L+1,ICOMP) = LIST(L,ICOMP)
          LIST(L+1,JCOMP) = LIST(L,ICOMP)
   90 CONTINUE
*
*       Set dominant component in first location and specify membership.
      LIST(2,ICOMP) = JCOMP
      LIST(2,JCOMP) = ICOMP
      LIST(1,ICOMP) = NNB1
      LIST(1,JCOMP) = NNB1
*
*       Initialize T0 & X0 for both components (note XVPRED for FPOLY 
C       and note prediction in intgrt.F - Aug.1998, P.Kroupa).
C Note in old version only the JCOMP component was treated.
      T0(JCOMP) = TIME
      T0(ICOMP) = TIME
      T0R(JCOMP) = TIME
      DO 95 K = 1,3
          X0(K,JCOMP) = X(K,JCOMP)
          X0(K,ICOMP) = X(K,ICOMP)
          X0DOT(K,JCOMP) = XDOT(K,JCOMP)
          X0DOT(K,ICOMP) = XDOT(K,ICOMP)
   95 CONTINUE
*
*       Form new force polynomials (skip triple, quad, merge & collision).
      IF (IPHASE.LT.4) THEN
*       Predict current coordinates & velocities for the neighbours.
          CALL XVPRED(ICOMP,NNB1)
*
*       Obtain new polynomials & steps.
          CALL FPOLY1(ICOMP,JCOMP,2)
          CALL FPOLY2(ICOMP,JCOMP,2)
*
*       See whether to include new c.m. or single components in NLIST.
          IF (T0(ICOMP) + STEP(ICOMP).LT.TLIST) THEN
              CALL NLMOD(ICOMP,1)
          END IF
          IF (T0(JCOMP) + STEP(JCOMP).LT.TLIST) THEN
              CALL NLMOD(JCOMP,1)
          END IF
      END IF
*
      RETURN
*
      END


