      SUBROUTINE KSREG
*
*
*       New KS regularization.
*       ----------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/FSAVE/  SAVEIT(6)
      REAL*8  SAVE(16)
      EXTERNAL RENAME
*
*
*       Ensure ICOMP < JCOMP.
      IF (ICOMP.GT.JCOMP) THEN
          ISAVE = ICOMP
          ICOMP = JCOMP
          JCOMP = ISAVE
      END IF
*
      IMAJOR = ICOMP
      MINOR = JCOMP
*       Select major component.
      IF (BODY(ICOMP).GE.BODY(JCOMP)) THEN
          IMAJOR = ICOMP
          IMINOR = JCOMP
      ELSE
          IMAJOR = JCOMP
          IMINOR = ICOMP
      END IF
*
*       Save updated regular force & derivative for new KSINIT procedure.
      IX = IMAJOR
      IF (D1R(1,IMAJOR).EQ.0.0D0.AND.ICH.GT.0) IX = ICH
*       Switch to chain c.m. after path from REDUCE or CHTERM2 directly.
*     IF (D1R(1,IMAJOR).EQ.0.0D0) THEN
*         WRITE (6,95)  ICH, IX, LIST(1,IX), D1R(1,IX)
*  95     FORMAT (' ZERO WARNING  ICH IX ',2I6,I5,1P,E10.2)
*     END IF
      DO 96 K = 1,3
          SAVEIT(K) = FR(K,IX)
          SAVEIT(K+3) = D1R(K,IX)
   96 CONTINUE
*
      LI = LIST(1,IMAJOR) + 1
      LJ = LIST(1,IMINOR) + 1
   98 IF (LI.EQ.1.AND.LJ.EQ.1) GO TO 99
      IF (LJ.GT.1.AND.(LI.EQ.1.OR.LIST(LI,IMAJOR).LT.
     &                                            LIST(LJ,IMINOR))) THEN
*     IF (LI.EQ.1.OR.LIST(LI,IMAJOR).LT.LIST(LJ,IMINOR)) THEN
*       Subtract LJ particle.
          I = LIST(LJ,IMINOR)
          LJ = LJ - 1
          IF (I.EQ.IMAJOR) GO TO 98
          A3 = -1.0
      ELSE IF (LI.GT.1.AND.(LJ.EQ.1.OR.LIST(LI,IMAJOR).GT.
     &                                            LIST(LJ,IMINOR))) THEN
*       Add LI particle.
          I = LIST(LI,IMAJOR)
          LI = LI - 1
          IF (I.EQ.IMINOR) GO TO 98
          A3 = +1.0
      ELSE
          LI = LI - 1
          LJ = LJ - 1
          GO TO 98
      END IF
*     CALL XVPRED(I,0)
*     RIJ2 = (X(1,I)-X(1,IMINOR))**2 + (X(2,I)-X(2,IMINOR))**2 +
*    &       (X(3,I)-X(3,IMINOR))**2
*     RIJ = SQRT(RIJ2)
*     RDOT = (X(1,I)-X(1,IMINOR))*(XDOT(1,I)-XDOT(1,IMINOR)) +
*    &       (X(2,I)-X(2,IMINOR))*(XDOT(2,I)-XDOT(2,IMINOR)) +
*    &       (X(3,I)-X(3,IMINOR))*(XDOT(3,I)-XDOT(3,IMINOR))
*     A3 = A3*A2*BODY(I)/RIJ**3
*     DO K = 1,3
*         SAVEIT(K) = SAVEIT(K) + A3*(X(K,I)-X(K,IMINOR))
*         SAVEIT(K+3) = SAVEIT(K+3) + A3*((XDOT(K,I)-XDOT(K,IMINOR)) +
*    &                            3.0*(X(K,I)-X(K,IMINOR))*RDOT/RIJ2)
*     END DO
      GO TO 98
   99 CONTINUE
*
*       Replace #JCOMP by arbitrary body in case it is the only neighbour.
      NNB = LIST(1,IMAJOR)
      DO L = 2,NNB+1
          J = LIST(L,IMAJOR)
          IF (J.NE.IMINOR.AND.BODY(J).GT.0d0) THEN
              GO TO 101
          ENDIF
      END DO
*
      DO J = IFIRST+2,NTOT
          IF (J.NE.IMAJOR.AND.J.NE.IMINOR.AND.BODY(J).GT.0d0) THEN
              LIST(1,IMAJOR) = LIST(1,IMAJOR) + 1
              LIST(2,IMAJOR) = J
              GO TO 101
          END IF
      END DO
  101 CONTINUE
*
*       Copy neighbour list of #ICOMP without JCOMP.
      NNB1 = 1
      DO 1 L = 1,NNB
          IF (LIST(L+1,IMAJOR).EQ.IMINOR) GO TO 1
          NNB1 = NNB1 + 1
          ILIST(NNB1) = LIST(L+1,IMAJOR)
    1 CONTINUE
      ILIST(1) = NNB1 - 1
*
*       Save basic variables for components unless in correct location.
      DO 10 KCOMP = 1,2
*       Treat the first & second component in turn.
          IF (KCOMP.EQ.1) THEN
              I = ICOMP
          ELSE
              I = JCOMP
          END IF
          J = 2*NPAIRS + KCOMP
          IF (I.EQ.J) GO TO 10
*
          DO 2 K = 1,3
              SAVE(K) = X(K,I)
              SAVE(K+3) = X0DOT(K,I)
    2     CONTINUE
          SAVE(7) = BODY(I)
          SAVE(8) = RS(I)
          SAVE(9) = RADIUS(I)
          SAVE(10) = TEV(I)
          SAVE(11) = BODY0(I)
          SAVE(12) = TEV0(I)
          SAVE(13) = EPOCH(I)
          SAVE(14) = SPIN(I)
          SAVE(15) = ZLMSTY(I)
          SAVE(16) = TNEW(I)
          NAMEI = NAME(I)
          KSI = KSTAR(I)
*
*       Exchange first & second single particle with ICOMP & JCOMP.
          DO 4 K = 1,3
              X(K,I) = X(K,J)
              X0(K,I) = X0(K,J)
              X0DOT(K,I) = X0DOT(K,J)
              XDOT(K,I) = XDOT(K,J)
              F(K,I) = F(K,J)
              FDOT(K,I) = FDOT(K,J)
              FI(K,I) = FI(K,J)
              FIDOT(K,I) = FIDOT(K,J)
              D0(K,I) = D0(K,J)
              D1(K,I) = D1(K,J)
              D2(K,I) = D2(K,J)
              D3(K,I) = D3(K,J)
              FR(K,I) = FR(K,J)
              FRDOT(K,I) = FRDOT(K,J)
              D0R(K,I) = D0R(K,J)
              D1R(K,I) = D1R(K,J)
              D2R(K,I) = D2R(K,J)
              D3R(K,I) = D3R(K,J)
              X(K,J) = SAVE(K)
              X0DOT(K,J) = SAVE(K+3)
    4     CONTINUE
*
          BODY(I) = BODY(J)
          RS(I) = RS(J)
          RADIUS(I) = RADIUS(J)
          TEV(I) = TEV(J)
          TEV0(I) = TEV0(J)
          BODY0(I) = BODY0(J)
          EPOCH(I) = EPOCH(J)
          SPIN(I) = SPIN(J)
          ZLMSTY(I) = ZLMSTY(J)
          NAME(I) = NAME(J)
          KSTAR(I) = KSTAR(J)
          STEP(I) = STEP(J)
          STEPR(I) = STEPR(J)
          T0(I) = T0(J)
          T0R(I) = T0R(J)
          K = LIST(1,J) + 1
          DO 5 L = 1,K
              LIST(L,I) = LIST(L,J)
    5     CONTINUE
          BODY(J) = SAVE(7)
          RS(J) = SAVE(8)
          RADIUS(J) = SAVE(9)
          TEV(J) = SAVE(10)
          BODY0(J) = SAVE(11)
          TEV0(J) = SAVE(12)
          EPOCH(J) = SAVE(13)
          SPIN(J) = SAVE(14)
          ZLMSTY(J) = SAVE(15)
          TNEW(J) = SAVE(16)
          NAME(J) = NAMEI
          KSTAR(J) = KSI
   10 CONTINUE
*
*       Save neighbour list of first component for RENAME & FPOLY.
      DO 15 L = 1,NNB1
          LIST(L,IFIRST) = ILIST(L)
   15 CONTINUE
*
*       Increase pair index, total number & single particle index.
      NPAIRS = NPAIRS + 1
      NTOT = N + NPAIRS
      IFIRST = 2*NPAIRS + 1
*
*       Update all relevant COMMON list arrays.
      CALL RENAME
*
*       Check replacing of single KS component by corresponding c.m.
      DO 30 I = IFIRST-2,NTOT-1
   20     IF (LIST(1,I).GT.0.AND.LIST(2,I).LT.IFIRST) THEN
              J2 = LIST(2,I)
              J = KVEC(J2) + N
              IF (LIST(1,I).EQ.1) THEN
                  LIST(2,I) = J
              ELSE
                  L = 2
   22             JNEXT = LIST(L+1,I)
                  IF (JNEXT.LT.J) THEN
                      LIST(L,I) = JNEXT
                      L = L + 1
                      IF (L.LE.LIST(1,I)) GO TO 22
                      LIST(L,I) = J
                  ELSE IF (JNEXT.EQ.J) THEN
                      DO 25 LL = L,LIST(1,I)
                          LIST(LL,I) = LIST(LL+1,I)
   25                 CONTINUE
                      LIST(1,I) = LIST(1,I) - 1
                  ELSE
                      LIST(L,I) = J
                  END IF
*       Check again until first neighbour > ICOMP.
                  GO TO 20
              END IF
          END IF
   30 CONTINUE
*
*       Copy neighbour list for second component & c.m. (NNB1 = LIST(1,I)+1).
      I = 2*NPAIRS - 1
      DO 40 L = 1,NNB1
          LIST(L,I+1) = LIST(L,I)
          LIST(L,NTOT) = LIST(L,I)
   40 CONTINUE
*
*       Initialize the regularized solution.
      CALL KSINIT
      TNEW(NTOT) = T0(NTOT) + STEP(NTOT)
*
*       Check optional binary analysis after merger or multiple collision.
*     IF (KZ(4).GT.0.AND.IPHASE.GT.3) THEN
*         CALL EVOLVE(NPAIRS,-1)
*     END IF
*
*       Check updating of global index for chain c.m.
      IF (NCH.GT.0) THEN
          CALL CHFIND
      END IF
*
      RETURN
*
      END
