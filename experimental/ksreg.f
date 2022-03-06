      SUBROUTINE KSREG
*
*
*       New KS regularization.
*       ----------------------
*
      INCLUDE 'common6.h'
      REAL*8  SAVE(13)
      INTEGER ILINK1(NMAX),ILINK2(NMAX)
      EXTERNAL RENAME
      COMMON/LINKLIST/LINK(NMAX,NMAX)
*
*
*       Replace #JCOMP by arbitrary body in case it is the only neighbour.
      NNB = LIST(1,ICOMP)
      IF (NNB.EQ.1.AND.LIST(2,ICOMP).EQ.JCOMP) THEN
          LIST(2,ICOMP) = JCOMP + 1
          IF (JCOMP + 1.GT.NTOT) THEN
              LIST(2,ICOMP) = MAX(ICOMP-1,IFIRST+2)
          END IF
      END IF
*
*       Copy neighbour list of #ICOMP without JCOMP.
      NNB1 = 1
      DO 11 L = 1,NNB
          IF (LIST(L+1,ICOMP).EQ.JCOMP) GO TO 11
          NNB1 = NNB1 + 1
          ILIST(NNB1) = LIST(L+1,ICOMP)
   11 CONTINUE
      ILIST(1) = NNB1 - 1
*
*       Copy neighbour link list of #ICOMP without JCOMP.
      NNBX = LINK(1,ICOMP)
      NNBX1 = 1
      DO 12 L = 1,NNBX
          IF (LINK(L+1,ICOMP).EQ.JCOMP) GO TO 12
          NNBX1 = NNBX1 + 1
          ILINK1(NNB1) = LINK(L+1,ICOMP)
   12 CONTINUE
      ILINK1(1) = NNBX1 - 1
*
*       Copy neighbour link list of #JCOMP without ICOMP.
      NNBX = LINK(1,JCOMP)
      NNBX2 = 1
      DO 1 L = 1,NNBX
          IF (LINK(L+1,JCOMP).EQ.ICOMP) GO TO 1
          NNBX2 = NNBX2 + 1
          ILINK2(NNB1) = LINK(L+1,JCOMP)
    1 CONTINUE
      ILINK2(1) = NNBX2 - 1
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
          K = LINK(1,J) + 1
          DO 6 L =1,K
              LINK(L,I) = LINK(L,J)
    6     CONTINUE
          BODY(J) = SAVE(7)
          RS(J) = SAVE(8)
          RADIUS(J) = SAVE(9)
          TEV(J) = SAVE(10)
          BODY0(J) = SAVE(11)
          TEV0(J) = SAVE(12)
          EPOCH(J) = SAVE(13)
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
      PRINT*,' KSREG for ICOMP,JCOMP,IPAIR=',ICOMP,JCOMP,NPAIRS
*       Check replacing of single KS component by corresponding c.m.
      I = 2*NPAIRS - 1
   20 IF (LIST(2,I).LT.IFIRST) THEN
          PRINT*,' Found A Particle ',I,' having NB ',LIST(2,I)
          J = 0.5001*(LIST(2,I) + 1) + N
          DO 30 L = 2,NNB1
              IF (L.LT.NNB1.AND.LIST(L+1,I).LT.J) THEN
                  LIST(L,I) = LIST(L+1,I)
              ELSE
                  LIST(L,I) = J
                  GO TO 20
              END IF
   30     CONTINUE
*       Check again until first neighbour > 2*NPAIRS.
          GO TO 20
      END IF
*
      IF(TIME.GT.0.D0)THEN
      NNBX1 = ILINK1(1) + 1
      DO 455 K = 2,NNBX1
          I = ILINK1(K)
 255  IF (LIST(2,I).LT.IFIRST) THEN
          PRINT*,' Found B Particle ',I,' having NB ',LIST(2,I)
          J = 0.5001*(LIST(2,I) + 1) + N
          DO 355 L = 2,NNB1
              IF (L.LT.NNB1.AND.LIST(L+1,I).LT.J) THEN
                  LIST(L,I) = LIST(L+1,I)
              ELSE
                  LIST(L,I) = J
                  GO TO 255
              END IF
  355     CONTINUE
*       Check again until first neighbour > 2*NPAIRS.
          GO TO 255
      END IF
  455 CONTINUE
      END IF
*      
*     PRINT*,' ICOMP NNB Bef=',ICOMP,LIST(1,ICOMP),
*    *                     (LIST(K,ICOMP),K=2,LIST(1,ICOMP))
*     PRINT*,' JCOMP NNB Bef=',JCOMP,LIST(1,JCOMP),
*    *                     (LIST(K,JCOMP),K=2,LIST(1,JCOMP))
*                 RS0 = RS(ICOMP)
*                 CALL NBLIST(NTOT,RS0)
*     PRINT*,' IPAIR NNB Aft=',IPAIR,LIST(1,NTOT),
*    *                     (LIST(K,NTOT),K=2,LIST(1,NTOT))
*       Copy neighbour list for second component & c.m.
*     DO 40 L = 1,NNB1
*         LIST(L,ICOMP) = LIST(L,NTOT)
*         LIST(L,JCOMP) = LIST(L,NTOT)
*  40 CONTINUE
*
      I = 2*NPAIRS - 1
      DO 40 L = 1,NNB1
          LIST(L,I+1) = LIST(L,I)
          LIST(L,NTOT) = LIST(L,I)
   40 CONTINUE
*       Initialize the regularized solution.
      CALL KSINIT
*
*       Check optional binary analysis after merger or multiple collision.
      IF (KZ(4).GT.0.AND.IPHASE.GT.3) THEN
          CALL EVOLVE(NPAIRS,-1)
      END IF
*
*       Check updating of global index for chain c.m.
      IF (NCH.GT.0) THEN
          CALL CHFIND
      END IF
*
      RETURN
*
      END
