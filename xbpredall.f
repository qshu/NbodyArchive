      subroutine xbpredall
*
*
*     Predict x and xdot. (L.WANG)

      INCLUDE 'common6.h'
      INCLUDE 'omp_lib.h'
      COMMON/XPRED/ TPRED(NMAX),predall
      LOGICAL PREDALL
      REAL*8 TPRED

      IF (PREDALL) RETURN

      NNPRED = NNPRED + 1
!     $omp parallel do default(shared)
!     $omp& private(J,S,S1,S2)
      DO 40 J = IFIRST,NTOT
*     IF(TPRED(J).NE.TIME) THEN
         S = TIME - T0(J)
         S1 = 1.5*S
         S2 = 2.0*S
         X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S +X0(1,J)
         X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S +X0(2,J)
         X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S +X0(3,J)
         XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
         XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
         XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
         TPRED(J) = TIME
*     END IF
 40   CONTINUE
!     $omp end parallel do
      
      IF (NPAIRS.GT.0) THEN
*     Resolve perturbed KS pairs with c.m. prediction after NBSORT.
*     JJ = -1
*     *!$omp parallel do IF(NPAIRS.GE.ITHREAD) default(shared)
*     *!$omp& private(J,JPAIR,S,J1,J2,A1,A2,A3,A4,IMOD,DTU,DTU1,DTU2,
*     *!$omp&        Q1,Q2,Q3,UI,RDOT,V,RI,RINV)
*     *!$omp& reduction(+:NBPRED)
         DO 45 JPAIR = 1,NPAIRS
*     IF (.NOT.NOPRED(JJ).OR..NOT.NOPRED(JJ+1)) THEN
            IF (LIST(1,2*JPAIR - 1).GT.0) THEN
*     Ignore c.m. prediction after full N loop (all active KS needed).
               J = N + JPAIR
*     Predict ALL binaries even unperturbed ones for parallel code (R.Sp.)
               ZZ = 1.0
*     Distinguish between low and high-order prediction of U & UDOT.
               IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
               CALL KSRES2(JPAIR,J1,J2,ZZ,TIME)
            END IF
 45      CONTINUE
      end if

      PREDALL=.true.

      return

      end
      
