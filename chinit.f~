      SUBROUTINE CHINIT(ISUB)
*
*
*       Initialization of chain system.
*       -------------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &            NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK,ANG(3),FIRR(3),FD(3)
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
      COMMON/ECHAIN/  ECH
*
*
*       Define chain membership.
      CALL SETSYS
*
*       Initialize c.m. variables.
      DO 2 K = 1,7
          CM(K) = 0.0D0
    2 CONTINUE
*
*       Transform to the local c.m. reference frame.
      DO 4 L = 1,NCH
          J = JLIST(L)
*       Place the system in first single particle locations.
          CM(7) = CM(7) + M(L)
          DO 3 K = 1,3
              X4(K,L) = X(K,J)
              XDOT4(K,L) = XDOT(K,J)
              CM(K) = CM(K) + M(L)*X4(K,L)
              CM(K+3) = CM(K+3) + M(L)*XDOT4(K,L)
    3     CONTINUE
    4 CONTINUE
*
*       Set c.m. coordinates & velocities of subsystem.
      DO 5 K = 1,6
          CM(K) = CM(K)/CM(7)
    5 CONTINUE
*
*       Specify initial conditions for chain regularization.
      LK = 0
      RSUM2 = 0.0
      DO 8 L = 1,NCH
          DO 7 K = 1,3
              LK = LK + 1
              X4(K,L) = X4(K,L) - CM(K)
              XDOT4(K,L) = XDOT4(K,L) - CM(K+3)
              XCH(LK) = X4(K,L)
              VCH(LK) = XDOT4(K,L)
              RSUM2 = RSUM2 + X4(K,L)**2
    7     CONTINUE
    8 CONTINUE
*
*       Calculate internal energy and and save in chain energy.
      CALL CONST(XCH,VCH,M,NCH,ENERGY,ANG,GAM)
      ECH = ENERGY
*
*       Find sum of mass products and define gravitational radius for CHLIST.
      SUM = 0.0D0
      DO 10 L = 1,NCH-1
          DO 9 K = L+1,NCH
              SUM = SUM + M(L)*M(K)
    9     CONTINUE
   10 CONTINUE
      RGRAV = SUM/ABS(ENERGY)
*
*       Introduce characteristic size for initial perturber list.
      RSUM = SQRT(RSUM2)
*
*       Set global index of c.m. body and save name (SUBSYS sets NAME = 0).
      ICH = ICOMP
      NAME0 = NAME(ICH)
*
*       Define subsystem indicator (ISYS = 1, 2, 3 for triple, quad, chain).
      ISYS(NSUB+1) = 3
*
*       Form ghosts and initialize c.m. motion in ICOMP (= JLIST(1)).
      CALL SUBSYS(NCH,CM)
*
*       Copy neighbour list for ghost removal.
      NNB = LIST(1,ICH)
      DO 12 L = 2,NNB+1
          JPERT(L-1) = LIST(L,ICH)
   12 CONTINUE
*
*       Remove ghosts (saved in JLIST) from perturber neighbour lists.
      CALL NBREM(ICH,NCH,NNB)
*
*       Also remove ghosts from list of ICOMP (use NTOT as dummy here).
      JPERT(1) = ICOMP
      CALL NBREM(NTOT,NCH,1)
*
*       Initialize perturber list for integration of chain c.m.
      CALL CHLIST(ICH)
*
*       Perform differential F & FDOT corrections due to perturbers.
      DO 15 K = 1,3
          FIRR(K) = 0.0D0
          FD(K) = 0.0
   15 CONTINUE
      CALL CHFIRR(ICH,0,X(1,ICH),XDOT(1,ICH),FIRR,FD)
      DO 20 K = 1,3
          F(K,ICH) = F(K,ICH) + 0.5*FIRR(K)
          FDOT(K,ICH) = FDOT(K,ICH) + ONE6*FD(K)
          FIDOT(K,ICH) = FIDOT(K,ICH) + FD(K)
          D1(K,ICH) = D1(K,ICH) + FD(K)
   20 CONTINUE
*
*       Take maximum integration interval equal to c.m. step.
      TMAX = STEP(ICOMP)
*
*       Check next treatment time of perturbers & output time.
      CALL TCHAIN(NSUB,TSMIN)
      TMAX = MIN(TMAX,TSMIN)
      TMAX = MIN(TMAX,TADJ - TIME)
*
*       Copy total energy and output & capture option for routine CHAIN.
      CM(8) = BE(3)
      KZ27 = KZ(27)
      KZ30 = KZ(30)
*
*       Assign new subsystem index and begin chain regularization.
      ISUB = NSUB
      NCHAIN = NCHAIN + 1
*
*       Set phase indicator < 0 to ensure new NLIST in routine INTGRT.
      IPHASE = -1
*
      RETURN
*
      END
