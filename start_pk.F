      SUBROUTINE START
*
*
*       Initialization of data & polynomials.
*       ------------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Initialize global scalars, counters & useful constants.
      CALL ZERO
*
*       Read input parameters.
      CALL INPUT
*
*       Set initial conditions: BODY(I), X(K,I), XDOT(K,I); I=1,N & K=1,3.
      CALL DATA
*
*       Scale initial conditions to new units.
      CALL SCALEA
*
*       Set total mass in case routines DATA & SCALE are not used.
      ZMASS = 0.0D0
      DO 10 I = 1,N
          ZMASS = ZMASS + BODY(I)
   10 CONTINUE
*
*       Define mean mass in scaled units and solar mass conversion factor.
      BODYM = ZMASS/FLOAT(N)
      ZMBAR = ZMBAR/BODYM
*
*       Check option for external force.
      IF (KZ(14).GT.0) THEN
          CALL XTRNL0
      ELSE
C New (Aug. 1998, P.Kroupa, following nb5):
*       Define nominal time scale & velocity scale (RBAR = 1 & ZMBAR = 1).
          TSCALE_pk = 14.94/SQRT(FLOAT(N))
          VSTAR_pk = 0.0655*SQRT(FLOAT(N))
C Old:
          TSCALE = 2.13/TCR
          DO 20 K = 1,4
              TIDAL(K) = 0.0D0
   20     CONTINUE
      END IF
*
*       Check optional scaling to hot system.
      IF (KZ(29).GT.0) THEN
          CALL HOTSYS
      END IF
*
*       Check option for initial binaries.
      IF (KZ(8).EQ.1.OR.KZ(8).EQ.3) THEN
          CALL BINPOP
      END IF
*
*       Check optional initialization for tidal two-body capture.
      IF (KZ(27).NE.0) THEN
          CALL INTIDE
      END IF
*
*       Set sequential name, maximum mass & primary velocity.
      BODY1 = 0.0
      DO 30 I = 1,N
          NAME(I) = I
          BODY1 = MAX(BODY1,BODY(I))
          DO 25 K = 1,3
              X0DOT(K,I) = XDOT(K,I)
   25     CONTINUE
   30 CONTINUE
*
*       Set optional stellar evolution parameters.
      IF (KZ(19).GT.2) THEN
          CALL INSTAR
      END IF
*
*       Initialize optional cloud parameters.
      IF (KZ(13).NE.0) THEN
          CALL CLOUD0
      END IF
*
*       Set initial neighbour list & corresponding radius.
      RS0 = RC
      NNB0 = 0
      DO 40 I = 1,N
          CALL NBLIST(I,RS0)
          NNB0 = NNB0 + LIST(1,I)
   40 CONTINUE
*
*       Initialize fixed block steps (32 levels).
      CALL IBLOCK
*
*       Obtain force & first derivative.
      call cputim(tt1)
#ifdef PARALLEL
      CALL FPOLY1_MPI(1,N,0)
#else
      CALL FPOLY1(1,N,0)
#endif
      call cputim(tt2)
      if(rank.eq.0)print*,' fpoly1 time=',(tt2-tt1)*60.
*
*       Obtain second & third force derivatives and set time-steps.
      call cputim(tt1)
#ifdef PARALLEL
      CALL FPOLY2_MPI(1,N,0)
#else
      CALL FPOLY2(1,N,0)
#endif
      call cputim(tt2)
      if(rank.eq.0)print*,' fpoly2 time=',(tt2-tt1)*60.
*
*       Regularize any primordial binaries (assume sequential ordering).
      IF (NBIN0.GT.0) THEN
          DO 50 IPAIR = 1,NBIN0
              ICOMP = 2*IPAIR - 1
              JCOMP = 2*IPAIR
C New (Aug. 1998, P.Kroupa, following nb5):
              RIJ2 = 0.0
*       Include standard distance criterion.
              DO 45 K = 1,3
                  RIJ2 = RIJ2 + (X(K,ICOMP) - X(K,JCOMP))**2
   45         CONTINUE
              IF (RIJ2.LT.RMIN**2) THEN
                  CALL KSREG
              END IF
C Old:
c             CALL KSREG
C------
   50     CONTINUE
      END IF
*
*       Initialize the time-step list used to find next body.
      DTLIST = 100.0
      DO 60 I = 1,NTOT
          DTLIST = MIN(DTLIST,STEP(I))
   60 CONTINUE
*
*       Set initial time-step list interval twice the smallest step.
      DTLIST = 2.0*DTLIST
   70 NNB = 1
      TLIST = TLIST + DTLIST
*
*       Select all members due in the interval (0,TLIST).
      DO 80 J = 1,NTOT
          IF (T0(J) + STEP(J).LT.TLIST) THEN
              NNB = NNB + 1
              NLIST(NNB) = J
          END IF
   80 CONTINUE
*
*       Check whether membership range is acceptable.
      IF (NNB.EQ.1) GO TO 70
      IF (NNB.GT.LMAX) THEN
          TLIST = TLIST - DTLIST
          DTLIST = 0.66*DTLIST
          GO TO 70
      END IF
*
*       Reduce new DTLIST to prevent early crowding and set membership.
      DTLIST = 0.2*DTLIST
      NLIST(1) = NNB - 1
*
*       Check the average neighbour number.
      ZNB = FLOAT(NNB0)/FLOAT(N)
      IF (ZNB.LT.0.25*ZNBMAX) THEN
          if(rank.eq.0)WRITE (6,90)  ZNB
   90     FORMAT (/,12X,'WARNING!   SMALL NEIGHBOUR NUMBERS   <NNB> =',
     &                                                             F5.1)
      END IF
*
      RETURN
*
      END
