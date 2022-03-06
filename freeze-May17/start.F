      SUBROUTINE START
*
*
*       Initialization of data & polynomials.
*       ------------------------------------
*
      INCLUDE 'common6.h'
      EXTERNAL SCALE
      PARAMETER  (NS=12)
*
*
*       Initialize global scalars, counters & useful constants.
      CALL ZERO
*
*       Read input parameters.
      CALL INPUT
*
*       Open all Files.
      CALL FILE_INIT(0)
*
*       Set initial conditions: BODY(I), X(K,I), XDOT(K,I); I=1,N & K=1,3.
      CALL DATA
*
*       Scale initial conditions to new units.
      CALL SCALE
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
*       Introduce scaling factors DAYS, YRS, SU, RAU, SMU, TSTAR & VSTAR.
      CALL UNITS
*
*       Check option for external force.
      IF (KZ(14).GT.0) THEN
          CALL XTRNL0
      END IF
*
*       Check optional scaling to hot system.
      IF (KZ(29).GT.0) THEN
          CALL HOTSYS
      END IF
*
*       Check option for initial binaries.
      IF (KZ(8).EQ.1.OR.KZ(8).GE.3) THEN
          CALL BINPOP
      END IF
*
*       Include stable primordial triples.
      IF (KZ(11).GT.1.AND.KZ(8).GT.0) THEN
          CALL HIPOP
      END IF
*
*       Check optional initialization for tidal two-body capture.
      IF (KZ(27).GT.0) THEN
          CALL INTIDE
      END IF
*
*       Set sequential name, maximum mass & primary velocity.
      BODY1 = 0.0
      DO 20 I = 1,N
          NAME(I) = I
          BODY1 = MAX(BODY1,BODY(I))
          DO 15 K = 1,3
              X0DOT(K,I) = XDOT(K,I)
   15     CONTINUE
   20 CONTINUE
*
*       Initialize fixed block steps (64 levels).
      CALL IBLOCK
*
*       Create table of inverse Stumpff coefficients.
      DO 30 I = 1,NS
          SCOEFF(I) = 1.0/((I + 1)*(I + 2))
   30 CONTINUE
*
*       Set optional stellar evolution parameters.
      IF (KZ(19).GT.2) THEN
          CALL INSTAR
      END IF
*
*       Initialize optional cloud parameters.
      IF (KZ(13).GT.0) THEN
          CALL CLOUD0
      END IF
*
*       Set initial neighbour list & corresponding radius.
      RS0 = RC
      NNB0 = 0
*
      DO 40 I = 1,N
          CALL NBLIST(I,RS0)
          NNB0 = NNB0 + LIST(1,I)
   40 CONTINUE
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
*       Regularize any hard primordial binaries (assume sequential ordering).
      IF (NBIN0.GT.0) THEN
          DO 50 IPAIR = 1,NBIN0
              ICOMP = 2*IPAIR - 1
              JCOMP = 2*IPAIR
              RIJ2 = 0.0
*       Include standard distance criterion.
              DO 45 K = 1,3
                  RIJ2 = RIJ2 + (X(K,ICOMP) - X(K,JCOMP))**2
   45         CONTINUE
              IF (RIJ2.LT.RMIN**2) THEN
                  CALL KSREG
              END IF
   50     CONTINUE
*
      NNBMAX = MIN(N/2,LMAX-3)
*     ZNBMIN = MAX(0.1*FLOAT(NNBMAX),1.0)
      ZNBMIN = MAX(0.01*FLOAT(NNBMAX),1.0)
      ZNBMAX = 0.9*FLOAT(NNBMAX)
*
*       Check initial neighbour lists again.
          DO 55 I = IFIRST,NTOT
              IF(I.GT.N)THEN
              ICOMP = 2*IPAIR - 1
              JCOMP = 2*IPAIR
                  RS0 = RS(ICOMP)
              ELSE
                  RS0 = RC
              END IF
                  CALL NBLIST(I,RS0)
   55     CONTINUE
*
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
*
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
*       Check option for writing the initial conditions on unit 10.
      if(rank.eq.0)then
      IF (KZ(22).EQ.1) THEN
          DO 85 I = 1,N
             WRITE (10,100) BODY(I),(X(K,I),K=1,3),(XDOT(K,I),K=1,3)
   85     CONTINUE
  100 FORMAT(1X,1P,7(1X,D12.5))
      END IF
      end if
*
      RETURN
*
      END
