      PROGRAM KSTEST
      INCLUDE 'common6.h'
      INCLUDE 'kscmn6.h'

      CALL ZERO
      CALL DEFINE
      SMAX = 0.125
      CALL IBLOCK

      N = 3

      NNBMAX = N
      RS0    = 10.0
      DTMIN  = 0.001
      RMIN   = 1.0
      ETAI   = 0.02
      ETAR   = 0.02
      ETAU   = 0.1
      ECLOSE = 1.0
      GMIN   = 1d-6
      GMAX   = 1d-3
      CLIGHT = 0.0
      TCRIT  = 100000.0
      QE     = 1e-2

      NTOT = N
!     NZERO = N
!     NBZERO = NNBMAX
!     ZNBMIN = 0.2*FLOAT(NNBMAX)
!     ZNBMAX = 0.9*FLOAT(NNBMAX)
      CMSEP2 = GMIN**(-0.666666667)
      RMIN2 = RMIN**2
      RSMIN = RS0
      RC = RS0

      KZ(10) = 0 !3

      SEMI = 0.1*RMIN
      ECC = 0.9999999

      BODY(1) = 0.5
      BODY(2) = 1.0 - BODY(1)
      BODY(3) = 1e-7
      ZMASS = 1.0
      BODYM = ZMASS/FLOAT(N)

      RB = SEMI*(1.0+ECC)
      VB = SQRT(2.0*(1.0-ECC)/RB)

      DO I = 1,N
          DO K = 1,3
              X(K,I) = 0.0
              XDOT(K,I) = 0.0
          END DO
      END DO
      X(1,1)    =  BODY(2)*RB
      XDOT(2,1) =  BODY(2)*VB
      X(1,2)    = -BODY(1)*RB
      XDOT(2,2) = -BODY(1)*VB
      X(2,3)    =  0.5*RS0
!     XDOT(1,3) =  1.0/SQRT(X(2,3))

      BODY1 = 0.0
      DO 20 I = 1,N
          NAME(I) = I
          BODY1 = MAX(BODY1,BODY(I))
          DO 15 K = 1,3
              X0DOT(K,I) = XDOT(K,I)
   15     CONTINUE
   20 CONTINUE

      DO 40 I = 1,N
          RS(I) = RS0
          LIST(1,I) = N-1
          K = 2
          DO J = 1,N
              IF (J.NE.I) THEN
                  LIST(K,I) = J
                  K = K+1
              END IF
          END DO
   40 CONTINUE

      CALL FPOLY1(1,N,0)
      CALL FPOLY2(1,N,0)

      ICOMP = 1
      JCOMP = 2
      IPHASE = 1
      CALL KSREG
      IPHASE = 0

      D0(1:3,1:NTOT) = 0.0
      D1(1:3,1:NTOT) = 0.0
      D2(1:3,1:NTOT) = 0.0
      D3(1:3,1:NTOT) = 0.0
      D0R(1:3,1:NTOT) = 0.0
      D1R(1:3,1:NTOT) = 0.0
      D2R(1:3,1:NTOT) = 0.0
      D3R(1:3,1:NTOT) = 0.0
      FI(1:3,1:NTOT) = 0.0
      FIDOT(1:3,1:NTOT) = 0.0
      FR(1:3,1:NTOT) = 0.0
      FRDOT(1:3,1:NTOT) = 0.0
      F(1:3,1:NTOT) = 0.0
      FDOT(1:3,1:NTOT) = 0.0

      TIME_CH = TTOT

      CALL ENERGY
      ETOT = ZKIN - POT + ETIDE + EPL
      ETOT = ETOT + EBIN + ESUB + EMERGE + ECOLL + EMDOT + ECDOT
      E0 = ETOT

      EPREV = 0.0
1     CALL XVPRED(IFIRST,NTOT)
      CALL ENERGY
      ETOT = ZKIN - POT + ETIDE + EPL
      ETOT = ETOT + EBIN + ESUB + EMERGE + ECOLL + EMDOT + ECDOT
!     write(6,*) ' '
!     write(6,*) 'T=',TIME,'E=',ETOT,'DT=',STEP(1),'DE=',ETOT-E0
!     write(6,*) 'X =',X(1:2,1:2)
!     write(6,*) 'XD=',XDOT(1:2,1:2)
      EPREV=ETOT

      IF (TIME.GT.TCRIT) THEN
          write(6,*) 'dE=',ETOT-E0
          STOP
      END IF

      TIME = T0(1) + STEP(1)
      ISTAT(1) = 0
      I1 = 1
      CALL KSINT(I1,1)
      TTOT = TIME

      GO TO 1

      END
