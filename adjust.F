      SUBROUTINE ADJUST
*
*
*       Parameter adjustment and energy check.
*       --------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/ECHAIN/  ECH
      SAVE  DTOFF
      DATA  DTOFF /100.0D0/
*
*
*       Predict X & XDOT for all particles (except unperturbed pairs).
      CALL XVPRED(IFIRST,NTOT)
*
*       Obtain the total energy at current time (resolve all KS pairs).
      CALL ENERGY
*
*       Initialize c.m. terms.
      DO 10 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   10 CONTINUE
*
*       Obtain c.m. & angular momentum integrals and Z-moment of inertia.
      AZ = 0.0D0
      ZM = 0.0D0
      ZMASS = 0.0D0
      DO 20 I = 1,N
          ZMASS = ZMASS + BODY(I)
          DO 15 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   15     CONTINUE
          RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                                   (X(3,I) - RDENS(3))**2
          AZ = AZ + BODY(I)*(X(1,I)*XDOT(2,I) - X(2,I)*XDOT(1,I))
          ZM = ZM + BODY(I)*(X(1,I)**2 + X(2,I)**2)
   20 CONTINUE
*
*       Form c.m. coordinates & velocities (vectors & scalars).
      DO 25 K = 1,3
          CMR(K) = CMR(K)/ZMASS
          CMRDOT(K) = CMRDOT(K)/ZMASS
   25 CONTINUE
*
      CMR(4) = SQRT(CMR(1)**2 + CMR(2)**2 + CMR(3)**2)
      CMRDOT(4) = SQRT(CMRDOT(1)**2 + CMRDOT(2)**2 + CMRDOT(3)**2)
*
*       Subtract the kinetic energy of c.m. due to possible cloud effects.
      IF (KZ(13).GT.0) ZKIN = ZKIN - 0.5*ZMASS*CMRDOT(4)**2
*
*       Form virial ratio using single particles & c.m. (isolated or tidal). 
      IF (KZ(14).EQ.0.OR.KZ(14).GT.2) THEN
          Q = ZKIN/POT
          E(3) = ZKIN - POT
      ELSE
*       Use Chandrasekhar eq. (5.535) for virial ratio.
          ZKVIR = ZKIN + 0.5*TIDAL(4)*AZ
          Q = ZKVIR/(POT - 2.0*ETIDE)
*       Modify angular momentum integral using Chandrasekhar eq. (5.530).
          AZ = AZ + 0.5*TIDAL(4)*ZM
          E(3) = ZKIN - POT + ETIDE
      END IF
*
*       Define crossing time using single particle energy (cf. option 14).
      TCR = ZMASS**2.5/(2.0*ABS(E(3)))**1.5
      IF (Q.GT.1.0) THEN
          TCR = TCR*SQRT(2.0*Q)
      END IF
*       Form provisional total energy.
      ETOT = ZKIN - POT + ETIDE
*
*       Include KS pairs, triple, quad, mergers, collisions & chain.
      ETOT = ETOT + EBIN + ESUB + EMERGE + ECOLL + EMDOT + ECDOT
      IF (NCH.GT.0) THEN
          ETOT = ETOT + ECH
      END IF
*
*       Update energies and form the relative error (divide by ZKIN or E(3)).
      IF (TIME.LE.0.0D0) THEN
          DE = 0.0D0
          BE(1) = ETOT
          BE(3) = ETOT
          DELTA1 = 0.0D0
      ELSE
          BE(2) = BE(3)
          BE(3) = ETOT
          DE = BE(3) - BE(2)
          DELTA1 = DE
          DETOT = DETOT + DE
          DE = DE/MAX(ZKIN,ABS(E(3)))
*       Save sum of relative energy error for main output and accumulate DE.
          ERROR = ERROR + DE
          ERRTOT = ERRTOT + DE
      END IF
*
*       Set provisional half-mass radius.
      RSCALE = 0.5*ZMASS**2/POT
*
*       Determine average neighbour number and smallest neighbour sphere.
      NNB = 0
      RS0 = RSCALE
      DO 40 I = IFIRST,NTOT
          NNB = NNB + LIST(1,I)
          IF (LIST(1,I).GT.0) RS0 = MIN(RS0,RS(I))
   40 CONTINUE
      NNB = FLOAT(NNB)/FLOAT(N - NPAIRS)
*
*       Use current value if minimum neighbour sphere not implemented.
      IF (RSMIN.EQ.0.0D0) RSMIN = RS0
*
*       Find density centre & core radius (Casertano & Hut, Ap.J. 298, 80).
      IF (N-NPAIRS.GE.20.AND.KZ(29).EQ.0) THEN
          CALL CORE
      ELSE
          NC = N
          ZMC = ZMASS
          RHOD = 1.0
          RHOM = 1.0
          RC = RSCALE
          RC2 = RC**2
          RC2IN = 1.0/RC2
      END IF
*
*       Check optional sorting of Lagrangian radii & half-mass radius.
      IF (KZ(7).GT.0) THEN
          CALL LAGR(RDENS)
      END IF
*
*       Scale average & maximum core density by the mean value.
      RHOD = 4.0*TWOPI*RHOD*RSCALE**3/(3.0*ZMASS)
      RHOM = 4.0*TWOPI*RHOM*RSCALE**3/(3.0*ZMASS)
*
*       Adopt density contrasts of unity for hot system.
      IF (KZ(29).GT.0.AND.ZKIN.GT.POT) THEN
          RHOD = 1.0
          RHOM = 1.0
      END IF
*
*       Check optional determination of regularization parameters.
      IF (KZ(16).GT.0) THEN
          RMIN0 = RMIN
*
*       Form close encounter distance from scale factor & density contrast.
          RMIN = 4.0*RSCALE/(FLOAT(N)*RHOD**0.3333)
*       Use harmonic mean to reduce fluctuations (avoid initial value).
          IF (TIME.GT.0.0D0) RMIN = SQRT(RMIN0*RMIN)
*       Impose maximum value for sufficient perturbers.
          RMIN = MIN(RMIN,RSMIN*GMIN**0.3333)
*       Define scaled DTMIN by RMIN & <M> and include ETAI for consistency.
          DTMIN = 0.04*SQRT(ETAI/0.02D0)*SQRT(RMIN**3/BODYM)
*       Specify binding energy per unit mass of hard binary (impose Q = 0.5).
          ECLOSE = 4.0*MAX(ZKIN,ABS(ZKIN - POT))/ZMASS
*       Adopt central velocity as upper limit (avoids large kick velocities).
          IF (2.0*ZKIN/ZMASS.GT.VC**2) ECLOSE = 2.0*VC**2
          IF (Q.GT.0.5) THEN
              ECLOSE = 0.5*ECLOSE/Q
          END IF
          ECLOSE = MIN(ECLOSE,1.0D0)
      END IF
*
*       Check optional modification of DTMIN, ECLOSE & TCR for hot system.
      IF (KZ(29).GT.0.AND.Q.GT.1.0) THEN
          DTMIN = 0.04*SQRT(ETAI/0.02D0)*SQRT(RMIN**3/BODYM)
          SIGMA2 = 2.0*ZKIN/ZMASS
          VP2 = 4.0*BODYM/RMIN
          DTMIN = DTMIN*SQRT((VP2 + SIGMA2/Q)/(VP2 + 2.0D0*SIGMA2))
          ECLOSE = SIGMA2
          TCR = 2.0*RSCALE/SQRT(SIGMA2)
      END IF
*
*       Set useful scalars for the integrator.
      SMIN = 2.0*DTMIN
      RMIN2 = RMIN**2
      RMIN22 = 4.0*RMIN2
      EBH = -0.5*BODYM*ECLOSE
*       Specify square close encounter force from heaviest body at 5*RMIN.
      FCRIT2 = (BODY1/(25.0*RMIN2))**2
*       Adopt 2*RSMIN for neighbour sphere volume factor in routine REGINT.
      RSFAC = MAX(25.0/TCR,3.0D0*VC/(2.0D0*RSMIN))
*
*       Update density contrast factor for neighbour sphere modification.
      IF (TIME.LE.0.0D0.OR.KZ(40).EQ.0) THEN
          ALPHA = FLOAT(NNBMAX)*SQRT(0.08D0*RSCALE**3/FLOAT(N-NPAIRS))
      END IF
*       Include optional stabilization to increase neighbour number.
      IF (KZ(40).GT.0.AND.FLOAT(NNB).LT.0.5*NNBMAX) THEN
          FAC = 1.0 + (0.5*NNBMAX - NNB)/(0.5*FLOAT(NNBMAX))
          ALPHA = FAC*ALPHA
      END IF
*
*       Define tidal radius for isolated system (2*RTIDE used in ESCAPE).
      IF (KZ(14).EQ.0) RTIDE = 10.0*RSCALE
*
*       Print energy diagnostics & KS parameters.
      ICR = TTOT/TCR
      WRITE (6,50)  TTOT, Q, DE, BE(3),RMIN, DTMIN, ECLOSE, ICR, DELTA1
   50 FORMAT (/,' ADJUST:  TIME =',F8.2,'  Q =',F5.2,'  DE =',1P,E10.2,
     &          '  E =',0P,F10.6,'  RMIN =',1P,E8.1,'  DTMIN =',E8.1,
     &          '  ECLOSE =',0PF5.2,'  TC =',I5,'  DELTA =',1P,E9.1)
      CALL FLUSH(6)
*
*       Perform automatic error control (RETURN on restart with KZ(2) > 1).
      CALL CHECK(DE)
      IF (ABS(DE).GT.5.0*QE) GO TO 70
*
*       Check for escaper removal.
      IF (KZ(23).GT.0) THEN
          CALL ESCAPE
      END IF
*
*       Check correction for c.m. displacements.
      IF (KZ(31).GT.0) THEN
          CALL CMCORR
      END IF
*
*       See whether standard output is due (allow for setting TIME = TPREV).
      IF (TIME.GE.TNEXT - 20.0*DTMIN) THEN
          CALL OUTPUT
      END IF
*
*       Update time for next adjustment.
      TADJ = TADJ + DTADJ
*
*       Obtain elapsed CPU time and update total since last output/restart.
      CALL CPUTIM(TCOMP)
      CPUTOT = CPUTOT + TCOMP - CPU0
      CPU0 = TCOMP
*
*       Save COMMON after energy check (skip TRIPLE, QUAD, CHAIN).
      TDUMP = TIME
      IF (KZ(2).GE.1.AND.NSUB.EQ.0) CALL MYDUMP(1,2)
*
*       Check termination criteria (TIME > TCRIT & N <= NCRIT).
      IF (TTOT.GT.TCRIT - 20.0*DTMIN.OR.N.LE.NCRIT) THEN
*       Terminate after optional COMMON save.
          WRITE (6,60)  TTOT, CPUTOT/60.0, ERRTOT, DETOT
   60     FORMAT (//,9X,'END RUN',3X,'TIME =',F7.1,'  CPUTOT =',F7.1,
     &                  '  ERRTOT =',F10.6,'  DETOT =',F10.6)
          IF (KZ(1).GT.0.AND.NSUB.EQ.0) CALL MYDUMP(1,1)
          STOP
      END IF
*
*       Check optional truncation of time.
      IF (KZ(35).GT.0.AND.TIME.GE.DTOFF) THEN
          CALL OFFSET(DTOFF)
      END IF
*
   70 RETURN
*
      END
