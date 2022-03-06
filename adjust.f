      SUBROUTINE ADJUST
*
*
*       Parameter adjustment and energy check.
*       --------------------------------------
*
      INCLUDE 'common4.h'
      INTEGER I,K,JPAIR,ICR
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
*       Form virial ratio using single particles & c.m. (isolated or tidal). 
      IF (KZ(14).EQ.0) THEN
          Q = ZKIN/POT
      ELSE 
*       Form the generalized virial ratio (cf. Fukushige & Heggie, MN 318).
          IF (KZ(14).EQ.1) THEN
*       Use Chandrasekhar eq. (5.535) for virial ratio (rotating frame only).
              Q = ZKIN/(POT - 2.0*(ETIDE + 0.5*TIDAL(4)*AZ))
*       Modify angular momentum integral by Chandrasekhar eq. (5.530).
              AZ = AZ + 0.5*TIDAL(4)*ZM
          ELSE
              Q = ZKIN/(POT - 2.0*ETIDE)
          END IF
      END IF
*
*       Define crossing time and save single particle energy.
      ETOT = ZKIN - POT + ETIDE + EBLCKHL
*       Add Coriolis term to Jacobi energy integral (exclude standard case).
      IF (KZ(14).GE.2) ETOT = ETOT - 0.5*TIDAL(4)*AZ
      TCR = ZMASS**2.5/(2.0*ABS(ETOT))**1.5
      IF (Q.GT.1.0) THEN
          TCR = TCR*SQRT(2.0*Q)
      END IF
      E(3) = ETOT
*
*       Include KS, triple & quad, mergers, mass loss, collisions & chain.
      ETOT = ETOT + EBIN + ESUB + EMERGE + EMDOT + ECOLL + ECDOT
      IF (NCH.GT.0) THEN
          ETOT = ETOT + ECH
      END IF
*
*       Update energies and form the relative error (divide by ZKIN or E(3)).
      IF (TIME.LE.0.0D0) THEN
          DE = 0.0D0
          BE(1) = ETOT
          BE(3) = ETOT
      ELSE
          BE(2) = BE(3)
          BE(3) = ETOT
          DE = BE(3) - BE(2)
          DETOT = DETOT + DE
*         DE = DE/MAX(ZKIN,ABS(ETOT))
          DE = DE/MAX(ZKIN,ABS(E(3)))
*       Save sum of relative energy error for main output and accumulate DE.
          ERROR = ERROR + DE
          ERRTOT = ERRTOT + DE
      END IF
*
*       Set provisional half-mass radius.
          RSCALE = 0.5*ZMASS**2/POT
*
*       Find density centre & core radius (Casertano & Hut, Ap.J. 298, 80).
      IF (N.GT.50.AND.KZ(29).EQ.0) THEN
*       Send all single particles & c.m. to GRAPE for neighbour list.
          CALL GPSEND
*         IF (TIME.GE.TPRINT - 20.0*DTMIN) THEN
*             CALL GPSEND
              CALL CORE
*         ELSE
*             RHOD = 1.0
*             RHOM = 1.0
*         END IF
      ELSE
*       Adopt density centre at zero for small N.
          NC = N
          ZMC = ZMASS
          RC = RSCALE
          RHOD = 1.0
          RHOM = 1.0
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
*       Introduce equilibrium half-mass radius (predicted or actual).
          IF (KZ(7).EQ.0.AND.KZ(27).LE.0) THEN
*       Include merger energy to avoid small denominator (avoid dissipation).
              REQ = 0.25*ZMASS**2/ABS(ZKIN - POT + EMERGE)
          ELSE
              REQ = RSCALE
          END IF
*
*       Form close encounter distance from scale factor & density contrast.
*         RMIN = 4.0*RSCALE/(FLOAT(N)*RHOM**0.3333)
          RMIN = 4.0*RSCALE/(FLOAT(N)*RHOD**0.3333)
*         RMIN = 4.0*RSCALE/(FLOAT(N)*MAX(RHOM,RHOD)**0.3333)
*         RMIN = MAX(RMIN,5.0D-05)
*       Use harmonic mean to reduce fluctuations (avoid initial value).
          IF (TIME.GT.0.0D0) RMIN = SQRT(RMIN0*RMIN)
*       Impose maximum value for increased efficiency (small N limit).
          RMIN = MIN(RMIN,0.01*REQ)
*       Define scaled DTMIN by RMIN & <M> and include ETA for consistency.
          DTMIN = 0.04*SQRT(ETA/0.02D0)*SQRT(RMIN**3/BODYM)
*       Specify binding energy per unit mass of hard binary (impose Q = 0.5).
          ECLOSE = 4.0*MAX(ZKIN,ABS(ZKIN - POT))/ZMASS
*       Adopt central velocity as upper limit (avoids large kick velocities).
          IF (2.0*ZKIN/ZMASS.GT.VC**2) ECLOSE = 2.0*VC**2
          IF (Q.GT.0.5) THEN
              ECLOSE = 0.5*ECLOSE/Q
          END IF
      END IF
*
*       Check optional modification of DTMIN, ECLOSE & TCR for hot system.
      IF (KZ(29).GT.0.AND.Q.GT.1.0) THEN
          DTMIN = 0.04*SQRT(ETA/0.02D0)*SQRT(RMIN**3/BODYM)
          SIGMA2 = 2.0*ZKIN/ZMASS
          VP2 = 4.0*BODYM/RMIN
          DTMIN = DTMIN*SQRT((VP2 + SIGMA2/Q)/(VP2 + 2.0D0*SIGMA2))
          ECLOSE = SIGMA2
          TCR = 2.0*RSCALE/SQRT(SIGMA2)
      END IF
*
*       Set useful scalars for the integrator and hard binary energy.
      SMIN = 2.0*DTMIN
      RMIN2 = RMIN**2
      RMIN22 = 4.0*RMIN2
      EBH = -0.5*BODYM*ECLOSE
*     IF (NZERO.LE.1000) EBH = -0.25*(BODY1 + BODYM)*ECLOSE
      IF (TIME.LE.0.0D0) STEPJ = 0.01*(60000.0/FLOAT(N))**1.5
*
*       Determine mass factor (2/3 power) for possible use.
      ZMB = BODYM
      DO 30 JPAIR = 1,NPAIRS
          IF (BODY(N+JPAIR).GT.0) THEN
              ZMB = MIN(BODY(N+JPAIR),ZMB)
          END IF
   30 CONTINUE
      BODY23 = (2.0*BODY1/ZMB)**0.6667
*
*       Define tidal radius for isolated system (2*RTIDE used in ESCAPE).
      IF (KZ(14).NE.-1.AND.TIDAL(1).EQ.0.0) RTIDE = 10.D0*RSCALE
*
*       Print energy diagnostics & KS parameters.
      ICR = TTOT/TCR
      WRITE (6,50)  TTOT, Q, DE, BE(3), RMIN, DTMIN, ECLOSE, ICR
   50 FORMAT (/,' ADJUST:  TIME =',F8.2,'  Q =',F5.2,'  DE =',F10.6,
     &          '  E =',F10.6,'  RMIN =',1PE8.1,'  DTMIN =',E8.1,
     &          '  ECLOSE =',0PF5.2,'  TC =',I4)
      CALL FLUSH(3)
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
*     WRITE(79,*)TTOT,N,ZMASS*ZMBAR,RSCALE*RBAR
*
*       See whether standard output is due (allow for setting TIME = TPREV).
      IF (TIME.GE.TPRINT - 20.0*DTMIN) THEN
          WRITE (40,55)  TTOT, NC, RC, ZMC, NHIVEL, NPAIRS, EBIN
   55     FORMAT (' T NC RC ZMC NHIVEL KS EB ',
     &              F8.1,I5,2F7.3,I11,I4,F10.6)
          CALL FLUSH(40)
          CALL OUTPUT
*
*       Reset GRAPE and send all particles (temporary measure). 
*
*         CALL GPWIPE(GPID,TIME)
*         CALL GPSEND
      END IF
*
*       Check optional truncation of time.
      IF (KZ(35).GT.0.AND.TIME.GE.DTOFF) THEN
          CALL OFFSET(DTOFF)
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
*       Save COMMON after satisfactory energy check.
      TDUMP = TIME
      IF (KZ(2).GE.1.AND.NSUB.EQ.0) CALL MYDUMP(1,2)
*
*       Check termination criteria (TPHYS > TCRIT, TIME > TCRIT, N <= NCRIT).
      IF (KZ(14).GT.1.AND.TPHYS.GT.TCRIT) THEN
          TCRIT = TTOT
      END IF
      IF (TTOT.GT.TCRIT - 20.0*DTMIN.OR.N.LE.NCRIT) THEN
*       Terminate after optional COMMON save.
          WRITE (6,60)  TTOT, CPUTOT/60.0, ERRTOT, DETOT
   60     FORMAT (//,9X,'END RUN',3X,'TIME =',F7.1,'  CPUTOT =',F7.1,
     &                  '  ERRTOT =',F10.6,'  DETOT =',F10.6)
          IF (KZ(1).GT.0.AND.NSUB.EQ.0) CALL MYDUMP(1,1)
          CALL gpfree
          STOP
      END IF
*
   70 RETURN
*
      END
