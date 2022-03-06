      SUBROUTINE OUTPUT
*
*
*       Output and data save.
*       ---------------------
*
      INCLUDE 'common6.h'
      COMMON/ECHAIN/  ECH
      REAL*4  XS(3,NMAX),VS(3,NMAX),BODYS(NMAX),RHOS(NMAX),AS(20)
      REAL*4  XNS(NMAX),PHI(NMAX)
*
*         Call Computation of Moments of Inertia (with Chr. Theis)
          IF(KZ(39).EQ.2)CALL ELLAN
*       Obtain energy error in case routine ADJUST not called recently.
      IF (TIME.GE.TADJ.OR.TIME.LE.0.0D0) GO TO 10
*
*       Predict X & XDOT for all particles (except unperturbed pairs).
      CALL XVPRED(IFIRST,NTOT)
*
*       Obtain the total energy at current time (resolve all KS pairs).
#ifdef PARALLEL
      CALL ENERGY_MPI
#else
      CALL ENERGY
#endif
*
*       Include KS pairs, triple, quad, mergers, collisions & chain.
      ETOT = ZKIN - POT + ETIDE + EBIN + ESUB + EMERGE + ECOLL + EMDOT +
     &       ECH
*
*       Update energies and form the relative error (divide by ZKIN or ETOT).
      BE(2) = BE(3)
      BE(3) = ETOT
      DE = BE(3) - BE(2)
      DETOT = DETOT + DE
      DE = DE/MAX(ZKIN,ABS(ETOT))
*       Save sum of relative energy error for main output and accumulate DE.
      ERROR = ERROR + DE
*
*       Find density centre & core radius (Casertano & Hut, Ap.J. 298, 80).
      IF (N.GE.20.AND.KZ(29).EQ.0) THEN
          CALL CORE
      END IF
*
*       Check optional sorting of Lagrangian radii & half-mass radius.
      IF (KZ(7).GT.0) THEN
          CALL LAGR(RDENS)
      END IF
*
*       Initialize diagnostic variables.
   10 NP = 0
      IUNP = 0
      AMIN = 100.0
*
*       Find smallest semi-major axis and count unperturbed KS pairs.
      DO 20 IPAIR = 1,NPAIRS
          NP = NP + LIST(1,2*IPAIR-1)
          SEMI = -0.5*BODY(N+IPAIR)/H(IPAIR)
          IF (SEMI.GT.0.0) AMIN = MIN(AMIN,SEMI)
          IF (LIST(1,2*IPAIR-1).EQ.0) IUNP = IUNP + 1
   20 CONTINUE
*
*       Perform time-step & neighbour statistics.
      DTI = 0.0
      DTRI = 0.0
      CNNB = 0.0
      CMAX = 0.0
      NNB = 0
      DO 30 I = IFIRST,NTOT
          DTI = DTI + 1.0/STEP(I)
          DTRI = DTRI + 1.0/STEPR(I)
          CNNB = CNNB + LIST(1,I)/STEP(I)
          RHO = LIST(1,I)/RS(I)**3
          CMAX = MAX(CMAX,RHO)
          NNB = NNB + LIST(1,I)
   30 CONTINUE
*
*       Estimate relative cost & effective neighbour number of AC scheme.
      COST = CNNB/(FLOAT(N - NPAIRS)*DTRI)
      CNNB = CNNB/DTI
*       Scale maximum particle density contrast by the mean value.
      CMAX = 2.0*CMAX*RSCALE**3/FLOAT(N)
*
*       Set average neighbour number & density centre displacement.
      NNB = FLOAT(NNB)/FLOAT(N - NPAIRS)
      RD = SQRT(RDENS(1)**2 + RDENS(2)**2 + RDENS(3)**2)
*
*       Check print frequency indicator & optional model counter.
      NPRINT = NPRINT + 1
      IF (NPRINT.GT.NFIX.OR.TIME.LE.0.0D0) THEN
          NPRINT = 1
          IF (KZ(3).GT.0) MODEL = MODEL + 1
      END IF
*
*       Form binary & merger energy ratios.
      EB = EBIN/(ZKIN - POT)
      EM = EMERGE/(ZKIN - POT)
*
*       Print main output diagnostics.
      I6 = TSCALE*TIME
*
*     IF(rank.EQ.0)THEN
      WRITE (6,40)  rank,
     &     TIME, N, NNB, NPAIRS, NMERGE, NSUB, NSTEPI, NSTEPR,
     &              NSTEPU, ERROR, BE(3)
   40 FORMAT (//,I3,' T =',F7.1,'  N =',I5,'  <NB> =',I3,'  KS =',I3,
     &           '  NM =',I2,'  NS =',I2,'  NSTEPS =',I11,I10,I11,
     &                                     '  DE =',F15.6,'  E =',F17.6)
*     END IF
*
      IF (KZ(21).NE.0) THEN
          CALL CPUTIM(TCOMP)
		  TCOMP = (TCOMP-TTOT)*60.
          IF (VC.EQ.0.0D0) VC = RSCALE/TCR
          TRC = 1.02*FLOAT(NC)**2*BODYM/(VC**3*LOG(FLOAT(NC)))
          DMIN1 = MIN(DMIN1, DMIN2, DMIN3, DMIN4)
*
      IF(rank.EQ.0)THEN
          WRITE (6,45)  NRUN, MODEL, TCOMP, TRC, DMIN1, DMIN2, DMIN3,
     &                  DMIN4, AMIN, RMAX, RSMIN
   45     FORMAT (/,' NRUN =',I3,'  M# =',I3,'  CPU =',F8.1,'  TRC =',
     &                        F5.1, '  DMIN =',1P4E8.1,'  AMIN =',E8.1,
     &                              '  RMAX =',E8.1,'  RSMIN =',0PF5.2)
      END IF
      END IF
*
      IF(rank.EQ.0)THEN
      WRITE (6,50)
   50 FORMAT (/,'    <R>  RTIDE  RDENS   RC    NC   MC   RHOD   RHOM',
     &                    '  CMAX   <Cn>  Ir/R  UN  NP    RCM    VCM',
     &                       '         AZ     EB/E   EM/E   TCR     T6')
*
      WRITE (6,55)  RSCALE, RTIDE, RD, RC, NC, ZMC, RHOD, RHOM, CMAX,
     &              CNNB, COST, IUNP, NP, CMR(4), CMRDOT(4), AZ, EB, EM,
     &              TCR, I6
   55 FORMAT (' #1',F5.2,F6.1,F7.2,F6.2,I5,F7.3,F6.0,F7.0,F6.0,F6.1,
     &                           F6.2,2I4,F9.3,F8.4,F11.6,2F7.3,F6.2,I6)
*
      WRITE (6,60)
   60 FORMAT (/,7X,'NNPRED    NBCORR  NBFULL  NBVOID  NRCONV',
     &         '    NICONV  NBSMIN  NBDIS  NBDIS2  NCMDER  NBDER',
     &           '  NFAST  NBFAST    NBLOCK     NBPRED')
      WRITE (6,65)  NNPRED, NBCORR, NBFULL, NBVOID, NRCONV, NICONV,
     &              NBSMIN, NBDIS, NBDIS2, NCMDER, NBDER, NFAST,
     &              NBFAST, NBLOCK, NBPRED
   65 FORMAT (' #2',I10,I10,3I8,I10,I8,I7,2I8,2I7,I8,I10,I11)
*
      WRITE (6,70)
   70 FORMAT (/,5X,'NKSTRY  NKSREG  NKSHYP     NKSPER  NPRECT  NKSREF',
     &           '  NKSMOD  NTTRY  NTRIP  NQUAD  NCHAIN  NMERG',
     &           '  NSTEPT  NSTEPQ  NSTEPC    NBLCKR')
      WRITE (6,75)  NKSTRY, NKSREG,  NKSHYP, NKSPER, NPRECT, NKSREF,
     &              NKSMOD, NTTRY, NTRIP, NQUAD, NCHAIN, NMERG,
     &              NSTEPT, NSTEPQ, NSTEPC, NBLCKR
   75 FORMAT (' #3',3I8,I11,3I8,3I7,I8,I7,3I8,I10)
      END IF
*
*       Check output for mass loss or tidal capture.
      IF (KZ(19).GT.0.OR.KZ(27).GT.0) THEN
          CALL EVENTS
      END IF
*
*       Reset minimum encounter distances & maximum apocentre separation.
      DMIN2 = 100.0
      DMIN3 = 100.0
      DMIN4 = 100.0
      RSMIN = 100.0
      RMAX = 0.0
*
*       Check integer overflows (2?{32} or 2.1 billion).
      IF (NSTEPI.GT.2000000000) THEN
          NSTEPI = 0
      END IF
      IF (NSTEPU.GT.2000000000) THEN
          NSTEPU = 0
      END IF
      IF (NBPRED.GT.2000000000) THEN
          NBPRED = 0
      END IF
*
*       Ensure NLIST does not become large for block-step version.
      IF (TIME.LT.TBLOCK) TLIST = 0.0
*
*       Exit if error exceeds restart tolerance (TIME < TADJ means no CHECK).
      IF (ABS(ERROR).GT.5.0*QE.AND.TIME.LT.TADJ) GO TO 100
*
*       Check optional analysis & output of KS binaries.
      IF (KZ(8).NE.0.AND.NPAIRS.GT.0) THEN
          CALL BINOUT
      END IF
*
*       Include optional diagnostics of block-steps.
      IF (KZ(34).GT.0.AND.TIME.GT.0.D0) THEN
          CALL LEVELS
      END IF
*
*       Check optional output of single bodies & binaries.
      CALL BODIES
*
*       See whether to write data bank of binary diagnostics on unit 9.
      IF (KZ(8).GE.3.AND.NPAIRS.GT.0) THEN
          CALL BINDAT
      END IF
*
*       Check optional writing of data on unit 3 (frequency NFIX).
      IF (KZ(3).EQ.0.OR.NPRINT.NE.1) GO TO 100
*
*
*       Convert masses, coordinates & velocities to single precision.
      DO 90 I = 1,NTOT
          BODYS(I) = BODY(I)
          RHOS(I) = RHODBL(I)
          XNS(I) = XNDBL(I)
          PHI(I) = PHIDBL(I)
          DO 85 K = 1,3
              XS(K,I) = X(K,I)
              VS(K,I) = XDOT(K,I)
   85     CONTINUE
   90 CONTINUE
*
*       Split into WRITE (3) NTOT & WRITE (3) ..  if disc instead of tape.
#ifdef PARALLEL
      IF(rank.EQ.0)THEN
#endif
*
C This part was edited in Aug.1998, P.Kroupa
C RHOS=mass density in sphere containing 5 nearest stars 
C   (-> Casertano & Hut 1985, ApJ 298, 80)
C XNS=number density in same volume
C PHI=potential at position of particle J
*
        if (BK(3).EQ.0.OR.BK(3).EQ.2) then
          AS(1) = TIME
          AS(2) = FLOAT(NPAIRS)
          AS(3) = RBAR
          AS(4) = ZMBAR
          AS(5) = RTIDE
          AS(6) = TIDAL(4)
          AS(7) = RDENS(1)
          AS(8) = RDENS(2)
          AS(9) = RDENS(3)
          AS(10) = TIME/TCR
          AS(11) = I6
          AS(12) = FLOAT(NZERO)
          AS(13) = RC
          AS(14) = NC
          AS(15) = VC
          AS(16) = RHOM
          AS(17) = CMAX
          AS(18) = RSCALE
          AS(19) = RSMIN
          AS(20) = DMIN1
*
          WRITE (3000+rank)  NTOT, MODEL, NRUN
          WRITE (3000+rank)  (AS(K),K=1,20), 
     &         (BODYS(J),J=1,NTOT),(RHOS(J),J=1,NTOT),(XNS(J),J=1,NTOT),
     &         ((XS(K,J),K=1,3),J=1,NTOT), ((VS(K,J),K=1,3),J=1,NTOT),
     &         (PHI(J),J=1,NTOT),(NAME(J),J=1,NTOT)
        end if
*
C Added for binary star analysis (P.Kroupa, Aug.1998):
        if (BK(3).EQ.1.OR.BK(3).EQ.2) then
*
C Length Scale: l*RBAR=lp, l in model units, lp in physical units [pc]
C Mass scale: m*ZMBAR=mp, mp in Msun
C Time scale: t*TSCALE_pk = tp, tp in Myr
C Velocity scale: v*VSTAR_pk = vp, vp in km/s
C RSCALE is half-mass radius in model units
*
          AS(1) = TIME
          AS(2) = FLOAT(NPAIRS)
          AS(3) = RBAR
          AS(4) = ZMBAR
          AS(5) = RTIDE
          AS(6) = TIDAL(4)
          AS(7) = RDENS(1)
          AS(8) = RDENS(2)
          AS(9) = RDENS(3)
          AS(10) = TIME/TCR
          AS(11) = TSCALE_pk
          AS(12) = VSTAR_pk
          AS(13) = RC
          AS(14) = NC
          AS(15) = VC
          AS(16) = RHOM
          AS(17) = CMAX
          AS(18) = RSCALE
          AS(19) = RSMIN
          AS(20) = DMIN1
*
          NK = 20
          WRITE (40+rank)  NTOT, MODEL, NRUN, NK
          WRITE (40+rank)  (AS(K),K=1,NK), 
     &         (BODY(J),J=1,NTOT),(RHOS(J),J=1,NTOT),
     &         (XNS(J),J=1,NTOT),
     &         ((X(K,J),K=1,3),J=1,NTOT),((XDOT(K,J),K=1,3),J=1,NTOT),
     &         (PHI(J),J=1,NTOT),(NAME(J),J=1,NTOT)
*
          write(6,*)
          write(6,'(a,I3,a,I7)')' Wrote to fort.',40+rank,
     &                          '; model#= ',MODEL
          write(6,'(a,F10.4,a)')' at time= ',TIME*TSCALE_pk,' Myr'
          write(6,*)
        end if
#ifdef PARALLEL
      END IF
#endif
*
*       Update next output interval and initialize the corresponding error.
  100 TNEXT = TNEXT + DELTAT
      ERROR = 0.0D0
*
      RETURN
*
      END


