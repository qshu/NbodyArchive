      SUBROUTINE OUTPUT
*
*
*       Output and data save.
*       ---------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX)
      COMMON/ECHAIN/  ECH
      REAL*8  X1(3,4),V1(3,4),UI(4),VI(4),XREL2(3),VREL2(3)
      REAL*4  XS(3,NMAX),VS(3,NMAX),BODYS(NMAX),RHOS(NMAX),AS(20)
      REAL*4  XNS(NMAX),PHI(NMAX)
      LOGICAL  FIRST
      SAVE  FIRST
      DATA  FIRST /.TRUE./
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
      ETOT = ZKIN - POT + ETIDE + EBIN + ESUB + EMERGE + ECOLL + EMDOT
     &                                                         + ECDOT
      IF (NCH.GT.0) THEN
          ETOT = ETOT + ECH
      END IF
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
      MULT = 0
*
*       Find smallest semi-major axis and count unperturbed KS pairs.
      DO 20 IPAIR = 1,NPAIRS
          NP = NP + LIST(1,2*IPAIR-1)
          SEMI = -0.5*BODY(N+IPAIR)/H(IPAIR)
          IF (SEMI.GT.0.0) AMIN = MIN(AMIN,SEMI)
          IF (LIST(1,2*IPAIR-1).EQ.0) IUNP = IUNP + 1
          IF (NAME(N+IPAIR).LT.-2*NZERO) MULT = MULT + 1
   20 CONTINUE
*
*       Perform time-step & neighbour statistics (NS is # single stars).
      DTI = 0.0
      DTRI = 0.0
      CNNB = 0.0
      CMAX = 0.0
      NNB = 0
      NS = 0
      DO 30 I = IFIRST,NTOT
          DTI = DTI + 1.0/STEP(I)
          DTRI = DTRI + 1.0/STEPR(I)
          CNNB = CNNB + LIST(1,I)/STEP(I)
          RHON = LIST(1,I)/RS(I)**3
          CMAX = MAX(CMAX,RHON)
          NNB = NNB + LIST(1,I)
          IF (I.LE.N.AND.BODY(I).GT.0.0D0) NS = NS + 1
   30 CONTINUE
      NS = NS - NSUB
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
      IF(rank.EQ.0)THEN
      WRITE (6,40)  rank,
     &     TIME, N, NNB, NPAIRS, NMERGE, MULT, NS, NSTEPI, NSTEPB,
     &              NSTEPR, NSTEPU, ERROR, BE(3)
   40 FORMAT (//,I3,' T =',F7.1,'  N =',I5,'  <NB> =',I3,'  KS =',I5,
     &       '  NM =',I2,' MM =',I2,' NS =',I6,' NSTEPS =',
     &          I11,2I10,I11,'  DE =',E15.6,'  E =',F17.6)
      END IF
*
      IF (KZ(21).GT.0) THEN
          CALL CPUTIM(TCOMP)
		  TCOMP = (TCOMP-TTOT)*60.
          IF (VC.EQ.0.0D0) VC = RSCALE/TCR
          TRC = 1.02*FLOAT(NC)**2*BODYM/(VC**3*LOG(FLOAT(NC)))
          DMIN1 = MIN(DMIN1, DMIN2, DMIN3, DMIN4, DMINC)
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
   50 FORMAT (/,'    <R>  RTIDE  RDENS   RC      NC   MC   RHOD   RHOM',
     &          '  CMAX   <Cn>  Ir/R    UN    NP    RCM    VCM',
     &          '         AZ     EB/E   EM/E   TCR     T6')
*
      WRITE (6,55)  RSCALE, RTIDE, RD, RC, NC, ZMC, RHOD, RHOM, CMAX,
     &              CNNB, COST, IUNP, NP, CMR(4), CMRDOT(4), AZ, EB, EM,
     &              TCR, I6
   55 FORMAT (' #1',F5.2,F6.1,F7.2,F6.2,I7,F7.3,F6.0,F7.0,F6.0,F6.1,
     &                           F6.2,2I6,F9.3,F8.4,F11.6,2F7.3,F6.2,I6)
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
      DMINC = 100.0
      RSMIN = 100.0
      RMAX = 0.0
*
*       Check integer overflows (2^{32} or 2.1 billion).
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
      IF (KZ(8).GT.0.AND.NPAIRS.GT.0) THEN
          CALL BINOUT
      END IF
*
*       Include optional diagnostics of block-steps.
      IF (KZ(33).GT.0.AND.TIME.GT.0.D0) THEN
          CALL LEVELS
      END IF
*
*       Check optional output of single bodies & binaries.
      CALL BODIES
*
*       See whether to write data bank of binary diagnostics on unit 9.
      IF (KZ(8).GE.3.AND.NPAIRS.GT.0) THEN
          CALL BINDAT
          IF (KZ(8).GT.3) THEN
              CALL HIDAT
          END IF
      END IF
*
*       Check optional diagnostics of evolving stars.
      IF (KZ(12).GT.0.AND.TIME.GT.TPLOT) THEN
          CALL HRPLOT
      END IF
*
*       Check optional writing of data on unit 3 (frequency NFIX).
      IF (KZ(3).EQ.0.OR.NPRINT.NE.1) GO TO 100
*
      DO 80 IPAIR = 1,NPAIRS
          RHO(2*IPAIR-1) = RHO(N+IPAIR)
          RHO(2*IPAIR) = RHO(N+IPAIR)
 80   CONTINUE 
*
*       Convert masses, coordinates & velocities to single precision.
      DO 90 I = 1,NTOT
          BODYS(I) = BODY(I)
          RHOS(I) = RHO(I)
          XNS(I) = XNDBL(I)
          PHI(I) = PHIDBL(I)
          DO 85 K = 1,3
              XS(K,I) = X(K,I)
              VS(K,I) = XDOT(K,I)
   85     CONTINUE
   90 CONTINUE
*
*       Replace any ghosts by actual M, R & V (including 2 binaries).
      DO 95 JPAIR = 1,NPAIRS
          J2 = 2*JPAIR
          J1 = J2 - 1
          ICM = N + JPAIR
*       Determine merger & ghost index for negative c.m. name.
          IF (NAME(ICM).LT.0.AND.BODY(ICM).GT.0.0) THEN
              CALL FINDJ(J1,J,IM)
*       Note: J is ghost index and IM is merger index.
              IF (J.LE.0) GO TO 95
              BODYS(J1) = CM(1,IM)
              BODYS(J) = CM(2,IM)
              ZMB = CM(1,IM) + CM(2,IM)
*       Form global coordinates and velocities from c.m. with XREL & VREL.
              DO K = 1,3
                  X1(K,1) = X(K,J1) + CM(2,IM)*XREL(K,IM)/ZMB
                  X1(K,2) = X(K,J1) - CM(1,IM)*XREL(K,IM)/ZMB
                  V1(K,1) = XDOT(K,J1) + CM(2,IM)*VREL(K,IM)/ZMB
                  V1(K,2) = XDOT(K,J1) - CM(1,IM)*VREL(K,IM)/ZMB
*
                  XS(K,J1) = X1(K,1)
                  XS(K,J)  = X1(K,2)
                  VS(K,J1) = V1(K,1)
                  VS(K,J)  = V1(K,2)
              END DO
*       Look for ghosts of possible second (i.e. outer) merged binary.
              IF (NAMEG(J).GT.NZERO) THEN
                  ICM2 = 0
                  DO  JJ = N+1,NTOT
                      IF (NAME(JJ).EQ.NAMEG(J)) ICM2 = JJ
                  END DO
*       Treat the second binary using inactive KS variables.
                  IF (ICM2.GT.0) THEN
                      IPAIR = ICM2 - N
                      I1 = 2*IPAIR - 1
                      I2 = I1 + 1
                      BODYS(I1) = CM(3,IM)
                      BODYS(I2) = CM(4,IM)
*       Copy KS variables to local scalars.
                      DO K = 1,4
                          UI(K) = U(K,IPAIR)
                          VI(K) = UDOT(K,IPAIR)
                      END DO
*       Transform to physical variables and multiply by 4 (momentum formula).
                      CALL KSPHYS(UI,VI,XREL2,VREL2)
                      ZM = CM(3,IM) + CM(4,IM)
                      DO K = 1,3
                          VREL2(K) = 4.0*VREL2(K)
                          X1(K,3) = X(K,J2) + CM(4,IM)*XREL2(K)/ZM
                          X1(K,4) = X(K,J2) - CM(3,IM)*XREL2(K)/ZM
                          V1(K,3) = XDOT(K,J2) + CM(4,IM)*VREL2(K)/ZM
                          V1(K,4) = XDOT(K,J2) - CM(3,IM)*VREL2(K)/ZM
*
                          XS(K,I1) = X1(K,3)
                          XS(K,I2)  = X1(K,4)
                          VS(K,I1) = V1(K,3)
                          VS(K,I2)  = V1(K,4)
                          XS(K,ICM2) = X(K,J2)
                          VS(K,ICM2) = XDOT(K,J2)
                      END DO
                  END IF
              END IF
          END IF
   95 CONTINUE
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
          NK = 20
*
      IF (FIRST) THEN
          OPEN (UNIT=3,STATUS='NEW',FORM='UNFORMATTED',FILE='OUT3')
          FIRST = .FALSE.
      END IF
*
          WRITE (3)  NTOT, MODEL, NRUN, NK
          WRITE (3)  (AS(K),K=1,NK),
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
C Time scale: t*TSCALE = tp, tp in Myr
C Velocity scale: v*VSTAR = vp, vp in km/s
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
          AS(11) = TSCALE
          AS(12) = VSTAR
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
      IF (FIRST) THEN
          OPEN (UNIT=4,STATUS='NEW',FORM='UNFORMATTED',FILE='OUT4')
          FIRST = .FALSE.
      END IF
          WRITE (4)  NTOT, MODEL, NRUN, NK
          WRITE (4)  (AS(K),K=1,NK), 
     &         (BODY(J),J=1,NTOT),(RHOS(J),J=1,NTOT),
     &         (XNS(J),J=1,NTOT),
     &         ((X(K,J),K=1,3),J=1,NTOT),((XDOT(K,J),K=1,3),J=1,NTOT),
     &         (PHI(J),J=1,NTOT),(NAME(J),J=1,NTOT)
*
          write(6,*)
          write(6,'(a,I3,a,I7)')' Wrote to fort.',4,
     &                          '; model#= ',MODEL
          write(6,'(a,F10.4,a)')' at time= ',TIME*TSCALE,' Myr'
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
