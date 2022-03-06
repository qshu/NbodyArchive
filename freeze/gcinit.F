      SUBROUTINE GCINIT
*
*      Galactic Centre and Bulge Forces
*      Dynamical Friction    
*      Initialization Routine
*      R.Sp./ O.G. Nov. 2001
*
      INCLUDE 'common6.h'
      COMMON /GC/ XG0(3),VG0(3),AG0(3),ADOTG0(3),
     *      XG(3),VG(3),AG(3),ADOTG(3),AG2(3),AG3(3),
     *      XGI(3),VGI(3),AGI(3),ADOTGI(3),FTOT(3),FDTOT(3),
     *      ADF0(3),ADOTF0(3),ADF(3),ADOTF(3),
     *      XMBH,CC,CFAC,XCOUL,XMCL,TCEN0,STEPCC,EG0,TGCL,
     *      DEDFF,ECLS,ECLS0,ETAG,RGAL,VGAL,XMGAL,IMEM(NMAX)
      REAL*8 FREG(3),FD(3),XLG(3),STEPOL(NMAX)
      DATA PI/3.1415927D0/
*
#ifdef PARALLEL
#define MPIINIT 1
#else
#ifdef ENSEMBLE
#define MPIINIT 1
#else
#define MPIINIT 0
#endif
#endif
*       Read & print the cluster, dy.fr. and gal. centre parameters.
         if(rank.eq.0)then
         OPEN (55,STATUS='OLD',FILE='galcen.dat',ERR=999)
         READ (55,*) ETAG,XMBH,XMGAL,VCIRC,RGAL,CFAC,GCOUL
         end if
*
#if MPIINIT
      CALL MPI_BCAST(ETAG,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(XMBH,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(XMGAL,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(VCIRC,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RGAL,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(CFAC,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(GCOUL,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif
*
*      Normalize masses to total initial cluster mass (=1.0D0)
*
      TINY = 1.D-30
*
      XMCTOT = ZMASS*ZMBAR
      XMBH = XMBH*ZMASS
      XMGAL = XMGAL*ZMASS
      RGAL = RGAL*RSCALE
      CC = XMGAL*RGAL*RBAR
      XCOUL = DLOG(GCOUL*N)
      XMCL = ZMASS
      VGAL = VCIRC*DSQRT(XMBH/RGAL+CC)
*
      if(rank.eq.0)then
      WRITE (6,10)
   10    FORMAT (/,5X,'ETAG',9X,'XMCL',9X,'XMBH',8X,' XMGAL(at RG)',
     *   7X,'VG',11X,'RG',9X,'CFAC',8X,'GCOUL ')
      WRITE (6,12) ETAG,ZMASS,XMBH,CC,VGAL,RGAL,CFAC,GCOUL
   12    FORMAT (1P8E13.2)
      WRITE (6,11)
   11    FORMAT (8X,' [Msol]: XMCL',5X,'XMBH',8X,' XMGAL(at RG)',
     *    1X,'VG[km/s]',6X,'RG[pc]',8X,'CFAC',8X,'ln(gN)')
      WRITE (6,14) XMCTOT,XMBH*ZMBAR,CC*ZMBAR,VGAL*VSTAR,
     *     RGAL*RBAR,CFAC,DLOG(GCOUL*N)
   14    FORMAT (13X,1P8E13.2)
      end if
*
*     Initialize Time and Positions/Velocities
*
      TCEN0 = 0.D0
      DR2I = 1.D0/RGAL**2
*
      XG(1) = RGAL
      XG(2) = 0.D0
      XG(3) = 0.D0
      VG(1) = 0.D0
      VG(2) = VGAL
      VG(3) = 0.D0
*      Force and derivative at cluster centre due to gravity of galactic centre
*
       RIJ2 = XG(1)**2 + XG(2)**2 + XG(3)**2
       DR2I = 1.0D0/RIJ2
       DRDV = XG(1)*VG(1) + XG(2)*VG(2) + XG(3)*VG(3)
*
*      First: central black hole
*
       DR3I = XMBH*DR2I*DSQRT(DR2I)
       DRDP = 3.0D0*DRDV*DR2I
*
       DO 2 K = 1,3
          AG(K) = -XG(K)*DR3I
          ADOTG(K) = -(VG(K) - XG(K)*DRDP)*DR3I
 2     CONTINUE
*
*      Second: bulge
*
       DR3I = CC*DR2I
       DRDP = 2.0D0*DRDV*DR2I
*
       DO 3 K = 1,3
          AG(K) = AG(K) - XG(K)*DR3I
          ADOTG(K) = ADOTG(K) - (VG(K) - XG(K)*DRDP)*DR3I
 3     CONTINUE
*
*      Force and derivative cluster centre due to dynamical friction
*
       RIJ2 = XG(1)**2 + XG(2)**2 + XG(3)**2
       VIJ2 = VG(1)**2 + VG(2)**2 + VG(3)**2
       DR2I = 1.0D0/RIJ2
       DV2I = 1.0D0/MAX(VIJ2,TINY)
       DV3I = DV2I*DSQRT(DV2I)
       DFAC = -CFAC*XCOUL*XMCL*DR2I*DV3I*VIJ2
       DRDV = XG(1)*VG(1) + XG(2)*VG(2) + XG(3)*VG(3)
       DADV = AG(1)*VG(1) + AG(2)*VG(2) + AG(3)*VG(3)
       PRINT*,rank,' Initial TDF=',1.D0/DFAC,' TORB=',2.D0*PI*RGAL/VGAL
*
       DRDP = DADV*DV2I
*
       DO 4 K = 1,3
*
       ADF(K) = DFAC*VG(K)
       ADOTF(K) = DFAC*(AG(K) - VG(K)*DADV*DV2I
     *                                  - 2.D0*VG(K)*DRDV*DR2I)
 4     CONTINUE
*
      DO 6 K = 1,3
      XG0(K) = XG(K)
      VG0(K) = VG(K)
      AG0(K) = AG(K)
      ADOTG0(K) = ADOTG(K)
      ADF0(K) = ADF(K)
      ADOTF0(K) = ADOTF(K)
      FTOT(K) = AG(K) + ADF(K)
      FDTOT(K) = ADOTG(K) + ADOTF(K)
 6    CONTINUE
*
      FI2 = FTOT(1)**2 + FTOT(2)**2 + FTOT(3)**2
      FD2 = FDTOT(1)**2 + FDTOT(2)**2 + FDTOT(3)**2
*
      DTCEN = ETAG**2*DSQRT(FI2/MAX(FD2,TINY))
*
      PRINT*,rank,' Call STEPK with DTCEN=',DTCEN,(FTOT(K),K=1,3)
      CALL FLUSH(6)
*       Convert predicted step to nearest block time-step (truncated down).
          CALL STEPK(DTCEN,DTN)
*
      PRINT*,rank,' Return STEPK with DTN=',DTN
      CALL FLUSH(6)
*
      STEPCC = DTN
*
      XXG = DSQRT(XG(1)**2 + XG(2)**2 + XG(3)**2)
      XVG = DSQRT(VG(1)**2 + VG(2)**2 + VG(3)**2)
      XAG = DSQRT(AG(1)**2 + AG(2)**2 + AG(3)**2)
      XADOTG = DSQRT(ADOTG(1)**2 + ADOTG(2)**2 + ADOTG(3)**2)
      XADF = DSQRT(ADF(1)**2 + ADF(2)**2 + ADF(3)**2)
      XADOTF = DSQRT(ADOTF(1)**2 + ADOTF(2)**2 + ADOTF(3)**2)
*
      WRITE(77,177)TIME,TCEN0,STEPCC,XXG,XVG,XAG,XADOTG,XADF,XADOTF,
     *    XG(1),XG(2),XG(3)
 177  FORMAT(1X,1P,9D13.3,3D17.5)
      WRITE(78,177)TIME,(AG(K),K=1,3),(ADOTG(K),K=1,3)
*
      PRINT*,rank,' AG=',(AG(K),K=1,3),' ||=',
     *                  DSQRT(AG(1)**2+AG(2)**2+AG(3)**2)
      PRINT*,rank,' ADOTG=',(ADOTG(K),K=1,3),' ||=',
     *                  DSQRT(ADOTG(1)**2+ADOTG(2)**2+ADOTG(3)**2)
      PRINT*,rank,' ADF=',(ADF(K),K=1,3),' ||=',
     *                  DSQRT(ADF(1)**2+ADF(2)**2+ADF(3)**2)
      PRINT*,rank,' ADOTF=',(ADOTF(K),K=1,3),' ||=',
     *                  DSQRT(ADOTF(1)**2+ADOTF(2)**2+ADOTF(3)**2)
*
*      Include regular forces for all particles and timesteps
*
       DTMAXC = 0.D0
       DTMINC = 0.D0
       ILARGE = 0
       ISMALL = 0
       ITSTEP = 0
       DEDF1 = 0.0D0
       DEDF2 = 0.0D0
       DEDF3 = 0.0D0
       DEDF4 = 0.0D0
       DEDFF = 0.0D0
       TGCL = 0.0D0
*
       DO 200 I = 1,N
*
*      Relative positions and velocities to g.c. using pred. pos. of particle
*
       DO 201 K = 1,3
          XGI(K) = XG(K) + X(K,I)
          VGI(K) = VG(K) + XDOT(K,I)
 201   CONTINUE
*
       RIJ2 = XGI(1)**2 + XGI(2)**2 + XGI(3)**2
       DR2I = 1.0/RIJ2
       DRDV = XGI(1)*VGI(1) + XGI(2)*VGI(2) + XGI(3)*VGI(3)
*
*      First: central black hole
*
       DR3I = XMBH*DR2I*DSQRT(DR2I)
       DRDP = 3.0D0*DRDV*DR2I
*
*      Force and derivative at particle position due to g.c.
*
*
       DO 202 K = 1,3
          AGI(K) = -XGI(K)*DR3I
          ADOTGI(K) = -(VGI(K) - XGI(K)*DRDP)*DR3I
 202   CONTINUE
*
*      Second: bulge
*
       DR3I = CC*DR2I
       DRDP = 2.0D0*DRDV*DR2I
*
*      Force and derivative at particle position due to g.c.
*
       DO 203 K = 1,3
          AGI(K) = AGI(K) - XGI(K)*DR3I
          ADOTGI(K) = ADOTGI(K) - (VGI(K) - XGI(K)*DRDP)*DR3I
 203   CONTINUE
*
*      Subtraction of apparent (pred.) force on cluster centre
*
       DO 204 K = 1,3
          AGI(K) = AGI(K) - AG(K)
          ADOTGI(K) = ADOTGI(K) - ADOTG(K)
 204   CONTINUE
*
*       Obtain old time-steps using DT = ETA*F/FDOT (re-instated 9/98).
          FR2 = FR(1,I)**2 + FR(2,I)**2 + FR(3,I)**2
          FDR2 = D1R(1,I)**2 + D1R(2,I)**2 + D1R(3,I)**2
*       Include precaution for small velocities (i.e. DT = ETA*TCR)
          IF (FDR2.LT.FR2) FDR2 = FR2/TCR**2
          DTROLD = ETAR*SQRT(FR2/FDR2)
*
       DO 205 K = 1,3
*      Subtract external forces (pred.) for non-cluster members
*         IF(IMEM(I).EQ.0)THEN
*            AGI(K) = AGI(K) - ADF(K)
*            ADOTGI(K) = ADOTGI(K) - ADOTF(K)
*         END IF
*
          FR(K,I) = FR(K,I) + AGI(K)
          FRDOT(K,I) = FRDOT(K,I) + ADOTGI(K)
          D1R(K,I) = FRDOT(K,I)
*
 205   CONTINUE
*
*       Obtain new time-steps using DT = ETA*F/FDOT (re-instated 9/98).
          FR2 = FR(1,I)**2 + FR(2,I)**2 + FR(3,I)**2
          FDR2 = D1R(1,I)**2 + D1R(2,I)**2 + D1R(3,I)**2
*       Include precaution for small velocities (i.e. DT = ETA*TCR)
          IF (FDR2.LT.FR2) FDR2 = FR2/TCR**2
          DTR = ETAR*ETAR*SQRT(FR2/FDR2)
*       Determine irregular and regular steps by the general criterion.
*         DTR = TSTEP(FR(1,I),D1R(1,I),D2R(1,I),D3R(1,I),ETAR)
*
          DTREL = DABS(DTR-DTROLD)/DTROLD
*
          IF(DTR.GT.DTROLD)THEN
              ILARGE = ILARGE + 1
              IF(DTREL.GT.DTMAXC)THEN
                  IPARTP = I
                  DTMAXC = DTREL
              END IF
          ELSE
              ISMALL = ISMALL + 1
              IF(DTREL.GT.DTMINC)THEN
                  IPARTM = I
                  DTMINC = DTREL
              END IF
          END IF
*
*       Convert predicted step to nearest block time-step (truncated down).
          CALL STEPK(DTR,DTN)
*
          STEPR(I) = DTN
          IF (STEP(I).GT.STEPR(I))THEN
              ITSTEP = ITSTEP + 1
              STEP(I) = STEPR(I)
          END IF
*
*      Start Energy Checks
         RIJ2 = XGI(1)**2 + XGI(2)**2 + XGI(3)**2
         DR2I = 1.0/RIJ2
         DR3I = XMBH*DR2I*DSQRT(DR2I)
          PHIEXT1 = -DR3I*RIJ2
          PHIEXT2 = CC*DLOG(DSQRT(RIJ2)/RGAL)
*
          DEDF3 = DEDF3 + BODY(I)*PHIEXT1
          DEDF4 = DEDF4 + BODY(I)*PHIEXT2
          DEDF1 = DEDF1 + BODY(I)*(VG(1)**2+VG(2)**2+VG(3)**2)/2.D0
          DEDF2 = DEDF2 + BODY(I)*(VG(1)*XDOT(1,I) + VG(2)*XDOT(2,I) +
     *            VG(3)*XDOT(3,I))
*
 200   CONTINUE
*
       ECLK = XMCL*(VG(1)**2+VG(2)**2+VG(3)**2)/2.D0
       RIJ2 = XG(1)**2 + XG(2)**2 + XG(3)**2
         DR2I = 1.0/RIJ2
         DR3I = XMBH*DR2I*DSQRT(DR2I)
         PHICLS = -DR3I*RIJ2 + CC*DLOG(DSQRT(RIJ2)/RGAL)
*
       ECLS0 = ECLK + XMCL*PHICLS
*
          PRINT*,rank,' ITSTEP=',ITSTEP
          EG0 = DEDF1 + DEDF2 + DEDF3 + DEDF4
*
      PRINT*,rank,' Time Step Situation after ',
     *   ' inclusion of external field: '
      PRINT*,rank,' Max down Change of DTR =',DTMINC,' for ',
     *    IPARTM,' ISMALL =',ISMALL
      PRINT*,rank,' Max up   Change of DTR =',DTMAXC,' for ',
     *    IPARTP,' ILARGE =',ILARGE
      PRINT*,rank,' Initial Energies DEDF1-4=',DEDF1,DEDF2,DEDF3,DEDF4,
     *  ' tot =',EG0,' ECLS=',ECLS0
*
       DEDF1 = 0.0D0
       DEDF2 = 0.0D0
       DEDF3 = 0.0D0
       DEDF4 = 0.0D0
*
      RETURN
*
 999  PRINT*,rank,' Error in Reading file unit 55, file name galcen.dat'
      STOP
*
      END
