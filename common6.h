*       COMMON6.
*       -------
*
      INCLUDE 'params.h'
      IMPLICIT REAL*8  (A-H,O-Z)
c      REAL*8         RADIUS,TEV,BODY0
      INTEGER        KVEC,KZ,BK,LISTR,LISTD,LISTV,KSLOW,NLIST,NAME,LIST,
     &               KSTAR,ILIST,JLIST,JPERT
      INCLUDE 'mpif.h'
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      SAVE
      COMMON/MPIDAT/group,rank,ierr,isize,status
*
      COMMON/NBODY/  X(ID,NMAX),X0(ID,NMAX),X0DOT(ID,NMAX),F(ID,NMAX),
     &               XN(ID,NMAX),XNDOT(ID,NMAX),
     &               FDOT(ID,NMAX),BODY(NMAX),RS(NMAX),XDOT(ID,NMAX),
     &               FI(ID,NMAX),D1(ID,NMAX),D2(ID,NMAX),D3(ID,NMAX),
     &               FR(ID,NMAX),D1R(ID,NMAX),D2R(ID,NMAX),D3R(ID,NMAX),
     &               STEP(NMAX),T0(NMAX),STEPR(NMAX),T0R(NMAX)
*
      COMMON/PAIRS/  TAU,HT,U(4,KMAX),U0(4,KMAX),UDOT(4,KMAX),
     &               FU(4,KMAX),D1U(4,KMAX),D2U(4,KMAX),D3U(4,KMAX),
     &               DTAU(KMAX),T0U(KMAX),T1U(KMAX),T2U(KMAX),T3U(KMAX),
     &               HDOT(KMAX),D1HDOT(KMAX),D2HDOT(KMAX),D3HDOT(KMAX),
     &               FUDOT(4,KMAX),TDOT2(KMAX),TDOT3(KMAX),
     &               H(KMAX),R(KMAX),R0(KMAX),GAMMA(KMAX)
*
      COMMON/NAMES/  N,NTOT,NPAIRS,NNBOPT,NCRIT,NFIX,NMERGE,NSUB,NCH,
     &               IPHASE,IFIRST,ICOMP,JCOMP,ICLOSE,JCLOSE,KSPAIR,
     &               NRUN,MODEL,NC,NZERO,NBIN0,IDUM1,NAME0,KZ(40),
     &               BK(10),
     &               KVEC(2*KMAX),LISTR(MLR),LISTD(MLD),LISTV(MLV),
     &               KSLOW(KMAX),NLIST(NMAX),NAME(NMAX),LIST(LMAX,NMAX)
*
      COMMON/COUNTS/ NSTEPI,NSTEPR,NSTEPU,NNPRED,NBCORR,NBFULL,NBVOID,
     &               NRCONV,NBSMIN,NLSMIN,NBDIS,NBDIS2,NCMDER,NBDER,
     &               NFAST,NBFAST,NBLOCK,NBPRED,NICONV,NCHAIN,NSTEPC,
     &               NKSTRY,NKSREG,NKSHYP,NKSPER,NPRECT,NKSREF,NKSMOD,
     &               NTTRY,NTRIP,NQUAD,NMERG,NSTEPT,NSTEPQ,NDISS,NTIDE,
     &               NCOLL,NSYNC,NSESC,NBESC,NMESC,NTIMER,NSTEPS,NPRINT,
     &               NDUMP,NPOP(10),NCOUNT(30),NBLCKR,IPE,NPES,NDUM(2)
*
      COMMON/PARAMS/ CPU,ETAI,ETAR,DTADJ,DELTAT,TCRIT,TCRITp,
     &               QE,RBAR,ZMBAR,
     &               DTMIN,RMIN,ETAU,CMSEP2,ECLOSE,GMIN,GMAX,ETA0,
     &               TWOPI,ONE3,ONE6,ONE9,ONE12,TCR0,TRH,BODYM,BODY1,
     &               SMIN,RMIN2,RMIN22,FCRIT2,ALPHA,ZNBMIN,ZNBMAX,
     &               TADJ,TNEXT,CPU0,CPUTOT,TLIST,DTLIST,STAB1,STAB2,
     &               TIME,ZMASS,RSCALE,TCR,TRC,BE(3),CMR(4),CMRDOT(4),
     &               ZKIN,POT,EBIN,EBIN0,ESUB,EMERGE,ECOLL,EDISS,ESYNC,
     &               E(12),ERROR,ERRTOT,DETOT,ETCORR,AZ,PCRIT,
     &               RTIDE,TSCALE,TSCALE_pk,
     &               TIDAL(4),ETIDE,RSFAC,RSPH2,BETA,
     &               RC,RC2,RC2IN,VC,ZMC,RDENS(3),RHOD,RHOM,RSMIN,RMAX,
     &               DMIN1,DMIN2,DMIN3,DMIN4,DMINC,SBCOLL,BBCOLL,CHCOLL,
     &               DELTAS,ORBITS(9),GPRINT(9),TLASTT,TLASTS,TLASTB(9),
     &               TDUMP,DUMMY(5)
*
      COMMON/STARS/  SIGMA0,RSYNC,EPOCH,TSTAR,VSTAR,VSTAR_pk,
     &               ZMRG,ZMHE,ZMRS,ZMWD,
     &               ZMSN,ZMDOT,EMDOT,TMDOT,TPHYS,DUMMY1(4),
     &               RADIUS(NMAX),TEV(NMAX),BODY0(500),
     &               NMDOT,NRG,NHE,NRS,NWD,NSN,NTZ,NBS,NDUM2(10),
     &               NTYPE(10),KSTAR(NMAX)
*
      COMMON/HERMIT/ FIDOT(ID,NMAX),D0(ID,NMAX),FRDOT(ID,NMAX),
     &               D0R(ID,NMAX),TIMENW(NMAX)
*
      COMMON/BLOCKS/ TPREV,TBLOCK,DTK(64)
*
      COMMON/LISTS/  ILIST(NMAX),JLIST(NMAX),JPERT(LMAX)
*        Common block to keep neighbour density and potential high prec (R.Sp.)
      COMMON/WORK2/RHODBL(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
*
      COMMON/TIMING/ttot,ttreg,ttirr,ttpre,ttinit,ttint,ttks,
     *  ttcomm,ttadj,ttmov,ttnbp,ttsub,ttsub2,ttfrc,
     *  isernb,iserreg











