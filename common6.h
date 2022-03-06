*       COMMON6.
*       -------
*
      INCLUDE 'params.h'
      IMPLICIT REAL*8  (A-H,O-Z)
      INTEGER BK
*
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
     &               FRG(ID,NMAX), FRGDOT(ID,NMAX),
     &               STEP(NMAX),T0(NMAX),STEPR(NMAX),T0R(NMAX)
*
      COMMON/PAIRS/  HT,U(4,KMAX),U0(4,KMAX),UDOT(4,KMAX),FU(4,KMAX),
     &               FUDOT(4,KMAX),FUDOT2(4,KMAX),FUDOT3(4,KMAX),
     &               H(KMAX),HDOT(KMAX),HDOT2(KMAX),HDOT3(KMAX),
     &               HDOT4(KMAX),DTAU(KMAX),TDOT2(KMAX),TDOT3(KMAX),
     &               R(KMAX),R0(KMAX),GAMMA(KMAX),SF(7,KMAX),H0(KMAX),
     &               FP0(4,KMAX),FD0(4,KMAX),TBLIST,DTB,KBLIST(10*KMAX),
     &               NNTB
*
      COMMON/NAMES/  N,NTOT,NPAIRS,NNBOPT,NCRIT,NFIX,NMERGE,NSUB,NCH,
     &               IPHASE,IFIRST,ICOMP,JCOMP,ICLOSE,JCLOSE,KSPAIR,
     &               NRUN,MODEL,NC,NZERO,NBIN0,IDUM1,NAME0,NHI0,KZ(40),
     &               BK(10),LSHORT(NMAX),NNBMAX,
     &               KVEC(2*KMAX),LISTR(MLR),LISTD(MLD),LISTV(MLV),
     &               KSLOW(KMAX),NLIST(NMAX),NAME(NMAX),LIST(LMAX,NMAX)
*
      COMMON/COUNTS/ NSTEPI,NSTEPR,NSTEPU,NNPRED,NBCORR,NBFULL,NBVOID,
     &               NRCONV,NBSMIN,NLSMIN,NBDIS,NBDIS2,NCMDER,NBDER,
     &               NFAST,NBFAST,NBLOCK,NBPRED,NICONV,NCHAIN,NSTEPC,
     &               NKSTRY,NKSREG,NKSHYP,NKSPER,NPRECT,NKSREF,NKSMOD,
     &               NTTRY,NTRIP,NQUAD,NMERG,NSTEPT,NSTEPQ,NDISS,NTIDE,
     &               NCOLL,NSYNC,NSESC,NBESC,NMESC,NTIMER,NSTEPS,NPRINT,
     &               NDUMP,NBPREV,NPOP(10),JCMAX,NEWHI,NBLCKR,IPE,NPES,
     &               NSTEPB,NBFLUX,NDUMMY(8)
*
      COMMON/PARAMS/ CPU,ETAI,ETAR,DTADJ,DELTAT,TCRIT,TCRITp,
     &               QE,RBAR,ZMBAR,
     &               DTMIN,RMIN,ETAU,CMSEP2,ECLOSE,GMIN,GMAX,ETA0,
     &               TWOPI,ONE3,ONE6,ONE9,ONE12,TCR0,TRH,BODYM,BODY1,
     &               SMIN,RMIN2,RMIN22,FCRIT2,ALPHA,ZNBMIN,ZNBMAX,EBH,
     &               TADJ,TNEXT,CPU0,CPUTOT,TLIST,DTLIST,STAB1,STAB2,
     &               TIME,ZMASS,RSCALE,TCR,TRC,BE(3),CMR(4),CMRDOT(4),
     &               ZKIN,POT,EBIN,EBIN0,ESUB,EMERGE,ECOLL,EDISS,ESYNC,
     &               E(12),ERROR,ERRTOT,DETOT,ETCORR,AZ,PCRIT,EBCH0,
     &               RTIDE,TSCALE,TIDAL(4),ETIDE,RSFAC,RSPH2,BETA,
     &               RC,RC2,RC2IN,VC,ZMC,RDENS(3),RHOD,RHOM,RSMIN,RMAX,
     &               DMIN1,DMIN2,DMIN3,DMIN4,DMINC,SBCOLL,BBCOLL,CHCOLL,
     &               DELTAS,ORBITS(9),GPRINT(9),TLASTT,TLASTS,TLASTB(9),
     &               TDUMP,SCOEFF(12),TOFF,TTOT,DUMMY
*
      COMMON/STARS/  EPOCH0,ZMRG,ZMHE,ZMRS,ZMWD,ZMSN,ZMNH,ZMBH,ZMDOT,
     &               AU,PC,GM,DAYS,YRS,SU,SMU,RAU,TSTAR,VSTAR,STEPX,
     &               TMDOT,TPHYS,TURN,EMDOT,ECDOT,EKICK,TPLOT,DUMMY1,
     &               XHYD,YHEL,ZMET,ZPARS(20),
     &               RADIUS(NMAX),TEV(NMAX),TEV0(NMAX),BODY0(NMAX),
     &               EPOCH(NMAX),NMDOT,NRG,NHE,NRS,NWD,NNH,NSN,NBS,NTZ,
     &               NBH,NKICK,NBKICK,KSAVE(2),ITYPE(10),KSTAR(NMAX)
*
      COMMON/HERMIT/ FIDOT(ID,NMAX),D0(ID,NMAX),FRDOT(ID,NMAX),
     &               D0R(ID,NMAX),TIMENW(NMAX)
*
      COMMON/BLOCKS/ TPREV,TBLOCK,DTK(64)
*
      COMMON/LISTS/  ILIST(NMAX),JLIST(NMAX),JPERT(5*LMAX)
*        Common block to keep neighbour density and potential high prec (R.Sp.)
      COMMON/WORK2/RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
*
      COMMON/TIMING/ttota,ttreg,ttirr,ttpre,ttinit,ttint,ttks,
     *  ttcomm,ttadj,ttmov,ttnbp,ttsub,ttsub2,ttfrc, ttgrcomm, ttgrcalc,
     *  isernb,iserreg,xtsub1,xtsub2, ipipe, np_use, npmax











