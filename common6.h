*       common6.
*       -------
*
      INCLUDE 'params.h'
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8  MP,MP0,MPDOT
*
      INCLUDE 'mpif.rainer.h'
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      COMMON/MPIDAT/group,rank,ierr,isize,status
*     SAVE
*
      COMMON/NAMES/  NTOT,NPAIRS,NTTOT,N,NNBMAX,NCRIT,NFIX,NMERGE,NSUB,
     &               IPHASE,IFIRST,ICOMP,JCOMP,ICLOSE,JCLOSE,JCMAX,
     &               KSPAIR,NRUN,MODEL,NC,NZERO,NBZERO,NBIN0,NHI0,
     &               NAME0,NCH,NCHAOS,IDUM1,KZ(50),NNBOPT,NEW2(8)

      COMMON/PARAMS/ ETAI,ETAR,DTADJ,DELTAT,TCRIT,QE,RBAR,ZMBAR,
     &               DTMIN,RMIN,ETAU,CMSEP2,ECLOSE,GMIN,GMAX,ETA0,
     &               TWOPI,ONE3,ONE6,ONE9,ONE12,TCR0,TRH,BODYM,BODY1,
     &               SMIN,RMIN2,RMIN22,STEPJ,ALPHA,ZNBMIN,ZNBMAX,EBH,
     &               TIME,TADJ,TNEXT,CPU,CPU0,CPUTOT,ZMASS,RSCALE,TCR,
     &               TRC,BE(3),CMR(4),CMRDOT(4),ZKIN,POT,EBIN,EBIN0,
     &               ESUB,EMERGE,ECOLL,EDISS,ESYNC,E(12),ERROR,ERRTOT,
     &               DETOT,ETCORR,AZ,PCRIT,EBCH0,RTIDE,TSCALE,TIDAL(4),
     &               HT,ETIDE,EGRAV,RSFAC,RSPH2,RC,RC2,RC2IN,VC,ZMC,
     &               RDENS(3),RHOD,RHOM,RSMIN,RMAX,DMIN1,DMIN2,DMIN3,
     &               DMIN4,DMINC,SBCOLL,BBCOLL,CHCOLL,DELTAS,ORBITS(9),
     &               GPRINT(9),TLASTT,TLASTS,TLASTB(9),TDUMP,
     &               SCOEFF(12),TOFF,TTOT,EBESC,EMESC,ESESC,CLIGHT,
     &               RZ,TINY,SMAX,WTOT,WTOT0,TCRITp,DUMMY(97)

      COMMON/COUNTS/ NSTEPI,NSTEPR,NSTEPU,NNPRED,NBCORR,NBFULL,NBVOID,
     &               NNTB,NBSMIN,NLSMIN,NBDIS,NBDIS2,NCMDER,NBDER,
     &               NFAST,NBFAST,NBLOCK,NBPRED,NICONV,NCHAIN,NSTEPC,
     &               NKSTRY,NKSREG,NKSHYP,NKSPER,NPRECT,NEWKS,NKSMOD,
     &               NTTRY,NTRIP,NQUAD,NMERG,NSTEPT,NSTEPQ,NDISS,NTIDE,
     &               NCOLL,NSYNC,NSESC,NBESC,NMESC,NTIMER,NSTEPS,NPRINT,
     &               NDUMP,NBPREV,NEWHI,NSTEPB,NBFLUX,NMTRY,NWARN,
     &               NIRECT,NURECT,NBRECT,NRRECT,KSMAG,NOFL(2),NPOP(10),
     &               NBLCKR,NDUMMY(99)
 
      COMMON/PLPOT/  MP,AP2,VIR,MP0,MPDOT,TDELAY,RTIDE0,QVIR,PLDUM(4)

      COMMON/BLOCKS/ TPREV,TBLOCK,DTK(64),KVEC(2*KMAX)

      COMMON/STARS/  EPOCH0,ZMRG,ZMHE,ZMRS,ZMWD,ZMSN,ZMNH,ZMBH,ZMDOT,
     &               AU,PC,GM,DAYS,YRS,SU,SMU,RAU,TSTAR,VSTAR,STEPX,
     &               TMDOT,TPHYS,TURN,EMDOT,ECDOT,EKICK,TPLOT,DTPLOT,
     &               XHYD,YHEL,ZMET,ZPARS(20),SPNFAC,IQCOLL,NAS,NBH,
     &               NBKICK,NBR,NBRK,NBS,NCHA,NCIRC,NCOAL,NCONT,NDD,
     &               NEMOD,NGB,NGLOB,NGLOB0,NHE,NHG,NHI,NHYP,NKICK,
     &               NMDOT,NMS,NNH,NRG,NRO,NROCHE,NRS,NRSAVE,NSHOCK,
     &               NSLP,NSN,NSP,NSPIR,INSTAB,NTZ,NWD,NCE,NHYPC,NBH0,
     &               ITYPE(5),KSAVE(4),KTYPE(0:14,0:14),NEINT,IBLUE,
     &               NGDUM(12),LISTR(MLR),LISTD(MLD),LISTV(MLV)


      COMMON/NBODY/  X(3,NMAX),X0(3,NMAX),X0DOT(3,NMAX),F(3,NMAX),
     &               XN(3,NMAX),XNDOT(3,NMAX),
     &               FDOT(3,NMAX),BODY(NMAX),RS(NMAX),XDOT(3,NMAX),
     &               FI(3,NMAX),D1(3,NMAX),D2(3,NMAX),D3(3,NMAX),
     &               FR(3,NMAX),D1R(3,NMAX),D2R(3,NMAX),D3R(3,NMAX),
     &               STEP(NMAX),T0(NMAX),STEPR(NMAX),T0R(NMAX),
     &               TIMENW(NMAX),RADIUS(NMAX),TEV(NMAX),TEV0(NMAX),
     &               BODY0(NMAX),EPOCH(NMAX),SPIN(NMAX),TPRED(NMAX),
     &               ZLMSTY(NMAX),FIDOT(3,NMAX),D0(3,NMAX),
     &               FRDOT(3,NMAX),D0R(3,NMAX),KSTAR(NMAX)

*
      COMMON/PAIRS/  U(4,KMAX),U0(4,KMAX),UDOT(4,KMAX),FU(4,KMAX),
     &               FUDOT(4,KMAX),FUDOT2(4,KMAX),FUDOT3(4,KMAX),
     &               H(KMAX),HDOT(KMAX),HDOT2(KMAX),HDOT3(KMAX),
     &               HDOT4(KMAX),DTAU(KMAX),TDOT2(KMAX),TDOT3(KMAX),
     &               R(KMAX),R0(KMAX),GAMMA(KMAX),SF(7,KMAX),H0(KMAX),
     &               FP0(4,KMAX),FD0(4,KMAX),TBLIST,DTB,KBLIST(10*KMAX),
     &               KSLOW(KMAX),NAME(NMAX),LIST(LMAX,NMAX)
*
      COMMON/LISTS/  ILIST(NMAX),JLIST(NMAX),JPERT(5*LMAX)
*        Common block to keep neighbour density and potential high prec (R.Sp.)
      COMMON/WORK2/RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX),PHII(NMAX),
     &               PHIR(NMAX),PHIR1(NMAX)
*
      COMMON/TIMING/ttota,ttreg,ttirr,ttpre,ttinit,ttint,ttks,
     *  ttcomm,ttadj,ttmov,ttnbp,ttsub,ttsub2,ttfrc,xtsub1,xtsub2,
     *  isernb,iserreg,iserpr

