# 1 "nbody6.F"
	PROGRAM NBODY6
*
*             N B O D Y 6++
*             *************
*
*       Regularized AC N-body code with triple & binary collisions.
*       --------------------------------------------------------
*
*       Hermite integration scheme with block-steps (V 4.0.0 April/99).
*       ------------------------------------------------------------------
*
*       Developed by Sverre Aarseth, IOA, Cambridge.
*       ............................................
*       Message Passing Version NBODY6++ for Massively Parallel Systems
*       Developed by Rainer Spurzem, ARI, Heidelberg
*
      
# 1 "./common6.h" 1 
*       COMMON6.
*       -------
*
      
# 1 "./params.h" 1 
*       NBODY6 parameters.
*       ------------------
*
*     PARAMETER  (NMAX=50000,KMAX=50,LMAX=64,MMAX=10,
*     PARAMETER  (NMAX=32768,KMAX=50,LMAX=64,MMAX=10,
*     PARAMETER  (NMAX=16384,KMAX=50,LMAX=64,MMAX=10,
      PARAMETER  (NMAX=10100,KMAX=3000,LMAX=256,MMAX=60,
*     PARAMETER  (NMAX=5100,KMAX=2000,LMAX=256,MMAX=60,
*     PARAMETER  (NMAX=1010,KMAX=50,LMAX=60,MMAX=1,
     *            MLD=60,MLR=60,MLV=60,MCL=1,NCMAX=10)
      PARAMETER  (ID=3)
      parameter (maxpe=1024)
*
*
*       ------------------------------------------------------
*       NMAX    Maximum number of single bodies & c.m.
*       KMAX    Maximum number of KS solutions.
*       LMAX    Maximum size of neighbour lists.
*       MMAX    Maximum number of merged binaries.
*       MLD     Maximum number of disrupted KS components.
*       MLR     Maximum number of recently regularized KS pairs.
*       MLV     Maximum number of high-velocity particles.
*       MCL     Maximum number of interstellar clouds.
*       NCMAX   Maximum members of chain members.
*       ------------------------------------------------------
*       ID      First dimension of 3d-vectors (usually = 3)
*               Has to be set to 4 on CRAY T3D for hardware reasons
# 5 "./common6.h" 2 
      IMPLICIT REAL*8  (A-H,O-Z)
      INTEGER BK
*
      
# 1 "./mpif.h" 1 
!
!  
!  (C) 1993 by Argonne National Laboratory and Mississipi State University.
!      All rights reserved.  See COPYRIGHT in top-level directory.
!
!
! user include file for MPI programs, with no dependencies 
!
! It really is not possible to make a perfect include file that can
! be used by both F77 and F90 compilers, but this is close.  We have removed
! continuation lines (allows free form input in F90); systems whose
! Fortran compilers support ! instead of just C or * for comments can
! globally replace a C in the first column with !; the resulting file
! should work for both Fortran 77 and Fortran 90.
!
! If your Fortran compiler supports ! for comments, you can run this 
! through sed with
!     sed -e 's/^C/\!/g'
!
! We have also removed the use of contractions (involving the single quote)
! character because some users use .F instead of .f files (to invoke the
! cpp preprocessor) and further, their preprocessor is determined to find
! matching single quote pairs (and probably double quotes; given the
! different rules in C and Fortran, this sounds like a disaster).  Rather than
! take the position that the poor users should get a better system, we
! have removed the text that caused problems.  Of course, the users SHOULD
! get a better system...
!
! return codes 
      INTEGER MPI_SUCCESS,MPI_ERR_BUFFER,MPI_ERR_COUNT,MPI_ERR_TYPE
      INTEGER MPI_ERR_TAG,MPI_ERR_COMM,MPI_ERR_RANK,MPI_ERR_ROOT
      INTEGER MPI_ERR_GROUP
      INTEGER MPI_ERR_OP,MPI_ERR_TOPOLOGY,MPI_ERR_DIMS,MPI_ERR_ARG
      INTEGER MPI_ERR_UNKNOWN,MPI_ERR_TRUNCATE,MPI_ERR_OTHER
      INTEGER MPI_ERR_INTERN,MPI_ERR_IN_STATUS,MPI_ERR_PENDING
      INTEGER MPI_ERR_REQUEST, MPI_ERR_LASTCODE
      PARAMETER (MPI_SUCCESS=0,MPI_ERR_BUFFER=1,MPI_ERR_COUNT=2)
      PARAMETER (MPI_ERR_TYPE=3,MPI_ERR_TAG=4,MPI_ERR_COMM=5)
      PARAMETER (MPI_ERR_RANK=6,MPI_ERR_ROOT=7,MPI_ERR_GROUP=8)
      PARAMETER (MPI_ERR_OP=9,MPI_ERR_TOPOLOGY=10,MPI_ERR_DIMS=11)
      PARAMETER (MPI_ERR_ARG=12,MPI_ERR_UNKNOWN=13)
      PARAMETER (MPI_ERR_TRUNCATE=14,MPI_ERR_OTHER=15)
      PARAMETER (MPI_ERR_INTERN=16,MPI_ERR_IN_STATUS=17)
      PARAMETER (MPI_ERR_PENDING=18,MPI_ERR_REQUEST=19)
      PARAMETER (MPI_ERR_LASTCODE=1073741823)
!
      INTEGER MPI_UNDEFINED
      parameter (MPI_UNDEFINED = (-32766))
!
      INTEGER MPI_GRAPH, MPI_CART
      PARAMETER (MPI_GRAPH = 1, MPI_CART = 2)
      INTEGER  MPI_PROC_NULL
      PARAMETER ( MPI_PROC_NULL = (-1) )
!
      INTEGER MPI_BSEND_OVERHEAD
      PARAMETER ( MPI_BSEND_OVERHEAD = 512 )

      INTEGER MPI_SOURCE, MPI_TAG, MPI_ERROR
      PARAMETER(MPI_SOURCE=2, MPI_TAG=3, MPI_ERROR=4)
      INTEGER MPI_STATUS_SIZE
      PARAMETER (MPI_STATUS_SIZE=4)
      INTEGER MPI_MAX_PROCESSOR_NAME, MPI_MAX_ERROR_STRING
      PARAMETER (MPI_MAX_PROCESSOR_NAME=256)
      PARAMETER (MPI_MAX_ERROR_STRING=512)
      INTEGER MPI_MAX_NAME_STRING
      PARAMETER (MPI_MAX_NAME_STRING=63)
!
      INTEGER MPI_COMM_NULL
      PARAMETER (MPI_COMM_NULL=0)
!
      INTEGER MPI_DATATYPE_NULL
      PARAMETER (MPI_DATATYPE_NULL = 0)
      
      INTEGER MPI_ERRHANDLER_NULL
      PARAMETER (MPI_ERRHANDLER_NULL = 0)
      
      INTEGER MPI_GROUP_NULL
      PARAMETER (MPI_GROUP_NULL = 0)
      
      INTEGER MPI_KEYVAL_INVALID
      PARAMETER (MPI_KEYVAL_INVALID = 0)
      
      INTEGER MPI_REQUEST_NULL
      PARAMETER (MPI_REQUEST_NULL = 0)
! 
      INTEGER MPI_IDENT, MPI_CONGRUENT, MPI_SIMILAR, MPI_UNEQUAL
      PARAMETER (MPI_IDENT=0, MPI_CONGRUENT=1, MPI_SIMILAR=2)
      PARAMETER (MPI_UNEQUAL=3)
!
!     MPI_BOTTOM needs to be a known address; here we put it at the
!     beginning of the common block.  The point-to-point and collective
!     routines know about MPI_BOTTOM, but MPI_TYPE_STRUCT as yet does not.
!
!     MPI_STATUS_IGNORE and MPI_STATUSES_IGNORE are similar objects
!     Until the underlying MPI library implements the C version of these
!     (a null pointer), these are declared as arrays of MPI_STATUS_SIZE
!
!     The types MPI_INTEGER1,2,4 and MPI_REAL4,8 are OPTIONAL.
!     Their values are zero if they are not available.  Note that
!     using these reduces the portability of code (though may enhance
!     portability between Crays and other systems)
!
      INTEGER MPI_TAG_UB, MPI_HOST, MPI_IO
      INTEGER MPI_BOTTOM
      INTEGER MPI_STATUS_IGNORE(MPI_STATUS_SIZE)
      INTEGER MPI_STATUSES_IGNORE(MPI_STATUS_SIZE)
      INTEGER MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION 
      INTEGER MPI_COMPLEX, MPI_DOUBLE_COMPLEX,MPI_LOGICAL
      INTEGER MPI_CHARACTER, MPI_BYTE, MPI_2INTEGER, MPI_2REAL
      INTEGER MPI_2DOUBLE_PRECISION, MPI_2COMPLEX, MPI_2DOUBLE_COMPLEX
      INTEGER MPI_UB, MPI_LB
      INTEGER MPI_PACKED, MPI_WTIME_IS_GLOBAL
      INTEGER MPI_COMM_WORLD, MPI_COMM_SELF, MPI_GROUP_EMPTY
      INTEGER MPI_SUM, MPI_MAX, MPI_MIN, MPI_PROD, MPI_LAND, MPI_BAND
      INTEGER MPI_LOR, MPI_BOR, MPI_LXOR, MPI_BXOR, MPI_MINLOC
      INTEGER MPI_MAXLOC
      INTEGER MPI_OP_NULL
      INTEGER MPI_ERRORS_ARE_FATAL, MPI_ERRORS_RETURN
!
      PARAMETER (MPI_ERRORS_ARE_FATAL=119)
      PARAMETER (MPI_ERRORS_RETURN=120)
!
      PARAMETER (MPI_COMPLEX=23,MPI_DOUBLE_COMPLEX=24,MPI_LOGICAL=25)
C     PARAMETER (MPI_REAL=26,MPI_DOUBLE_PRECISION=27,MPI_INTEGER=28)
      PARAMETER (MPI_REAL=11,MPI_DOUBLE_PRECISION=27,MPI_INTEGER=28)
      PARAMETER (MPI_2INTEGER=29,MPI_2COMPLEX=30,MPI_2DOUBLE_COMPLEX=31)
      PARAMETER (MPI_2REAL=32,MPI_2DOUBLE_PRECISION=33,MPI_CHARACTER=1)
      PARAMETER (MPI_BYTE=3,MPI_UB=16,MPI_LB=15,MPI_PACKED=14)

      INTEGER MPI_ORDER_C, MPI_ORDER_FORTRAN 
      PARAMETER (MPI_ORDER_C=56, MPI_ORDER_FORTRAN=57)
      INTEGER MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_CYCLIC
      INTEGER MPI_DISTRIBUTE_NONE, MPI_DISTRIBUTE_DFLT_DARG
      PARAMETER (MPI_DISTRIBUTE_BLOCK=121, MPI_DISTRIBUTE_CYCLIC=122)
      PARAMETER (MPI_DISTRIBUTE_NONE=123)
      PARAMETER (MPI_DISTRIBUTE_DFLT_DARG=-49767)
      INTEGER MPI_MAX_INFO_KEY, MPI_MAX_INFO_VAL
      PARAMETER (MPI_MAX_INFO_KEY=255, MPI_MAX_INFO_VAL=1024)
      INTEGER MPI_INFO_NULL
      PARAMETER (MPI_INFO_NULL=0)

!
! Optional Fortran Types.  Configure attempts to determine these.  
!
      INTEGER MPI_INTEGER1, MPI_INTEGER2, MPI_INTEGER4, MPI_INTEGER8
      INTEGER MPI_INTEGER16
      INTEGER MPI_REAL4, MPI_REAL8, MPI_REAL16
      INTEGER MPI_COMPLEX8, MPI_COMPLEX16, MPI_COMPLEX32
      PARAMETER (MPI_INTEGER1=1,MPI_INTEGER2=4)
      PARAMETER (MPI_INTEGER4=6)
      PARAMETER (MPI_INTEGER8=13)
      PARAMETER (MPI_INTEGER16=0)
      PARAMETER (MPI_REAL4=10)
      PARAMETER (MPI_REAL8=11)
      PARAMETER (MPI_REAL16=0)
      PARAMETER (MPI_COMPLEX8=23)
      PARAMETER (MPI_COMPLEX16=24)
      PARAMETER (MPI_COMPLEX32=0)

      COMMON /MPIPRIV/ MPI_BOTTOM,MPI_STATUS_IGNORE,MPI_STATUSES_IGNORE 
!
!     Without this save, some Fortran implementations may make the common
!     dynamic!
!    
!     For a Fortran90 module, we might replace /MPIPRIV/ with a simple
!     SAVE MPI_BOTTOM
!
      SAVE /MPIPRIV/

      PARAMETER (MPI_MAX=100,MPI_MIN=101,MPI_SUM=102,MPI_PROD=103)
      PARAMETER (MPI_LAND=104,MPI_BAND=105,MPI_LOR=106,MPI_BOR=107)
      PARAMETER (MPI_LXOR=108,MPI_BXOR=109,MPI_MINLOC=110)
      PARAMETER (MPI_MAXLOC=111, MPI_OP_NULL=0)
!
      PARAMETER (MPI_GROUP_EMPTY=90,MPI_COMM_WORLD=91,MPI_COMM_SELF=92)
      PARAMETER (MPI_TAG_UB=80,MPI_HOST=82,MPI_IO=84)
      PARAMETER (MPI_WTIME_IS_GLOBAL=86)
!
      INTEGER MPI_ANY_SOURCE
      PARAMETER (MPI_ANY_SOURCE = (-2))
      INTEGER MPI_ANY_TAG
      PARAMETER (MPI_ANY_TAG = (-1))
!
      INTEGER MPI_VERSION, MPI_SUBVERSION
      PARAMETER (MPI_VERSION    = 1, MPI_SUBVERSION = 2)
!
!     There are additional MPI-2 constants 
      INTEGER MPI_ADDRESS_KIND, MPI_OFFSET_KIND
      PARAMETER (MPI_ADDRESS_KIND=4)
      PARAMETER (MPI_OFFSET_KIND=8)
!
!     All other MPI routines are subroutines
!     This may cause some Fortran compilers to complain about defined and
!     not used.  Such compilers should be improved.
!
!     Some Fortran compilers will not link programs that contain
!     external statements to routines that are not provided, even if
!     the routine is never called.  Remove PMPI_WTIME and PMPI_WTICK
!     if you have trouble with them.
!
      DOUBLE PRECISION MPI_WTIME, MPI_WTICK,PMPI_WTIME,PMPI_WTICK
      EXTERNAL MPI_WTIME, MPI_WTICK,PMPI_WTIME,PMPI_WTICK
!
!     The attribute copy/delete subroutines are symbols that can be passed
!     to MPI routines
!
      EXTERNAL MPI_NULL_COPY_FN, MPI_NULL_DELETE_FN, MPI_DUP_FN
# 9 "./common6.h" 2 
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      SAVE
      COMMON/MPIDAT/group,rank,ierr,isize,status
*
      COMMON/NBODY/  X(ID,NMAX),X0(ID,NMAX),X0DOT(ID,NMAX),F(ID,NMAX),
     *               XN(ID,NMAX),XNDOT(ID,NMAX),
     *               FDOT(ID,NMAX),BODY(NMAX),RS(NMAX),XDOT(ID,NMAX),
     *               FI(ID,NMAX),D1(ID,NMAX),D2(ID,NMAX),D3(ID,NMAX),
     *               FR(ID,NMAX),D1R(ID,NMAX),D2R(ID,NMAX),D3R(ID,NMAX),
     *               STEP(NMAX),T0(NMAX),STEPR(NMAX),T0R(NMAX)
*
      COMMON/PAIRS/  HT,U(4,KMAX),U0(4,KMAX),UDOT(4,KMAX),FU(4,KMAX),
     *               FUDOT(4,KMAX),FUDOT2(4,KMAX),FUDOT3(4,KMAX),
     *               H(KMAX),HDOT(KMAX),HDOT2(KMAX),HDOT3(KMAX),
     *               HDOT4(KMAX),DTAU(KMAX),TDOT2(KMAX),TDOT3(KMAX),
     *               R(KMAX),R0(KMAX),GAMMA(KMAX),SF(7,KMAX),H0(KMAX),
     *               FP0(4,KMAX),FD0(4,KMAX),TBLIST,DTB,KBLIST(10*KMAX),
     *               NNTB
*
      COMMON/NAMES/  N,NTOT,NPAIRS,NNBOPT,NCRIT,NFIX,NMERGE,NSUB,NCH,
     *               IPHASE,IFIRST,ICOMP,JCOMP,ICLOSE,JCLOSE,KSPAIR,
     *               NRUN,MODEL,NC,NZERO,NBIN0,IDUM1,NAME0,NHI0,KZ(40),
     *               BK(10),LSHORT(NMAX),NNBMAX,
     *               KVEC(2*KMAX),LISTR(MLR),LISTD(MLD),LISTV(MLV),
     *               KSLOW(KMAX),NLIST(NMAX),NAME(NMAX),LIST(LMAX,NMAX)
*
      COMMON/COUNTS/ NSTEPI,NSTEPR,NSTEPU,NNPRED,NBCORR,NBFULL,NBVOID,
     *               NRCONV,NBSMIN,NLSMIN,NBDIS,NBDIS2,NCMDER,NBDER,
     *               NFAST,NBFAST,NBLOCK,NBPRED,NICONV,NCHAIN,NSTEPC,
     *               NKSTRY,NKSREG,NKSHYP,NKSPER,NPRECT,NKSREF,NKSMOD,
     *               NTTRY,NTRIP,NQUAD,NMERG,NSTEPT,NSTEPQ,NDISS,NTIDE,
     *               NCOLL,NSYNC,NSESC,NBESC,NMESC,NTIMER,NSTEPS,NPRINT,
     *               NDUMP,NBPREV,NPOP(10),JCMAX,NEWHI,NBLCKR,IPE,NPES,
     *               NSTEPB,NBFLUX,NDUMMY(8)
*
      COMMON/PARAMS/ CPU,ETAI,ETAR,DTADJ,DELTAT,TCRIT,TCRITp,
     *               QE,RBAR,ZMBAR,
     *               DTMIN,RMIN,ETAU,CMSEP2,ECLOSE,GMIN,GMAX,ETA0,
     *               TWOPI,ONE3,ONE6,ONE9,ONE12,TCR0,TRH,BODYM,BODY1,
     *               SMIN,RMIN2,RMIN22,FCRIT2,ALPHA,ZNBMIN,ZNBMAX,EBH,
     *               TADJ,TNEXT,CPU0,CPUTOT,TLIST,DTLIST,STAB1,STAB2,
     *               TIME,ZMASS,RSCALE,TCR,TRC,BE(3),CMR(4),CMRDOT(4),
     *               ZKIN,POT,EBIN,EBIN0,ESUB,EMERGE,ECOLL,EDISS,ESYNC,
     *               E(12),ERROR,ERRTOT,DETOT,ETCORR,AZ,PCRIT,EBCH0,
     *               RTIDE,TSCALE,TIDAL(4),ETIDE,RSFAC,RSPH2,BETA,
     *               RC,RC2,RC2IN,VC,ZMC,RDENS(3),RHOD,RHOM,RSMIN,RMAX,
     *               DMIN1,DMIN2,DMIN3,DMIN4,DMINC,SBCOLL,BBCOLL,CHCOLL,
     *               DELTAS,ORBITS(9),GPRINT(9),TLASTT,TLASTS,TLASTB(9),
     *               TDUMP,SCOEFF(12),TOFF,TTOT,DUMMY
*
      COMMON/STARS/  EPOCH0,ZMRG,ZMHE,ZMRS,ZMWD,ZMSN,ZMNH,ZMBH,ZMDOT,
     *               AU,PC,GM,DAYS,YRS,SU,SMU,RAU,TSTAR,VSTAR,STEPX,
     *               TMDOT,TPHYS,TURN,EMDOT,ECDOT,EKICK,TPLOT,DUMMY1,
     *               XHYD,YHEL,ZMET,ZPARS(20),
     *               RADIUS(NMAX),TEV(NMAX),TEV0(NMAX),BODY0(NMAX),
     *               EPOCH(NMAX),NMDOT,NRG,NHE,NRS,NWD,NNH,NSN,NBS,NTZ,
     *               NBH,NKICK,NBKICK,KSAVE(2),ITYPE(10),KSTAR(NMAX)
*
      COMMON/HERMIT/ FIDOT(ID,NMAX),D0(ID,NMAX),FRDOT(ID,NMAX),
     *               D0R(ID,NMAX),TIMENW(NMAX)
*
      COMMON/BLOCKS/ TPREV,TBLOCK,DTK(64)
*
      COMMON/LISTS/  ILIST(NMAX),JLIST(NMAX),JPERT(5*LMAX)
*        Common block to keep neighbour density and potential high prec (R.Sp.)
      COMMON/WORK2/RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
*
      COMMON/TIMING/ttota,ttreg,ttirr,ttpre,ttinit,ttint,ttks,
     *  ttcomm,ttadj,ttmov,ttnbp,ttsub,ttsub2,ttfrc,
     *  isernb,iserreg











# 18 "nbody6.F" 2 
      COMMON/STSTAT/  TINIT,NIR,NIB,NRGL,NKS
      EXTERNAL MERGE
*
# 23

# 26





# 38

*
*       Initialize the timer.
      CALL CPUTIM(ttota)
*
*       Read start/restart indicator & CPU time.
      print*,' rank=',rank
      IF(rank.eq.0)READ (5,*)  KSTART, TCOMP, TCRITp,
     *    isernb,iserreg
*
# 58

*
      IF (KSTART.EQ.1) THEN
*
*       Read input parameters, perform initial setup and obtain output.
          CPU = TCOMP
          CALL START
          call cputim(tt7)
          CALL ADJUST
          call cputim(tt8)
          ttadj = ttadj + (tt8-tt7)*60.
      ELSE
*
*       Open unit 1 only.
          CALL FILE_INIT(1)
*       Read previously saved COMMON variables from tape/disc on unit 1.
          CALL MYDUMP(0,1)
*
*       Open all other files.
          CALL FILE_INIT(2)
*
          IF (NDUMP.GE.3) STOP
*       Safety indicator preventing repeated restarts set in routine CHECK.
          CPU = TCOMP
          CPU0 = 0.0 
*       Set IPHASE = -1 for new NLIST in routine INTGRT (Hermite version).
          IPHASE = -1
*
*       Initialize evolution parameters which depend on metallicity.
          CALL ZCNSTS(ZMET,ZPARS)
*
*       Check reading modified restart parameters (KSTART = 3, 4 or 5).
          IF (KSTART.GT.2) THEN
              CALL MODIFY(KSTART)
          END IF
      END IF
*
* (R.Sp.)Set time flag and step number flags for beginning of run
      TINIT = TTOT
      NIR = NSTEPI
      NIB = NSTEPB
      NRGL = NSTEPR
      NKS = NSTEPU
*
      call cputim(tt2)
      ttinit = ttinit + (tt2-ttota)*60.
*       Advance solutions until next output or change of procedure.
    1 CONTINUE
      call cputim(tt1)
*
      CALL INTGRT
*
      call cputim(tt2)
      ttint = ttint + (tt2-tt1)*60.
*
      IF (IPHASE.EQ.1) THEN
*       Prepare new KS regularization.
      call cputim(tt1)
          CALL KSREG
          CALL FLUSH(6)
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.2) THEN
*       Terminate KS regularization.
      call cputim(tt1)
          CALL KSTERM
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.3) THEN
*       Perform energy check & parameter adjustments and print diagnostics.
          call cputim(tt7)
          CALL ADJUST
          call cputim(tt8)
          ttadj = ttadj + (tt8-tt7)*60.
*
      ELSE IF (IPHASE.EQ.4) THEN
*       Switch to unperturbed three-body regularization.
      call cputim(tt1)
          ISUB = 0 
          CALL TRIPLE(ISUB)
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.5) THEN
*       Switch to unperturbed four-body regularization.
      call cputim(tt1)
          ISUB = 0
          CALL QUAD(ISUB)
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
*       Adopt c.m. approximation for inner binary in hierarchical triple.
      ELSE IF (IPHASE.EQ.6) THEN
      call cputim(tt1)
          CALL MERGE
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.7) THEN
*       Restore old binary in hierarchical configuration.
      call cputim(tt1)
          CALL RESET
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
*       Begin chain regularization.
      ELSE IF (IPHASE.EQ.8) THEN
      call cputim(tt1)
          ISUB = 0
          CALL CHAIN(ISUB)
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
      END IF
*
*       Continue integration.
      GO TO 1
*
      END
