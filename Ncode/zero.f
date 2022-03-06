      SUBROUTINE ZERO
*
*
*       Initialization of global scalars. 
*       ---------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Initialize parameters & counters and set useful constants.
      TIME = 0.0D0
      TOFF = 0.0D0
      TTOT = 0.0D0
      TADJ = 0.0D0
      TBLIST = 0.0D0
      DTB = 0.0D0
      TNEXT = 0.0D0
      TMDOT = 0.0D0
      ZMDOT = 0.0D0
      TLIST = 0.0D0
      TLASTS = 0.0D0
      TLASTT = 0.0D0
      RMAX = 0.0
      EPOCH0 = 0.0
      DMIN1 = 100.0
      DMIN2 = 100.0
      DMIN3 = 100.0
      DMIN4 = 100.0
      DMINC = 100.0
      STAB1 = 1.0
      STAB2 = 2.0
      CPUTOT = 0.0
      ERROR = 0.0D0
      ERRTOT = 0.0D0
      DETOT = 0.0D0
      ETCORR = 0.0
      EBIN0 = 0.0D0
      ESUB = 0.0D0
      EMERGE = 0.0D0
      EDISS = 0.0D0
      ECOLL = 0.0D0
      ECDOT = 0.0D0
      EMDOT = 0.0D0
      EKICK = 0.0
      ESYNC = 0.0D0
      ETIDE = 0.0D0
      SBCOLL = 0.0D0
      BBCOLL = 0.0D0
      CHCOLL = 0.0D0
      NPAIRS = 0
      NBPREV = 0
      NNTB = 0
      NMERGE = 0
      NEWHI = 0
      NSUB = 0
      NCH = 0
      NBIN0 = 0
      NHI0 = 0
      IFIRST = 1
      IPHASE = 0
      ICOMP = 0
      JCOMP = 0
      ICLOSE = 0
      JCLOSE = 0
      KSPAIR = 0
      NPRINT = 0
      MODEL = 0
      NDUMP = 0
      NSTEPS = 0
      NTIMER = 0
      NLIST(1) = 0
      LISTR(1) = 0
      LISTD(1) = 0
      LISTV(1) = 0
      NTAIL = 0
*
      K = 2*KMAX
      DO 10 J = 1,K
          KVEC(J) = 0.50001*FLOAT(J + 1)
   10 CONTINUE
*
      DO 20 J = 1,12
          E(J) = 0.0D0
   20 CONTINUE
*
      DO 30 J = 1,9
          TLASTB(J) = 0.0D0
   30 CONTINUE
*
      DO 40 J = 1,NMAX
          RADIUS(J) = 0.0D0
          BODY0(J) = 0.0
          EPOCH(J) = 0.0
          TEV(J) = 0.0
          TEV0(J) = 0.0
          KSTAR(J) = 0
   40 CONTINUE
*
      DO 50 K = 1,3
          RDENS(K) = 0.0D0
          TIDAL(K) = 0.0
   50 CONTINUE
      TIDAL(4) = 0.0
*
*       Set fractional constants & two PI.
      ONE3 = 1.0D0/3.0D0
      ONE6 = 1.0D0/6.0D0
      ONE9 = 1.0D0/9.0D0
      ONE12 = 1.0D0/12.0D0
      TWOPI = 8.0D0*ATAN(1.0D0)
*
      RETURN
*
      END
