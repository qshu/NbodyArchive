*       COMMON1.
*       --------
*
      INCLUDE 'params.h'
*      implicit REAL*8 (a-h,o-z)
      real *8 x,x0,x0dot,t0,f,fdot,body,xdot,step,phi
*
      COMMON/NBODY/  X(3,NMAX),X0(3,NMAX),X0DOT(3,NMAX),T0(NMAX),
     &               F(3,NMAX),FDOT(3,NMAX),BODY(NMAX),XDOT(3,NMAX),
     &               STEP(NMAX),phi(nmax)
      integer N,NSTEPN,KZ,NLIST,NBLOCK, nbh, nstepbh
*
      COMMON/NAMES/  N,NSTEPN,KZ(10),NLIST(1),NBLOCK,nbh,nstepbh
*
      real *8  TIME,CPU,ETA,DELTAT,TCRIT,QE,EPS2,
     &               ONE3,ONE6,ONE9,ONE12,
     &               TNEXT,TLIST,DTLIST,ZMASS,RSCALE,TCR,BE
      COMMON/PARAMS/ TIME,CPU,ETA,DELTAT,TCRIT,QE,EPS2,
     &               ONE3,ONE6,ONE9,ONE12,
     &               TNEXT,TLIST,DTLIST,ZMASS,RSCALE,TCR,BE(3)
*
      real*8  TPREV,TBLOCK,DTK,TIMENW
      COMMON/BLOCKS/ TPREV,TBLOCK,DTK(64),TIMENW(NMAX)
c      integer zzzzdummy(nmax*20)
c      common/zzzz/zzzzdummy
