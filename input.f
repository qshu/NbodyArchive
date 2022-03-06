


      SUBROUTINE INPUT
*
*
*       Parameter definition & main input.
*       ----------------------------------
*
      INCLUDE 'common1.h'
*
*
*       Input parameters
*       ****************
*
*       ------------------------------------------------------------------
*       KSTART  Control index (1: new run; >1: restart; 3: new params).
*       TCOMP   Maximum computing time in minutes (saved in CPU).
*
*       N       Total particle number.
*       NRAND   Random number sequence skip.
*
*       ETA     Time-step parameter for total force polynomial.
*       DELTAT  Output time interval in units of the crossing time.
*       TCRIT   Termination time in units of the crossing time.
*       QE      Energy tolerance (stop if DE/E > 5*QE & KZ(2) > 0).
*       CUTOFF  Softening parameter (square saved in EPS2).
*
*       KZ(J)   Non-zero options for alternative paths (see table).
*
*       ALPHAS  Power-law index for initial mass function (routine DATA).
*       BODY1   Maximum particle mass before scaling.
*       BODYN   Minimum particle mass before scaling.
*
*       Q       Virial ratio (routine SCALE; Q = 0.5 for equilibrium).
*       VXROT   XY-velocity scaling factor (> 0 for solid-body rotation).
*       VZROT   Z-velocity scaling factor (not used if VXROT = 0).
*       nbh     number of BHs (point-mass particles with zero eps)
*       ------------------------------------------------------------------
*
*
*       Options KZ(J)
*       *************
*
*       ------------------------------------------------------------------
*       1  Common save on unit 1 if TCOMP > CPU or if TIME > TCRIT.
*       2  Energy error check using tolerance QE (stop if DE/E > 20*QE).
*          2 : Change ETA according to QE
*          4 : restart if DE > QE*20
*       3  Common save on unit 2/3 at output (seems to be unused !!!!!)
*       4  Initial conditions read from unit 4 (BODY, X, Y, Z, Vx, Vy, Vz).
*                  =2 stoa firmat input file
*       5  Initial conditions (=0: uniform & isotropic; =1: Plummer).
*       6  Output of significant binaries.
*       7  Lagrandian shells output 2: unit6, 3: unut 7
*       8  Interval to print bodies (0, 1 : each output)
*       9  Individual bodies printed at output time (MIN(5**KZ9,N)).
*                  =99 stoa format output
*                  >=100 binary output
*                  if nbh > 0, always print out BHs but
*                  select only KZ(9)-100 particles      
*      10  No scaling of output intervals (routine SCALE).
*       ------------------------------------------------------------------
*
*
      integer nrand,j
      real * 8 cutoff
*       Read & print input parameters.
      READ (5,*)  N, nbh, NRAND
      READ (5,*)  ETA, DELTAT, TCRIT, QE, CUTOFF
      READ (5,*)  (KZ(J),J=1,10)
*
      WRITE (6,10)
   10 FORMAT (////,12X,'N  NBH NRAND   ETA   DELTAT   TCRIT    QE',
     &                                         '       CUTOFF')
      WRITE (6,20)  N, NBH,NRAND, ETA, DELTAT, TCRIT, QE, CUTOFF
   20 FORMAT (/,8X,I5,I4,I7,F7.2,F8.1,F8.1,F11.5,F6.2)
      WRITE (6,30)  (KZ(J),J=1,10)
   30 FORMAT (//,12X,'OPTIONS  ',10I4,/)
*
      EPS2 = CUTOFF**2
      NLIST(1) = NRAND
*       Random number sequence skip used by routine DATA.
*
      RETURN
*
      END
