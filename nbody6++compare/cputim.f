# 1 "cputim.F"
# 1 "<built-in>"
# 1 "<command line>"
# 1 "cputim.F"
      SUBROUTINE CPUTIM(TCOMP)
*
*
* CPU time.
* ---------
*
      INCLUDE 'common6.h'
      REAL*8 TCOMP
      COMMON/ICPU0/ ICPU
      REAL*4 tt,TARRAY(2),etime
*
* Initialize timer (first call) or obtain elapsed time.
* TCOMP = FLOAT(ITIME)/6000.0
* TCOMP = REAL(IRTC())*6.67D-9/60.
# 29 "cputim.F"
          TCOMP = ETIME(TARRAY)/60.0

* TCOMP = MCLOCK()/6000.
* Elapsed CPU time in minutes on VAX, SUN or MIPS & IBM RS/6000.
* and T3D
*
      RETURN
*
      END
