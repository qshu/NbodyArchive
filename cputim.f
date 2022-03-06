


      SUBROUTINE CPUTIM(TCOMP)
*
*
*       CPU time.
*       ---------
*
      integer icpu
      COMMON/ICPU0/ ICPU
      REAL*4  tt,TARRAY(2),etime
      real * 8 tcomp
*
*
*       Initialize timer (first call) or obtain elapsed time.

      IF (ICPU.EQ.0) THEN
          TT = ETIME(TARRAY)
          TCOMP = tt
          ICPU = 1
      ELSE
*       Elapsed CPU time in minutes on VAX.
          TT = ETIME(TARRAY)/60.0
          TCOMP = tt
*       Elapsed CPU time in minutes on SUN or MIPS.
      END IF
*

      RETURN
*
      END
      subroutine printcpu
      real*8 tcpu
      call cputim(tcpu)
      write(6,600) tcpu
600   format('CPU min = ',f18.6)
      end
      
