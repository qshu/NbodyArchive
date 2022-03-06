***
      SUBROUTINE gpinit(i) 
      INTEGER i
      CALL g6_open(i)
      CALL g6_set_tunit(45)
      CALL g6_set_xunit(45)
      RETURN
      END
*
      SUBROUTINE gpwipe(i,time)
      INTEGER i
      REAL*8 time
      CALL g6_reset(i)
      CALL g6_reset_fofpga(i)
      CALL gpfree
      CALL gpinit(i)
      CALL g6_set_ti(i,time)
      RETURN
      END
*
      SUBROUTINE gpfree
      INCLUDE 'common4.h'
      CALL g6_close(gpid)
      RETURN
      END
***
