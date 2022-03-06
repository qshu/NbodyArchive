      SUBROUTINE KSPERI(IPAIR)
*
*
*       Pericentre KS variables.
*       ------------------------
*
      INCLUDE 'common6.h'
*
      CALL KSPERIT(IPAIR,TIME1)
      TIME = TIME1
*
      RETURN
*
      END
