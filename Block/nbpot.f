      SUBROUTINE NBPOT(NB,NP,POTS)
*
*
*       Potential energy of subsystem.
*       ------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Obtain potential energy of subsystem JLIST(NB) & JPERT(NP).
      CALL NBPOT2(NB,NP,POTS,JLIST,JPERT)
*
      RETURN
*
      END
