      SUBROUTINE MYDUMP(I,J)
*
*
*       COMMON save or read.
*       --------------------
*
      INCLUDE 'params.h'
      PARAMETER  (NA=96*NMAX,NB=111*KMAX+7,
     &            NC=(LMAX+2)*NMAX+3*KMAX+MLR+MLD+MLV+64,
     &            NF=41*MMAX,NG=11*NMAX+126,NH=20*MCL+16,
     &            NO=24*NMAX,NP=2*NMAX+84)
      REAL*4  A,B,E,F,G,H,O,P,Q,PL
      INTEGER  IC,ID,IR
*
*
      COMMON/NBODY/  A(NA)
      COMMON/PAIRS/  B(NB)
      COMMON/NAMES/  IC(NC)
      COMMON/COUNTS/ ID(60)
      COMMON/PARAMS/ E(330)
      COMMON/BINARY/ F(NF)
      COMMON/STARS/  G(NG)
      COMMON/CLOUDS/ H(NH)
      COMMON/RAND2/  IR(99)
      COMMON/HERMIT/ O(NO)
      COMMON/BLOCKS/ P(NP)
      COMMON/GALAXY/ Q(40)
      COMMON/PLPOT/  PL(24)
*
*
*       Open unit #J by reading dummy and rewinding.
      REWIND J
      READ (J,ERR=10,END=10)  DUMMY
   10 REWIND J
*
*       Read or save all COMMON variables (valid for tape or disc).
      IF (I.EQ.0) THEN
          READ (J)   A, B, IC, ID, E, F, G, H, IR, O, P, Q, PL
      ELSE
          WRITE (J)  A, B, IC, ID, E, F, G, H, IR, O, P, Q, PL
          END FILE J
          CLOSE (UNIT=J)
      END IF
*
      RETURN
*
      END
