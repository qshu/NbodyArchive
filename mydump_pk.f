      SUBROUTINE MYDUMP(I,J)
*
*
*       COMMON save or read.
*       --------------------
*
C Note: updtated in Aug.1998 by P.Kroupa to account for additions
C in common6.h
*
      INCLUDE 'params.h'
      IMPLICIT REAL*8  (A-H,O-Z)
      PARAMETER  (NA=(16*ID+6)*NMAX,NB=47*KMAX+2,
     &            NC=(LMAX+2)*NMAX+3*KMAX+MLR+MLD+MLV+73,
     &            NF=40*MMAX,NG=2*NMAX+NMAX/2+533,NH=10*MCL+8,
     &            NO=(4*ID+1)*NMAX)
      INCLUDE 'mpif.h'
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      COMMON/MPIDAT/group,rank,ierr,isize,status
      REAL*8  A,B,E,F,G,H,O
      INTEGER  IC,IDD,IR
*
      COMMON/NBODY/  A(NA)
      COMMON/PAIRS/  B(NB)
      COMMON/NAMES/  IC(NC)
      COMMON/COUNTS/ IDD(90)
      COMMON/PARAMS/ E(152)
      COMMON/BINARY/ F(NF)
C Take note that G=real*8, and some members of it (in common6.h) 
C are integer, is. these only add as half words. 
      COMMON/STARS/  G(NG)
      COMMON/CLOUDS/ H(NH)
      COMMON/RAND/   IR(99)
      COMMON/HERMIT/ O(NO)
*
      if(rank.ne.0)return
*
*       Open unit #J by reading dummy and rewinding.
C     GLOBAL_UNIT = J
      REWIND J
      READ (J,ERR=10,END=10)  DUMMY
   10 REWIND J
*
*       Read or save all COMMON variables (valid for tape or disc).
      IF (I.EQ.0) THEN
          READ (J)  A, B, IC, IDD, E, F, G, H, IR, O
      ELSE
          if(rank.eq.0)then
          WRITE (J) A, B, IC, IDD, E, F, G, H, IR, O
          END FILE J
          CLOSE (UNIT=J)
		  end if
      END IF
*       Write transition data IC(62) = KZ(39)
*     IF(IC(62).EQ.1)THEN
*     IARG = 1
*     CALL WTRANS(IARG)
*     IARG = 2
*     CALL WTRANS(IARG)
*     END IF
*
      RETURN
*
      END
