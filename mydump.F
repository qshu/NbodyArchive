      SUBROUTINE MYDUMP(I,J)
*
*
*       COMMON save or read.
*       --------------------
*
      INCLUDE 'params.h'
      IMPLICIT REAL*8  (A-H,O-Z)
      PARAMETER  (NA=(16*ID+6)*NMAX,NB=65*KMAX+4,
     &            NC=(LMAX+3)*NMAX+3*KMAX+MLR+MLD+MLV+75,ND=71,
     &            NE=163,NF=21*MMAX,NG=5*NMAX+NMAX/2+63,NH=10*MCL+8,
     &            NO=(4*ID+1)*NMAX,NP=66,NQ=2*NMAX+5*LMAX,
     &            NR=3*NMAX,NS=17)
      INCLUDE 'mpif.h'
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      COMMON/MPIDAT/group,rank,ierr,isize,status
      REAL*8  A,B,E,F,G,H,O,P,R,S
      INTEGER  IC,IDD,IR,IQ
*
      COMMON/NBODY/  A(NA)
      COMMON/PAIRS/  B(NB)
      COMMON/NAMES/  IC(NC)
      COMMON/COUNTS/ IDD(ND)
      COMMON/PARAMS/ E(NE)
      COMMON/BINARY/ F(NF)
C Take note that G=real*8, and some members of it (in common6.h) 
C are integer, is. these only add as half words. 
      COMMON/STARS/  G(NG)
      COMMON/CLOUDS/ H(NH)
      COMMON/RAND/   IR(99)
      COMMON/HERMIT/ O(NO)
      COMMON/BLOCKS/ P(NP)
      COMMON/LISTS/  IQ(NQ)
      COMMON/WORK2/  R(NR)
      COMMON/TIMING/ S(NS)
*
*     Open unit #J by reading dummy and rewinding.
C     GLOBAL_UNIT = J
#ifdef PARALLEL
      if(rank.eq.0)then
#endif
         REWIND J
         READ (J,ERR=10,END=10)  DUMMY
 10      REWIND J
#ifdef PARALLEL
      end if
#endif
*
*       Read or save all COMMON variables (valid for tape or disc).
      IF (I.EQ.0) THEN
*
#ifdef PARALLEL
            if(rank.eq.0)then
#endif
            READ (J)  A, B, IC, IDD, E, F, G, H, IR, O, P, IQ, R, S
#ifdef PARALLEL
            end if
*
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
*
         CALL MPI_BCAST(A(1),NA,MPI_REAL,0,MPI_COMM_WORLD,ierr) 
         CALL MPI_BCAST(B(1),NB,MPI_REAL,0,MPI_COMM_WORLD,ierr) 
         CALL MPI_BCAST(IC(1),NC,MPI_INTEGER,0,MPI_COMM_WORLD,
     &        ierr) 
         CALL MPI_BCAST(IDD(1),ND,MPI_INTEGER,0,MPI_COMM_WORLD,
     &        ierr) 
         CALL MPI_BCAST(E(1),NE,MPI_REAL,0,MPI_COMM_WORLD,ierr) 
         CALL MPI_BCAST(F(1),NF,MPI_REAL,0,MPI_COMM_WORLD,ierr) 
         CALL MPI_BCAST(G(1),NG,MPI_REAL,0,MPI_COMM_WORLD,ierr) 
         CALL MPI_BCAST(H(1),NH,MPI_REAL,0,MPI_COMM_WORLD,ierr) 
         CALL MPI_BCAST(IR(1),99,MPI_INTEGER,0,MPI_COMM_WORLD,
     &        ierr)
         CALL MPI_BCAST(O(1),NO,MPI_REAL,0,MPI_COMM_WORLD,ierr) 
         CALL MPI_BCAST(P(1),NP,MPI_REAL,0,MPI_COMM_WORLD,ierr) 
         CALL MPI_BCAST(IQ(1),NQ,MPI_INTEGER,0,MPI_COMM_WORLD,
     &        ierr)
         CALL MPI_BCAST(R(1),NR,MPI_REAL,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(S(1),NS,MPI_REAL,0,MPI_COMM_WORLD,ierr)
*
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif         
*
      ELSE
*
#ifdef PARALLEL
          if(rank.eq.0)then
#endif
          WRITE (J) A, B, IC, IDD, E, F, G, H, IR, O, P, IQ, R, S
          END FILE J
          CLOSE (UNIT=J)
#ifdef PARALLEL
          end if
#endif
*
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
