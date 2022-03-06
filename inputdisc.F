      SUBROUTINE INPUTDISC
*
*
*       Parameter input.
*       ----------------
*
      INCLUDE 'common6.h'
      REAL*8 R_STAR, EBIND_I(NMAX), R_I(NMAX),INCL,RV2,
     &       ECS_I,BETTA,ETTA,
     &       BETTA1,ETTA1,BETTA2,ETTA2,BETTA3,ETTA3,BETTA4,ETTA4 
      INTEGER INDEXVECTOR(NMAX)
      LOGICAL E_BIND(NMAX)
      COMMON/ENERG/ E_BIND
*     EXTERNAL VERIFY
*
#ifdef PARALLEL
#define MPIINIT 1
#else
#ifdef ENSEMBLE
#define MPIINIT 1
#else
#define MPIINIT 0
#endif
#endif
*
*       Make a formal call to define input parameters & counters.
      CALL DEFINE
*
      DO I = 1,N
      INDEXVECTOR(I)=I
      R_I(I)=SQRT(X(1,I)**2.0+X(2,I)**2.0+X(3,I)**2.0) 
      END DO
      CALL SORT1(N,R_I,INDEXVECTOR)
*
*
      IF(rank.eq.0)THEN
*
      OPEN(15,FILE='drag.dat',STATUS='OLD',ERR=60)
      READ(15,*) Q_DRAG,HZ,SMBH,D_MASS,ALPHAR,D_DELTA
*
      DO I=1,N
      L=INDEXVECTOR(I)
*     IF(I.EQ.100) D_SCALE=R_I(I)
      IF(I.EQ.INT(N*SMBH)) D_SCALE=R_I(I)
      END DO
*
      D_SIGM=5.0*D_MASS/(4*TWOPI*SQRT(TWOPI)*HZ*D_SCALE**2.0)
      R_STAR=SQRT(1/SU**2.0)
      Q_REAL=(2.0*Q_DRAG*SU*SU)/TWOPI
      D_MU=D_MASS/ZMASS      
      BETTA=(0.8*TWOPI*D_SCALE**2.0)/(2.0*D_MU*Q_DRAG*N) 
      ETTA=(5.7*LOG10(0.38*N)*TWOPI*D_SCALE**2.0)/
     &     (2.0*D_MU*Q_DRAG*N**2.0) 
      BETTA1=(0.8*TWOPI*D_SCALE**2.0)/(2.0*D_MU*Q_DRAG*10000) 
      ETTA1=(5.7*LOG10(0.38*10000)*TWOPI*D_SCALE**2.0)/
     &     (2.0*D_MU*Q_DRAG*10000**2.0) 
      BETTA2=(0.8*TWOPI*D_SCALE**2.0)/(2.0*D_MU*Q_DRAG*20000) 
      ETTA2=(5.7*LOG10(0.38*20000)*TWOPI*D_SCALE**2.0)/
     &     (2.0*D_MU*Q_DRAG*20000**2.0) 
      BETTA3=(0.8*TWOPI*D_SCALE**2.0)/(2.0*D_MU*Q_DRAG*30000) 
      ETTA3=(5.7*LOG10(0.38*30000)*TWOPI*D_SCALE**2.0)/
     &     (2.0*D_MU*Q_DRAG*30000**2.0) 
      BETTA4=(0.8*TWOPI*D_SCALE**2.0)/(2.0*D_MU*Q_DRAG*50000) 
      ETTA4=(5.7*LOG10(0.38*50000)*TWOPI*D_SCALE**2.0)/
     &     (2.0*D_MU*Q_DRAG*50000**2.0) 
      WRITE (6,31)
   31 FORMAT(//,3X,'DISK: Q_DRAG    HZ        SMBH      D_SCALE'       
     &      '   D_MASS    ALPHAR      D_DELTA   D_SIGM0')
      WRITE (6,32) Q_DRAG,HZ,SMBH,D_SCALE,D_MASS,ALPHAR,D_DELTA,D_SIGM
   32 FORMAT (/,6X,1P9E10.1)
*
      WRITE (6,33)
   33 FORMAT(//,3X,'DISK SCALING:  Q_REAL    BETTA~    ETTA      D_MU'
     &             '      r*        SU')
      WRITE (6,34) Q_REAL,BETTA,ETTA,D_MU,R_STAR,SU
   34 FORMAT (/,15X,1P6E10.1)
      WRITE (6,35)
   35 FORMAT(//,3X,'BETTA&ETTA:    BETTA10K  ETTA10K'
     &             '   BETTA20K  ETTA20K'
     &             '   BETTA30K  ETTA30K'
     &             '   BETTA50K  ETTA50K')
      WRITE (6,36) BETTA1,ETTA1,BETTA2,ETTA2,BETTA3,ETTA3,BETTA4,ETTA4
   36 FORMAT (/,15X,1P8E10.1)
*     PRINT*, Q_DRAG, HZ, SMBH,  D_SCALE D_MASS, ALPHAR, D_DELTA,D_SIGM
      GO TO 70
   60 PRINT*,'ERROR NO DRAG.DAT FILE!!!'
   70 CONTINUE
*
      END IF
*
#if MPIINIT
*
      CALL MPI_BCAST(Q_DRAG,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(HZ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(SMBH,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(D_SCALE,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(D_MASS,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ALPHAR,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(D_DELTA,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(D_SIGM,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(DN_1,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif
*
      CALL ENERGY
      ETOT = ZKIN - POT + EBLCKHL 
      DO I = 1,N
      INDEXVECTOR(I)=I
      EBIND_I(I)=BODY(I)*(0.5D0*(XDOT(1,I)**2.0+XDOT(2,I)**2.0+ 
     &                           XDOT(3,I)**2.0)+ 
*     &          CMBLHOLE/DSQRT(X(1,I)**2.0+X(2,I)**2.0+X(3,I)**2.0))
     &           PHIDBL(I))
*     PHIDBL(I)-  comprises also EBLCKHL in energy.f
      R_I(I)=SQRT(X(1,I)**2.0+X(2,I)**2.0+X(3,I)**2.0) 
      END DO
*     CALL SORT1(N,R_I,INDEXVECTOR)
*2008     CALL SORT1(N,EBIND_I,INDEXVECTOR)
*
      DO I=1,N
      L=INDEXVECTOR(I)
*
      RV2=(X(2,L)*XDOT(3,L)-X(3,L)*XDOT(2,L))**2.0+
     &    (X(3,L)*XDOT(1,L)-X(1,L)*XDOT(3,L))**2.0+
     &    (X(1,L)*XDOT(2,L)-X(2,L)*XDOT(1,L))**2.0
      ECS_I=SQRT(1.0+2.0*EBIND_I(L)*RV2*BODY(L)/
     &         (SMBH+I*BODY(L))**2.0)
      ANGLE1= (X(1,L)*XDOT(2,L)-X(2,L)*XDOT(1,L))/
     &         SQRT(RV2)
      INCL=(360/TWOPI)*ACOS((X(1,L)*XDOT(2,L)-X(2,L)*XDOT(1,L))/
     &         SQRT(RV2))
*
      if(rank.eq.0) THEN
*     PRINT*,'L;I,E,R, ECS=',L,I,EBIND_I(I),R_I(L),INCL, 
*     &                           ECS_I,RV2, I*BODY(I)
      END IF
      END DO
*
      RETURN
*
      END
