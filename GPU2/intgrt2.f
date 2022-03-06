      SUBROUTINE INTGRT2
*
*
*       N-body integrator flow control.
*       -------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      PARAMETER (NIMAX=1024,NPMAX=16)
      REAL*8   H2I(NIMAX),XI(3,NIMAX),VI(3,NIMAX),GPUACC(3,NIMAX),
     &         GPUJRK(3,NIMAX),GPUPHI(NIMAX)
      REAL*8  GF(3,2000),GFD(3,2000)
      INTEGER  NXTLST(NMAX),IBL(LMAX),NBLIST(NMAX),LISTQ(NMAX),NL(20)
      INTEGER  IR(NMAX),LISTGP(LMAX,NIMAX)
      LOGICAL LOOP,LSTEPM
      SAVE IQ,ICALL,NQ,LQ,LOOP,LSTEPM,STEPM,ISAVE,JSAVE,ISTART,NNPREV
      DATA IQ,ICALL,LQ,LOOP,LSTEPM,STEPM /0,2,11,.TRUE.,.FALSE.,0.03125/
      DATA ISAVE,JSAVE,ISTART /0,0,0/
      SAVE CPRED,CPRED2,CNB,CNBP
      DATA CPRED,CPRED2,CNB,CNBP /0.0D0,0.0D0,0.0D0,0.0D0/
*
*
*   Open the GPU libraries on each new run (note nnbmax = NN is printed).
      IF (ISTART.EQ.0) THEN
          NN = N
          NNPREV = NN
          CALL GPUNB_OPEN(NN)
          CALL GPUIRR_OPEN(NN,LMAX)
          ISTART = 1
*       Initialize all variables, also include KS components for safety.
          DO 990 I = 1,NTOT
              CALL GPUIRR_SET_JP(I,X(1,I),XDOT(1,I),F(1,I),FDOT(1,I),
     &                                              BODY(I),T0(I))
              CALL GPUIRR_SET_LIST(I,LIST(1,I))
  990     CONTINUE
      END IF
      CALL GPUIRR_PRED_ALL(NTOT,TIME)
      DO 10 I = 1,N
          CALL GPUIRR_FIRR(I,GF(1,I),GFD(1,I))
   10 CONTINUE
      DO 20 I = 1,10
          WRITE (6,15)  FI(1,I), GF(1,I),FIDOT(1,I),GFD(1,I),LIST(1,I),
     &                                                       LIST(2,I)
   15     FORMAT (' FI GP FID GPD   ',1P,4E12.4,0P,2I5)
   20 CONTINUE
      CALL FLUSH(6)
*
*
      STOP
      END
