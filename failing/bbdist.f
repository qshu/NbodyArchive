      SUBROUTINE BBDIST
*
*
*       Diagnostics for binary domain decomposition
*       -------------------------------------------
*       plus estimation of parallel speed-up (R.Sp.)
*
      INCLUDE 'common6.h'
      PARAMETER (IBINT=10,NBMAX=10)
      PARAMETER (KNMAX=5000,NDIM=3,NSCELL=2**NDIM)
      COMMON/TREEDAT/PLL(NDIM,0:KNMAX),XNN(0:KNMAX),NAMET(0:KNMAX),
     * XMM(0:KNMAX),XMR(NDIM,0:KNMAX),XMV(NDIM,0:KNMAX),IPAR(0:KNMAX),
     * ILEAVE(0:KNMAX),ISUB(0:KNMAX),THETA2,SMAX
      DIMENSION IPRIM(maxpe)
      DATA IPRIM/2,3,5,7,11,13,17,19,23,29,31,37,(maxpe-12)*0/
*
*       Reset group membership index for all particles
      DO 50 I = 1,NTOT
 50   NGRMEM(I) = 0
*
      DO 100 IPAIR = 1,NPAIRS
*
      I1 = 2*IPAIR - 1
      ICM = N + IPAIR
*
      IF(LIST(1,I1).EQ.0)THEN
*       Unperturbed Binaries
      WRITE(87,888)BODY(ICM),(X(K,ICM),K=1,3),(XDOT(K,ICM),K=1,3),
     &         RHO(ICM),XNDBL(ICM),PHIDBL(ICM),NAME(ICM)
*
      ELSE
      WRITE(88,888)BODY(ICM),(X(K,ICM),K=1,3),(XDOT(K,ICM),K=1,3),
     &         RHO(ICM),XNDBL(ICM),PHIDBL(ICM),NAME(ICM)
      END IF
*
 100  CONTINUE
*
      CLOSE(87)
      CLOSE(88)
*
 888  FORMAT(1X,1P,10(1X,D12.5),I6)
*
      NGR1 = 0
      NGR2 = 0
      DO 1111 IPAIR = 1,NPAIRS
      I1 = 2*IPAIR - 1
      ICM = IPAIR + N
      RI2 = X(I1,ICM)**2 + X(2,ICM)**2 + X(3,ICM)**2
      IF(RI2.GT.4.D0)THEN
      NGR1 = NGR1 + 1
      NGRMEM(I1) = 1
      NGRMEM(ICM) = IPRIM(1)
      LISTB(NGR1+1,1) = ICM
      ELSE
      NGR2 = NGR2 + 1
      NGRMEM(I1) = 2
      NGRMEM(ICM) = IPRIM(2)
      LISTB(NGR2+1,2) = ICM
      END IF
 1111 CONTINUE
*
      LISTB(1,1) = NGR1
      LISTB(1,2) = NGR2
      NGROUP = 2
*
      GOTO 2222
*
      PRINT*,' BBDIST: Going to call TREE'
      CALL FLUSH(6)
*
      IC = 0
      SMAX = 64.D0
      call cputim(tt1)
      CALL TREE(IC)
      call cputim(tt2)
      tttt = (tt2-tt1)*60.
*
      PRINT*,' BBDIST: Returning from TREE with IC =',IC,' tt=',tttt
      CALL FLUSH(6)
*
*      Descend the TREE to collect binary groups
*
*        Initialize parent pointer, node pointer, level pointer
      IN=0
      ILEVEL=0
      NGROUP = 0
      IMEM = 0
      MMEM = 0
*
 1000 CONTINUE
*
*     PRINT*,' ILEVEL=',ILEVEL,' Examining node IN ',IN,' Par=',
*    *   IPAR(IN),' first sub=',ISUB(IN)
*     PRINT*,' Node members=',XNN(IN),' Parent ',XNN(IPAR(IN))
*
*      Members in cell less than desired group membership?
      IF(IMEM.EQ.MMEM.AND.IPAR(IN).NE.IPAR(IN-1))THEN
*     IF(IPAR(IN).NE.IPAR(IN-1).AND.IMEM.EQ.MMEM.AND.
*    *   XNN(IN).LE.NBMAX.AND.XNN(IPAR(IN)).GT.NBMAX)THEN
      NGROUP = NGROUP + 1
      IMEM = 0
      MMEM = XNN(IPAR(IN))
      LISTB(1,NGROUP) = XNN(IPAR(IN))
*     PRINT*,' New Group opened at level ',ILEVEL,
*    *   ' with membership ',XNN(IPAR(IN))
      END IF
*
      IF(ILEAVE(IN).GT.0)THEN
      IMEM = IMEM + 1
*     PRINT*,' ILEVEL=',ILEVEL,' Found leave part=',NAMET(IN),
*    *  ' put to group ',NGROUP,' IMEM=',IMEM
      JPAIR = NAMET(IN) - N
      J1 = 2*JPAIR - 1
      LISTB(IMEM+1,NGROUP) = N + JPAIR
      NGRMEM(N+JPAIR) = IPRIM(NGROUP)
      NGRMEM(J1) = NGROUP
*
*       Look for next node
*
 500  CONTINUE
*
      IF(IPAR(IN+1).EQ.IPAR(IN))THEN
      IN = IN + 1
      GOTO 1000
      END IF
*
*       Goto parent node
      IN = IPAR(IN)
*
      ILEVEL = ILEVEL - 1
      IF(ILEVEL.EQ.0)GOTO 999
*
      GOTO 500
*
      ELSE
*       Goto first subcell
*     PRINT*,' Node ',IN,' further subdivided goto ',ISUB(IN)
      IN = ISUB(IN)
      ILEVEL = ILEVEL + 1
      GOTO 1000
*
      END IF
*
 999  CONTINUE
*
 2222 CONTINUE
*
       DO 700 IGR = 1,NGROUP
*
       IEXT = 0
       NMEM = LISTB(1,IGR)
       NTEMP = NMEM
       DO 750 IMEM = 2,NMEM+1
       I = LISTB(IMEM,IGR)
*
*       Start browsing through perturber list
       IPAIR = I - N
       I1 = 2*IPAIR - 1
       NPERT = LIST(1,I1)
*      
       IF(NPERT.GT.0)THEN
       DO 751 IL = 2,NPERT+1
       JP = LIST(IL,I1)
*       Add perturbers to group if they are from other groups
       IF(MOD(NGRMEM(JP),IPRIM(IGR)).NE.0)THEN
       NTEMP = NTEMP + 1
       LISTB(NTEMP+1,IGR) = JP
       NGRMEM(JP)=NGRMEM(JP)*IPRIM(IGR)
       IEXT = IEXT + 1
       END IF
*
 751   CONTINUE
       END IF
*
       WRITE(91,899)IGR,IMEM,NGRMEM(I1),NGRMEM(I),I,NAME(I),BODY(I),
     *              (X(K,I),K=1,3),(XDOT(K,I),K=1,3)
 750   CONTINUE
*
       PRINT*,' Group ',IGR,' ext perturbers=',IEXT,
     *  'old NMEM=',NMEM,' new NMEM=',NTEMP
*
 700   CONTINUE
*
 899   FORMAT(1X,6I6,7(1X,D12.5))
*
       IBTOT = 0
       DO 555 IGR = 1,NGROUP
 555   IBTOT = IBTOT + LISTB(1,IGR)
*
       PRINT*,' Total Number of Binaries in Groups=',IBTOT,
     *   ' NPAIRS=',NPAIRS
*
*      Prepare final c.m. and c.m.v. of cells
       DO 600 ICC=0,IC
       IF(ILEAVE(ICC).EQ.0)THEN
       WRITE(89,889)XNN(ICC),XMM(ICC),(XMR(K,ICC),K=1,3),
     &  (XMV(K,ICC),K=1,3)
       ELSE
       WRITE(90,889)XNN(ICC),XMM(ICC),(XMR(K,ICC),K=1,3),
     &  (XMV(K,ICC),K=1,3)
       END IF
 600   CONTINUE
*
      CLOSE(89)
      CLOSE(90)
      CLOSE(91)
*
 889   FORMAT(1X,1P,8(1X,D12.5))
*
      CALL FLUSH(6)
       RETURN
*
      END

      


