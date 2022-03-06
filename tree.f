      SUBROUTINE TREE(IC)
*
      INCLUDE 'common6.h'
*
      PARAMETER (KNMAX=5000,NDIM=3,NSCELL=2**NDIM)
      COMMON/TREEDAT/PLL(NDIM,0:KNMAX),XNN(0:KNMAX),NAMET(0:KNMAX),
     * XMM(0:KNMAX),XMR(NDIM,0:KNMAX),XMV(NDIM,0:KNMAX),IPAR(0:KNMAX),
     * ILEAVE(0:KNMAX),ISUB(0:KNMAX),THETA2,SMAX
*
      LOGICAL L,LOCINS
      DIMENSION PMIN(NDIM),PMAX(NDIM),XPERM(NDIM,NSCELL),L(NDIM)
      SAVE ITREE
      DATA XPERM/0,0,0,
     *           1,0,0,
     *           0,1,0,
     *           0,0,1,
     *           1,1,0,
     *           1,0,1,
     *           0,1,1,
     *           1,1,1/
*              
*        Initialize parent pointer, node pointer, level pointer
      IN=0
      IC=0
      ILEAVE(IC)=0
      IPAR(0)=0
      ILEVEL=0
      ITREE = ITREE + 1
*        Initialize coordinates of lower left corner of root cell
*        Initialize root node mass and velocity
      XNN(IN)=0.D0
      XMM(IN)=0.D0
      DO 20 K=1,NDIM
      XMR(K,IC)=0.D0
      XMV(K,IC)=0.D0
 20   PLL(K,IC)=-SMAX/2.D0
*        Initialize cell size
      S=SMAX
*        Main loop
 1000 CONTINUE
*
*      First subcell index
      ISCELL=0
      ISUB(IN)=IC+1
*
*     PRINT*,' Level,Node,Cell=',ILEVEL,IN,IC
*     PRINT*,' Parent and 1st subcell=',
*    *   IPAR(IN),ISUB(IN)
*        Reduce cell size
      S=SMAX/2.D0**(ILEVEL+1)
*        Loop for subcells
 500  CONTINUE
*
      ISCELL=ISCELL+1
*
*      Determine minimum and maximum coordinates of cell
      DO 50 K=1,NDIM
      PMIN(K)=PLL(K,IN)+XPERM(K,ISCELL)*S
 50   PMAX(K)=PMIN(K)+S
*
*     PRINT*,' pmin=',(pmin(k),k=1,NDIM)
*     PRINT*,' pmax=',(pmax(k),k=1,NDIM)
*
*        Search for particles in cell
      ICELL = 0
      DO 100 IPAIR=1,NPAIRS
      I = N + IPAIR
      DO 101 K=1,NDIM
 101  L(K) = (X(K,I).GT.PMIN(K).AND.X(K,I).LT.PMAX(K))
      ILOCIN = 0
      DO 102 K=1,NDIM
 102  IF(L(K))ILOCIN=ILOCIN+1
      LOCINS = ILOCIN.EQ.NDIM
      IF(LOCINS)THEN
      ICELL=ICELL+1
*     PRINT*,' Part ',I,' LOCIN=',LOCIN,(L(K),K=1,NDIM)
      IF(ICELL.EQ.1)IPART=I
      END IF
      IF(ICELL.EQ.2)GOTO 105
 100  CONTINUE
*
 105  CONTINUE
*
*      PRINT*,' This cell has ICELL=',ICELL,' parent ',IN
      IF(ICELL.GT.2)ICELL=2
*
*        Are there any particles in cell?
      IF(ICELL.GT.0)THEN
*                                  
          IC = IC + 1
          ILEAVE(IC)=0
          IPAR(IC)=IN
          IPAR(IC+1)=0
*
*        Is there more than one particle?
*        Initialize m,x,v of new node and determine lower left corner
          IF(ICELL.GT.1)THEN
              XNN(IC)=0.D0
              XMM(IC)=0.D0
              DO 200 K=1,NDIM
              PLL(K,IC)=PMIN(K)
              XMV(K,IC)=0.D0
 200          XMR(K,IC)=0.D0
*         PRINT*,' ISCELL=',ISCELL,
*    *     ' This cell ',IC,' is a node pll=',
*    *    (pll(k,ic),k=1,NDIM)
          ELSE
*        There is only one particle in cell?
              ILEAVE(IC)=IPART
              XNN(IC)=1.D0
              NAMET(IC)=IPART
              XMM(IC)=BODY(IPART)
              DO 300 K=1,NDIM
              XMV(K,IC)=XDOT(K,IPART)
 300          XMR(K,IC)=X(K,IPART)
*         PRINT*,' ISCELL=',ISCELL,
*    *      ' This cell ',IC,' is a leaf part=',IPART
          END IF
*
       ELSE
*         PRINT*,' ISCELL=',ISCELL,' This cell has no particles in it'
       END IF 
*
       IF(ISCELL.LT.NSCELL)GOTO 500
*        Set pointer to first subcell
*      PRINT*,' Now we go from ',IN,' to first subcell ',ISUB(IN)
       IN = ISUB(IN)
*        Increase tree level 
       ILEVEL = ILEVEL + 1
*     PRINT*,' New Level=',ILEVEL
*        Loop over subcells
 2000  CONTINUE
*        Is this node a leave?
*      PRINT*,' Now at cell ',IN,' is it leave 0/1? ',
*    *    ILEAVE(IN)
*        If this node is not a leave go further down
       IF(ILEAVE(IN).EQ.0)GOTO 1000
*
 3000  CONTINUE
*        Add node mass, c.m. and c.m.v. to parent
       IPP=IPAR(IN)
       XNN(IPP)=XNN(IPP)+XNN(IN)
       XMM(IPP)=XMM(IPP)+XMM(IN)
       DO 400 K=1,NDIM
       XMV(K,IPP)=XMV(K,IPP)+XMM(IN)*XMV(K,IN)
 400   XMR(K,IPP)=XMR(K,IPP)+XMM(IN)*XMR(K,IN)
*      IF(IN.EQ.53.OR.IN.EQ.54)THEN
*      PRINT*,' Added to IPP XMV=',IPP,XMV(1,IPP),
*    *    ' ipart=',ILEAVE(IN),XMM(IN)*XMV(1,IN)
*      END IF
*      PRINT*,' Mass ',XMM(IN),' of node ',IN,' added to node ',IPP,
*    *    ' which by now has mass ',XMM(IPP)
*      Jump to next subcell with same parent
*      PRINT*,' We finished cell ',IN,' parent ',IPAR(IN)
       IN=IN+1
*      PRINT*,' Now we go to next cell ',IN,
*    *    ' parent ',IPAR(IN)
       IF(IPAR(IN).EQ.IPAR(IN-1))GOTO 2000
*
*      Now we have got all subcells of this parent
       IN = IPAR(IN-1)
       ILEVEL = ILEVEL - 1
*
*      PRINT*,' Now we go up to level ',ILEVEL,' cell ',IN
       IF(ILEVEL.GT.0)GOTO 3000
*      Prepare final c.m. and c.m.v. of cells
       DO 600 ICC=0,IC
       IF(ILEAVE(ICC).EQ.0)THEN
       DO 601 K=1,NDIM
       XMR(K,ICC) = XMR(K,ICC)/XMM(ICC)
 601   XMV(K,ICC) = XMV(K,ICC)/XMM(ICC)
*      ELSE
*      IPART = ILEAVE(ICC)
*     PRINT*,' Node/Particle ',ICC,ILEAVE(ICC)
*     PRINT*,' Pos=',(XMR(K,ICC),K=1,NDIM)
*     PRINT*,'   X=',(X(K,IPART),K=1,NDIM)
*     PRINT*,' Vel=',(XMV(K,ICC),K=1,NDIM)
*     PRINT*,'   V=',(XDOT(K,IPART),K=1,NDIM)
       END IF
 600   CONTINUE
*
*      PRINT*,' TREE ready IC,XMM=',IC,XMM(0)
*      PRINT*,' XMR=',(XMR(K,0),K=1,NDIM)
*      PRINT*,' XMV=',(XMV(K,0),K=1,NDIM)
*      PRINT*,' IC, XMM=',(ICC,XMM(ICC),ICC=0,IC)
*      PRINT*,' IC, XMM,XMV=',(ICC,XMM(ICC),XMV(1,ICC),ICC=0,IC)
*
       RETURN
       END
