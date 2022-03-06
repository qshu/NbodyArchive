      BLOCK DATA BLOCK
*
*
*       Run time initializations.
*       -------------------------
*
      INCLUDE 'params.h'
      REAL*8  EP,DSC,FACM,TFAC,RANGE
      COMMON/RAND/  IY,IFF,IR(97) 
      COMMON/ISAVE/  LI0,LI,NS,NSLIST(LMAX)
      COMMON/ICPU0/  ICPU
      COMMON/IND6/  IND(6)
      COMMON/BSSAVE/  EP(4),DSC,FACM,TFAC,ITFAC,JC
      COMMON/SLOW0/  RANGE,ISLOW(10)
*
*
*       Initialize COMMON indices & B-S data array.
      DATA  IFF,LI0,NS,ICPU  /0,-1,0,0/
      DATA  IND  /1,2,3,4,5,6/
      DATA  EP  /0.04D0,0.0016D0,0.64D-4,0.256D-5/
      DATA  RANGE  /50.0D0/
      DATA  ISLOW  /1,2,4,8,16,32,64,128,256,512/
*
      END
