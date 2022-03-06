      SUBROUTINE YSAVE(Y)
*
*       Save Y-array in COMMON.
*       -----------------------
*
      INCLUDE 'COMMON1.CH'
      INCLUDE 'COMMON2.CH'
      REAL*8 Y(1)
*
      NC=N-1
      L=0
      DO I=1,4*NC
      L=L+1
      Q(I)=Y(L)
      END DO
      DO I=1,3
      L=L+1
      CMX(I)=Y(L)
      END DO
      L=L+1
      ENERGY=Y(L)
      DO I=1,4*NC
      L=L+1
      P(I)=Y(L)
      END DO
      DO I=1,3
      L=L+1
      CMV(I)=Y(L)
      END DO
      L=L+1
      CHTIME=Y(L)
*
      RETURN
      END
