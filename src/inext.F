      subroutine inext(NQ,LISTQ,TMIN,length,nxtlst)
!VAMPINST include
*
*       Next time-step block.
*       ---------------------
*
      include 'common6.h'
      integer length,j,NQ,nxtlst(*),LISTQ(*)
*
*       Find all new times equal to TMIN.
!VAMPINST subroutine_start
      length = 0
      do 10 L = 1,NQ
         j = LISTQ(L)
         if (TIMENW(j) .eq. TMIN) then
            length = length + 1
            nxtlst(length) = j
         end if
 10   continue
*
!VAMPINST return
          RETURN
      return
!VAMPINST subroutine_end
      end
