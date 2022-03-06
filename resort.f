      subroutine resort(length, nxtlst, tm)
c
c     ------------------------------------------
c     maintain particle list sorted by next time
c     ------------------------------------------
c     length (i/o) length of the region to sort / length of
c     the region with minimum timenw
c     nxtlist (i/o) index array
c     tm (o) minimum time
      include 'common1.h'
      integer length,j,nxtlst(*),jj
      real * 8 tmin, tm
c
c     write(6,*)'resort', length, tm
      if (length .eq. 0) then
         length = N
         do 5 j = 1,N
            nxtlst(j) = j
 5       continue
      endif
c     write(6,*)' call sort:', length
      if (length .gt. 1) call sort3(timenw, nxtlst, length)
c     do 30 j = 1,N-1
c        if (timenw(nxtlst(j)) .gt. timenw(nxtlst(j+1))) then
c           write(6,*)'Err:', j, timenw(nxtlst(j)), timenw(nxtlst(j+1))
c        endif
c        write(6,*)j, nxtlst(j), timenw(nxtlst(j))
c30   continue
      tmin = timenw(nxtlst(1))
      length = N
      do 10 j = 1,N
         jj = j
         if (timenw(nxtlst(jj)) .gt. tmin) then
            length = jj - 1
            goto 20
         endif
 10   continue
 20   continue
      tm = tmin
      return
      end
