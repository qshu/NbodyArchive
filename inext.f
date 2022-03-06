      subroutine inext(i, length, iloc, nxtlst, timenw, tmin)
c
c     --------------------------------------------------
c     inext.f : next body finder for integrator
c
c     By J. Makino
c     Ver 0.00 92/5/23
c     --------------------------------------------------
c
      integer i,length,ipos,iloc,nxtlst(*)
      real*8 timenw(*),tmin
      save ipos
c
c     print*,' inext: i, length, ipos, tmin=',i,length,ipos,tmin
      if ((length .eq. 0) .or. (ipos .eq. length + 1)) then
         call resort(length,nxtlst,tmin)
         ipos = 1
      endif
      iloc = ipos
      i = nxtlst(iloc)
      tmin = timenw(i)
      ipos = ipos + 1
      return
      end
