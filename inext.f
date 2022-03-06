      subroutine inext(i, length, iloc, nxtlst, timenw, tmin)
c
c     --------------------------------------------------
c     inext.f : next body finder for integrator
c
c     By J. Makino
c     Ver 0.00 92/5/23
c
c     Normal usage : first call with length = 0
c     to initialize the sorted list. Then keep calling this
c     routine.
c
c     Non-standard usage to find all particles with minimum time:
c     First call: Length = 0.
c     2nd and later call : Length <- -length
c     return value of       
c     length gives the number of
c     particles that share the minimum time
c     
c     --------------------------------------------------
c
      integer i,length,ipos,iloc,nxtlst(*)
      real*8 timenw(*),tmin
      save ipos
c
      if ((length .le. 0) .or. (ipos .eq. length + 1)) then
         if(length .lt. 0) then
            length = -length
         endif
c         write(6,*)' before resort:', length, tmin
         call resort(length,nxtlst,tmin)
c         write(6,*)' after resort:', length, tmin
         ipos = 1
      endif
      iloc = ipos
      i = nxtlst(iloc)
      tmin = timenw(i)
      ipos = ipos + 1
      return
      end
