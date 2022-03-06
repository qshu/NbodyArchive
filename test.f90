program test
  implicit none
  integer, parameter :: isize = 4
  integer, parameter :: nl    = 55
  integer :: inl,idiff,irun,jsize,ix,inum(isize),ista(isize)
  integer :: rank,isend,irecv,ir,irank,istsen,icnt,istrec,icnt2

  inl = nl/isize
  jsize = isize*inl
  idiff = nl - jsize
  irun = 0
  do ix = 1,isize
    inum(ix)=inl
    if(ix.le.idiff) inum(ix) = inum(ix) + 1
    ista(ix) = irun+1
    if(ista(ix).gt.nl) inum(ix) = 0
    irun = irun + inum(ix)
    write(*,*) ix,inum(ix),ista(ix)
  end do

  do rank =0,isize-1
     isend = rank + 1
     if(isend.eq.isize) isend = 0
     irecv = rank - 1
     if(irecv.eq.-1)irecv = isize - 1
     write(*,"('rank ',i0,' sends to ',i0,' and receive from ',i0)") rank,isend,irecv
     do ir = 0,isize-2
        irank = rank - ir
        if(irank.lt.0)irank=irank+isize
        istsen=ista(irank+1)
        icnt = inum(irank+1)
        if(irank.eq.0)irank=isize
        istrec = ista(irank)
        icnt2 = inum(irank)
        write(*,"('sent message ',i0,' start at ',i0,' with size ',i0)") ir,istsen,icnt
        write(*,"('recv message ',i0,' start at ',i0,' with size ',i0)") ir,istrec,icnt2
     end do
   end do

end program test
