


      SUBROUTINE OLDMYDUMP(I,J)
*
*
*       COMMON save or read.
*       --------------------
*
* Updated on 94/08/06 for BHcode.
*
      INCLUDE 'params.h'
      integer na,nb,nd
      PARAMETER  (NA=22*NMAX,NB=16,ND=NMAX+66)
      REAL*8  A,C,D, dummy
      INTEGER IB,i,j,k, l1, l2,l
*
      COMMON/NBODY/  A(NA)
      COMMON/NAMES/  IB(NB)
      COMMON/PARAMS/ C(20)
      COMMON/BLOCKS/ D(ND)
*
*
*       Open unit #J by reading dummy and rewinding.
*
*     original scheme does not work on Alpha. I'm trying
*     explicit call to open here
*
c      READ (J,ERR=10,END=10)  DUMMY
c   10 REWIND J
*
*       Read or save all COMMON variables (valid for tape or disc).
      IF (I.EQ.0) THEN
          open(unit=j,access='SEQUENTIAL',form='UNFORMATTED')
c          READ (J)   A
          do l = 1, 22
             l1 = (l-1)*(nmax)+1
             l2 = l*nmax
             read (J)  (A(k),k=l1,l2)
          enddo
          READ (J)   IB
          READ (J)   C
          READ (J)   D
          close(j)
      ELSE
          if( j.eq.1) then
              call unlink('fort.1')
              open(unit=j,name='fort.1',access='SEQUENTIAL',
     $             status='UNKNOWN',
     $        form='UNFORMATTED')
          endif
          if( j.eq.2) then
              call unlink('fort.2')
              open(unit=j,name='fort.2',access='SEQUENTIAL',
     $             status='UNKNOWN',
     $        form='UNFORMATTED')
          endif
          if( j.eq.3) then
              call unlink('fort.3')
              open(unit=j,name='fort.3',access='SEQUENTIAL',
     $             status='UNKNOWN',
     $        form='UNFORMATTED')
          endif
          do l = 1, 22
             l1 = (l-1)*(nmax)+1
             l2 = l*nmax
             WRITE (J)  (A(k),k=l1,l2)
          enddo
          write(6,*)'Mydump - A written'
          call flush(6)
          WRITE (J)  IB
          write(6,*)'Mydump - B written'
          call flush(6)
          WRITE (J)  C
          write(6,*)'Mydump - C written'
          call flush(6)
          WRITE (J)  D
          write(6,*)'Mydump - D written'
          call flush(6)
          CLOSE (UNIT=J)
      END IF
*
      RETURN
*
      END



      SUBROUTINE MYDUMP(I,J)
*
*
*       COMMON save or read.
*       --------------------
*
* Updated on 94/08/06 for BHcode.
*
      INCLUDE 'params.h'
      integer na,nb,nd, nc
      PARAMETER  (NA=56*NMAX,NB=16,nc=40,ND=NMAX*2+132)
      REAL*4  A,C,D
      INTEGER IB,i,j,k, l1, l2,l, ii
*     
      COMMON/NBODY/  A(NA)
      COMMON/NAMES/  IB(NB)
      COMMON/PARAMS/ C(nc)
      COMMON/BLOCKS/ D(ND)
*
*     Read or save all COMMON variables (valid for tape or disc).
      call open_dump(j,i)
      IF (I.EQ.0) THEN
         do ii = 1,56
            call read_array(a(1+(ii-1)*nmax),nmax)
         enddo
         call read_array(ib,nb)
         call read_array(c,nc)
         call read_array(d,nd)
      ELSE
         do ii = 1,56
            call write_array(a(1+(ii-1)*nmax),nmax)
         enddo
         call write_array(ib,nb)
         call write_array(c,nc)
         call write_array(d,nd)
      END IF
      call close_dump
*     
      RETURN
*
      END



