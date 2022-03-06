      PROGRAM group
*
      include 'mpif.h'

      integer ierr, grank, gsize, colour
      integer yellow, green, red, black
      integer key, colour_comm
      integer size1,size2,size3,size4
      integer rank1,rank2,rank3,rank4,itest

*     Initialize MPI
*     -----------------------
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,grank,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,gsize,IERR)
      PRINT*,' This is grank=',grank,' gsize=',gsize,
     *                          '   MPI_COMM_WORLD'

      yellow = 1
      green  = 2
      red    = 3
      black  = 4

      key = grank  ! bestimmt die Reihenfolge der ranks in den Untergruppen
c                    (hier:gleiche Reihenfolge wie in der Gesamtgruppe
      

c     einfaerben der Prozessoren,
c     ----------------------------------------
      if((grank.ge.0).and.(grank.lt.gsize/2))then
       colour=yellow
      endif                

      if((grank.ge.gsize/2).and.(grank.lt.gsize))then
       colour=green
      endif

c     splitten von MPI_COMM_WORLD
c     -----------------------------------------------
      call MPI_COMM_SPLIT(MPI_COMM_WORLD,colour,key,
     $                           colour_comm, ierr)


c     fertig, es wurden  soviele Gruppen erzeugt wie es farben gibt

c     ##################################################
c     tests: Zugriff auf die Untergruppen
c     ###################################################

c     Untergruppe yellow
c     ---------------------

      IF(COLOUR.EQ.YELLOW) THEN

       CALL MPI_COMM_RANK(colour_comm,rank1,ierr)
       CALL MPI_COMM_SIZE(colour_comm,size1,ierr)
       print*,'size of yellow-group = ', size1
       print*,'rank of yellow-process=', rank1
     *       ,'   grank=',grank

c      test the collective communication in the group
       itest = 1
       if(rank1.eq.0) itest=1000
       call MPI_BCAST(itest,1,MPI_INTEGER,0,colour_comm,ierr)

c      Do some NBODY-Code with yellow Processors

      ENDIF

c     Untergruppe green
c     ---------------------

      IF(COLOUR.EQ.GREEN) THEN

       CALL MPI_COMM_SIZE(colour_comm,size2,ierr)
       CALL MPI_COMM_RANK(colour_comm,rank2,ierr)
       print*,'size of green-group= ', size2
       print*,'rank of green-process=', rank2
     *       ,'   grank=',grank

c      test the collective communication in the group
       itest = 1
       if(rank2.eq.0) itest=2000
       call MPI_BCAST(itest,1,MPI_INTEGER,0,colour_comm,ierr)

c      Do some NBODY-Code with green Processors

      ENDIF

c     folgende Barriere ist nur fuer synch des outputs
      call MPI_BARRIER(MPI_COMM_WOLD,ierr)

      print*,'itest= ',itest, grank

      call mpi_finalize(ierror)
      stop

*     Test: globale communication in MPI_COMM_WORLD
*     ================================================
      itest=1
      if(grank.eq.0) itest=1111
      call MPI_BCAST(itest,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      print*,'itest= ',itest, grank

*     Finalize MPI
*     -----------------------------------------
      CALL MPI_FINALIZE(ierror)

*
      STOP
      END
