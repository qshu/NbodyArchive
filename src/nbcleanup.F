      SUBROUTINE NBCLEANUP
      INCLUDE 'common6.h'

!VAMPINST include


#ifdef VISIT1
!      CALL FLVISIT_NBODY1_CLOSE()
#endif
!VAMPINST subroutine_start
#ifdef VISIT2
      CALL VISCONCLOSE()
#endif
        write (*,*) 'ending nbody'
#ifdef USEHPM
        call f_hpmstop(1);
        call f_hpmterminate_(taskID);
#endif
#ifdef PARALLEL
       write(*,*) 'before MPI_FINALIZE'
       

!VAMPINST region_end
!VTREGION END wfcleanup

!VAMPINST region_end
!VTREGION END nbody
       
       call MPI_FINALIZE(ierr)
#endif

   99 FORMAT (//,' WFBENCH myrank=',I3,' NPE=',I5,' INITTIME =',F14.8,
     &         ' LOOPTIME =',
     &         F14.8,' OVERALLTIME =',F14.8)

!VAMPINST subroutine_end
      END
