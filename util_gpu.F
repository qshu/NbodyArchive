C--------------------------------------------------------------------
C     Open the local GPU
C--------------------------------------------------------------------
      SUBROUTINE OPEN_THE_GPU

      INCLUDE 'common6.h'
      INCLUDE 'common_grape.h'
      
C---  open GPU and set some parameters

      g6_calls = 0

      nbaver = 0
      nbpipe = 0

      CALL g6_open(clusterid)
      
c      WRITE(6,*) 'OPEN start', RANK, clusterid
      

      npmax = g6_npipes()

c      WRITE(6,*) 'OPEN start', RANK, npmax

      
CCC      IF(npmax.NE.NPIPE) WRITE(STDOUT,*) '!!! NPIPE is NOT 48 !!!'
CCC      CALL flush(STDOUT)
      
      IF(ipipe.GT.npmax) ipipe = npmax
      
      np_use = ipipe

c      WRITE(6,*) 'OPEN start', RANK, np_use, ipipe, npmax


      
      CALL g6_set_tunit(new_tunit)
      CALL g6_set_xunit(new_xunit)

      CALL g6_set_overflow_flag_test_mode(aflag, jflag, pflag)

      tmp_i = g6_set_neighbour_list_sort_mode(nflag)

C---  open GPU and set some parameters
      
      RETURN
      END
C--------------------------------------------------------------------
C     Open the local GPU
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     Close the local GPU
C--------------------------------------------------------------------
      SUBROUTINE CLOSE_THE_GPU

      INCLUDE 'common6.h'
      INCLUDE 'common_grape.h'
            
C---  close the GPU

      CALL g6_close(clusterid) 

C---  close the GPU
      
      RETURN
      END
C--------------------------------------------------------------------
C     Open the local GPU
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     Set the ALL particles on local GPU memory
C--------------------------------------------------------------------
      SUBROUTINE SET_ALL_PART_ON_GPU

      INCLUDE 'common6.h'
      INCLUDE 'common_grape.h'
            
C---  SET the particles on local GPU memory

      tj0  = 0.0D0
      dtj0 = 2.0D0**(-25)
      
C      DO j=1,N
C         tj(j)  = 0.0
C         dtj(j) = 2.0**(-25)
C      ENDDO					! j
 
      DO k=1,ID
         a2by18(k) = 0.0D0
         a1by6(k)  = 0.0D0 
         aby2(k)   = 0.0D0
      ENDDO					! k

c      WRITE(6,*) 'SETALL start', RANK, IFIRST, NTOT, ID, 
c     &BODY(IFIRST), X(1,IFIRST), XDOT(1,IFIRST), 
c     &BODY(NTOT), X(1,NTOT), XDOT(1,NTOT)

      DO j=IFIRST,NTOT

         CALL g6_set_j_particle(clusterid, j-IFIRST, j,
     &   tj0, dtj0,
     &   BODY(j),
     &   a2by18(1), a1by6(1), aby2(1),
     &   XDOT(1,j), X(1,j)) 

      ENDDO					! j
      
      RETURN
      END
C--------------------------------------------------------------------
C     Set the ALL particles on local GPU memory
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     Set the ACT particles on local GPU memory
C--------------------------------------------------------------------
      SUBROUTINE SET_ACT_PART_ON_GPU(IREG,NREG)

      INCLUDE 'common6.h'
      INCLUDE 'common_grape.h'

      INTEGER*4 IREG(NMAX), NREG

            
C---  SET the particles on local GPU memory

      tj0  = 0.0D0
      dtj0 = 2.0D0**(-25)
      
C      DO j=1,N
C         tj(j)  = 0.0
C         dtj(j) = 2.0**(-25)
C      ENDDO					! j
 
      DO k=1,ID
         a2by18(k) = 0.0D0
         a1by6(k)  = 0.0D0 
         aby2(k)   = 0.0D0
      ENDDO					! k


      DO L=1,NREG

         j=IREG(L)

         CALL g6_set_j_particle(clusterid, j-IFIRST, j,
     &   tj0, dtj0,
     &   BODY(j),
     &   a2by18(1), a1by6(1), aby2(1),
     &   XDOT(1,j), X(1,j)) 

      ENDDO					! L

      
      RETURN
      END
C--------------------------------------------------------------------
C     Set the ACT particles on local GPU memory
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     Calculate the REG force on local GPU memory
C--------------------------------------------------------------------
      SUBROUTINE CALC_REG_ON_GPU(IREG,istart,iend)

      INCLUDE 'common6.h'
      INCLUDE 'common_grape.h'

      INTEGER*4 IREG(NMAX), istart, iend, LIST_M(NMAX)
      INTEGER*4 LIST_G(LMAX,NMAX)


C---  calculate the pot, a & jerk on the GPU

      DO L=istart,iend,np_use

c      WRITE(6,*) 'LOOP REG', RANK, L, LIST(1,1)

         nn = np_use
         IF( (L+nn).GT.iend ) nn = iend - L + 1



CCC   CALL NBMERGE!!!	Present version LIST_M == LIST


      DO ii=1,nn

      idi = IREG(L+ii-1)

        LIST_M(1) = LIST(1,idi)

        DO tmp_i=1,LIST(1,idi)
           LIST_M(tmp_i+1) = LIST(tmp_i+1,idi)
        ENDDO
	
      END DO




C---  load the 0 masses for the "i" particle neighbours on the local GPU's 

        tmp = 0.0D0

        DO tmp_i=1,LIST_M(1)

          j = LIST_M(tmp_i+1)

CCC          CALL g6_set_j_particle_mxonly(clusterid, j-IFIRST, j, 
CCC     &                                 tmp, X(1,j))

C          tj(j)  = 0.0
C          dtj(j) = 2.0**(-25)

          tj0  = 0.0D0
          dtj0 = 2.0D0**(-25)
 
          DO k=1,ID
            a2by18(k) = 0.0D0
            a1by6(k)  = 0.0D0 
            aby2(k)   = 0.0D0
          ENDDO					! k
     
         CALL g6_set_j_particle(clusterid, j-IFIRST, j,
     &   tj0, dtj0,
     &   tmp,
     &   a2by18(1), a1by6(1), aby2(1),
     &   XDOT(1,j), X(1,j)) 

        ENDDO					! j




C---  fill the GPU "working" arrays

         DO ii=1,nn

         idi       = IREG(L+ii-1)
ccc         ind_i(ii) = NAME(idi)
         ind_i(ii) = idi

CCC      Obtain irregular & regular force and determine current neighbours.
CCC      Take volume between inner and outer radius equal to basic sphere.

         RS2      = RS(idi)**2
	 RCRIT2   = 1.05*RS2
         h2_i(ii) = RCRIT2

C      WRITE(6,*) RANK, 'GPU UTIL h2', ii-1, idi, sqrt(h2_i(ii)), 
C     &NAME(idi), LIST(1,idi)

         DO k=1,3
            x_i(k,ii) = X(k,idi)
            v_i(k,ii) = XDOT(k,idi)
         ENDDO				! k

         p_i(ii) = PHIDBL(idi)

         DO k=1,3
            a_i(k,ii)    = 2.0D0*F(k,idi)
            jerk_i(k,ii) = 6.0D0*FDOT(k,idi)
         ENDDO				! k

         ENDDO					! ii



CCC GPU manual for NB search with smal pipeline number

      IF(ipipe.LT.npmax) THEN
      
      DO ii = ipipe+1, npmax

      h2_i(ii) = 0.0D0

      DO k=1,3
        x_i(k,ii) = X(k,IREG(L))
      ENDDO				! k
      
      END DO
      
      END IF




CCC      WRITE(6,*) ii-1, idi, h2_i(48), NAME(idi)


C---  GPU first call with first values of p_i, a_i & jerk_i

        CALL g6calc_firsthalf(clusterid, N-NPAIRS, nn, ind_i, 
     &	               x_i, v_i,
     &                 a_i, jerk_i, p_i, eps2_g6a, h2_i)


        CALL g6calc_lasthalf(clusterid, N-NPAIRS, nn, ind_i, 
     &                 x_i, v_i,
     &                 eps2_g6a, h2_i, a_i, jerk_i, p_i)


        g6_calls = g6_calls + 1


C---  before second GPU call "correct" the 0.0 jerk values

C        IF(jerk_i(1,ii).EQ.0.0) jerk_i(1,ii) = 1.0E-05
C        IF(jerk_i(2,ii).EQ.0.0) jerk_i(2,ii) = 1.0E-05
C        IF(jerk_i(3,ii).EQ.0.0) jerk_i(3,ii) = 1.0E-05


C---  GPU second call with new values of p_i, a_i & jerk_i

C        CALL g6calc_firsthalf(clusterid, N, nn, ind_i, 
C     &	               x_i, v_i,
C     &                 a_i, jerk_i, p_i, eps2_g6a, h2_i)


C        CALL g6calc_lasthalf(clusterid, N, nn, ind_i, 
C     &                 x_i, v_i,
C     &                 eps2_g6a, h2_i, a_i, jerk_i, p_i)




CCC#ifdef NBCALC

C---  read the i particles neighbours

        ret1 = g6_read_neighbour_list(clusterid)


        IF(ret1.NE.0) THEN
        WRITE(6,*) 'NB list ret1', clusterid, ret1
	END IF 

        IF(ret1.EQ.0) THEN

	  nbpipe = nbpipe + nn

          DO ii=1,nn

            ret2 = g6_get_neighbour_list(clusterid, ii-1, 
     &                                   nbm, nblen, nbl)



          IF(ret2.EQ.0) THEN

C--- read the nblist

            idi         = IREG(L+ii-1)
            LIST_G(1,idi) = nblen

      IF(nblen.GT.LMAX-1) WRITE(6,*)'WARNING util_grape!!',nblen,LMAX-1

             DO k=1,nblen
               LIST_G(k+1,idi) = nbl(k)
             ENDDO				! k

          nbaver = nbaver + nblen



C--- Double check the GPU nblist for particle "idi" !!!

      DO tmp_i=1,LIST_G(1,idi)

        j = LIST_G(tmp_i+1,idi)

        tmp = ( X(1,idi) - X(1,j) )**2 + 
     &	      ( X(2,idi) - X(2,j) )**2 + 
     &	      ( X(3,idi) - X(3,j) )**2

        IF(((sqrt(tmp)-sqrt(RCRIT2))/sqrt(RCRIT2)).GT.1.0D-06) THEN
	
      WRITE(6,*) 'NB list ', rank, g6_calls, L, idi, 
     &NAME(idi), LIST_G(1,idi), 
     &j, NAME(j),
     &sqrt(tmp), sqrt(RCRIT2), ((sqrt(tmp)-sqrt(RCRIT2))/sqrt(RCRIT2))
	
	ENDIF
	  
      ENDDO					! j




          ENDIF ! IF(ret2.EQ.0)



          IF(ret2.NE.0) THEN
          WRITE(6,*) 'NB list rank, ret2, nbm, nblen', 
     &	  rank, ret2, nbm, nblen
          END IF


          IF(ret2.EQ.1) THEN	    
c	  WRITE(6,*) 'Length of NB list > LMAX', rank, ret2
c          CALL flush(6);
	  
C--- Set the negative value as a ERROR flag !!!

          idi         = IREG(L+ii-1)
          LIST_G(1,idi) = -1

          DO k=1,nblen
              LIST_G(k+1,idi) = -1
            ENDDO				! k
	      
          ENDIF	! IF(ret2.EQ.1)
	    


          ENDDO					! ii
        

	nbaver_frac0 = DBLE(nbaver-nbaver_old)/DBLE(nn)
        nbaver_frac = DBLE(nbaver)/DBLE(nbpipe)
	nbaver_old = nbaver

c        IF(mod(g6_calls,1000).EQ.0) WRITE(6,*) ' NB aver ', 
c     &nbaver_frac0, nbaver_frac, ' nn aver ', 
c     &DBLE(nbpipe)/DBLE(g6_calls)


        ELSE 

          IF(ret1.GT.0) WRITE(6,*) 'NB list overflow in hardware', 
     &	  rank, ret1

          IF(ret1.LT.0) WRITE(6,*) 'NB list hardware error', 
     &	  rank, ret1
          CALL flush(6);

        ENDIF ! IF(ret1.EQ.0)


C---  read the i particles neighbours

CCC#endif





C---  get the pot, a & adot values from the GPU working arrays

        DO ii=1,nn

        idi      = IREG(L+ii-1)
        PHIDBL(idi) = p_i(ii)

        DO k=1,3
           FRG(k,idi)    = a_i(k,ii)
           FRGDOT(k,idi) = jerk_i(k,ii)
        ENDDO				! k


c      WRITE(6,*) 'FORCE REG', RANK, ii, p_i(ii), a_i(1,ii), jerk_i(1,ii)


#ifndef PARALLEL
c save new GPU neighbourlist on LIST
          DO 655 K = 1,LMAX
          KLISTX(K,idi) = LIST_G(K,idi)
c          IMPI(K,idi) = LIST(K,idi)
          IMPI(K,L) = LIST(K,idi)
c          LIST(K,idi) = LIST_G(K,idi)
 655      CONTINUE

C We after need to pass the correct LIST_G to regint!!!

          NBSUM = 0

       call cputim(tt61)
c       CALL REGINT(idi,IMPI(1,idi))
       CALL REGINT(idi,IMPI(1,L))       
       call cputim(tt62)
       ttgrcalc = ttgrcalc - (tt62-tt61)*60.0
       ttreg = ttreg + (tt62-tt61)*60.0


         DO 615 K = 1,3
              F(K,idi) = 0.5D0*(FI(K,idi) + FR(K,idi))
              FDOT(K,idi) = ONE6*(FIDOT(K,idi) + FRDOT(K,idi))
  615    CONTINUE

#else

*       Start PARALLEL section for regint

          call cputim(tt63)          
c save new GPU neighbourlist on LIST 
          DO 656 K = 1,LMAX
          KLISTX(K,idi) = LIST_G(K,idi)
c          IMPI(K,idi) = LIST(K,idi)
          IMPI(K,L) = LIST(K,idi)	  
c          LIST(K,idi) = LIST_G(K,idi)
 656      CONTINUE
*

c       CALL REGINT(idi,IMPI(1,idi))
       CALL REGINT(idi,IMPI(1,L))       

*
          YMPI(37,L) = STEPR(idi)
          YMPI(38,L) = STEP(idi)
          YMPI(39,L) = TIMENW(idi)
          YMPI(40,L) = RS(idi)
          YMPI(41,L) = T0R(idi)
*
         DO 814 K = 1,3
              YMPI(K,L) = XN(K,idi)
              YMPI(K+3,L) = XNDOT(K,idi)
              YMPI(K+6,L) = FI(K,idi)
              YMPI(K+9,L) = FR(K,idi)
              YMPI(K+12,L) = FIDOT(K,idi)
              YMPI(K+15,L) = FRDOT(K,idi)
              YMPI(K+18,L) = D1R(K,idi)
              YMPI(K+21,L) = D2R(K,idi)
              YMPI(K+24,L) = D3R(K,idi)
*       Save corrections of force polynomials from routine fpcorr
              IF (KZ(38).GT.0.OR.idi.GT.N) THEN
              YMPI(K+27,L) = D1(K,idi)
              YMPI(K+30,L) = D2(K,idi)
              YMPI(K+33,L) = D3(K,idi)
              END IF
  814    CONTINUE

       call cputim(tt64)
       ttgrcalc = ttgrcalc - (tt64-tt63)*60.0
       ttreg = ttreg + (tt64-tt63)*60.0

#endif      
      


CCC        DO k=1,ID
CCC          X(k,idi) = XN(k,idi)
CCC          XDOT(k,idi) = XNDOT(k,idi)
CCC        ENDDO					! k


C       Correction of reg force list LIST_M & LIST


        ENDDO					! ii




C---  load back the masses for the "i" particle neighbours on the local GPU's 

        DO tmp_i=1,LIST_M(1)

          j = LIST_M(tmp_i+1)
          
C          CALL g6_set_j_particle_mxonly(clusterid, j-IFIRST, j, 
C     &                                 BODY(j), X0(1,j))
          
C          tj(j)  = 0.0
C          dtj(j) = 2.0**(-25)

          tj0  = 0.0D0
          dtj0 = 2.0D0**(-25)

          DO k=1,ID
            a2by18(k) = 0.0D0
            a1by6(k)  = 0.0D0 
            aby2(k)   = 0.0D0
          ENDDO					! k

         CALL g6_set_j_particle(clusterid, j-IFIRST, j,
     &   tj0, dtj0,
     &   BODY(j),
     &   a2by18(1), a1by6(1), aby2(1),
     &   XDOT(1,j), X(1,j)) 

        ENDDO					! j



      ENDDO					! L


      RETURN
      END

C--------------------------------------------------------------------
C     Calculate the REG force on local GPU memory
C--------------------------------------------------------------------
 
