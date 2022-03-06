C
C     COMMON_GRAPE
C
      INTEGER*4 clusterid, NPIPE
      PARAMETER(clusterid=0, NPIPE=48)
      INTEGER*4 new_xunit, new_tunit
      PARAMETER(new_xunit=51, new_tunit=51)
      INTEGER*4 aflag, jflag, pflag, nflag
      PARAMETER(aflag=1, jflag=1, pflag=1, nflag=1)            
      
      INTEGER*4 g6_npipes, g6_set_neighbour_list_sort_mode

      INTEGER*4 g6_read_neighbour_list
      INTEGER*4 g6_get_neighbour_list

      REAl*8 ti, tj0, dtj0
C      REAl*8 tj(NMAX), dtj(NMAX)

      INTEGER*4 g6_calls, nbpipe, nbaver
      REAL*8 nbaver_frac, nbaver_frac0
      REAl*8 a2by18(ID), a1by6(ID), aby2(ID)
      REAl*8 eps2, h2_i(NPIPE)
      REAl*8 x_i(ID,NPIPE), v_i(ID,NPIPE)
      REAl*8 p_i(NPIPE), a_i(ID,NPIPE), jerk_i(ID,NPIPE)
      
      INTEGER*4 nn, nbm, nblen, nbl(LMAX), ret1, ret2, ind_i(NPIPE)
c      PARAMETER(nbm=LMAX-3)
      PARAMETER(nbm=1000)

C     see "verify.f" :-)
            
      INTEGER*4 i, j, k, ii, idi, tmp_i

      REAL*8 YMPI(41,NMAX)
      COMMON/PAR1/YMPI

      INTEGER IMPI(LMAX,NMAX)
      COMMON/PAR2/IMPI
