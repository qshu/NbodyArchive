***
*
* grape6.h
*
* header / common definition for GRAPE work area
*
***
      INTEGER ipmax
      PARAMETER(ipmax=192)
*
      INTEGER g6_npipes,g6_calc_lasthalf
      INTEGER g6_read_neighbour_list,g6_get_neighbour_list
      INTEGER gperr
      INTEGER gpaddr(nmax),gpindx(nmax)
      REAL*8 gpt0,gpdtj,gpmass
      REAL*8 gpacc(3,ipmax),gpjerk(3,ipmax)
      REAL*8 gppot(ipmax),gph2(ipmax)
***
