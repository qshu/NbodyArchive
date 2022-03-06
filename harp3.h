c
c harp3.h
c
c header / common definition for harp3 work area
c
c Copyright Jun Makino 1995
c
c Version 0 95/1/3
      integer IPMAX, KPMAX, G6NCMAX
      parameter (ipmax =192, kpmax = 50, g6ncmax = 4)
      integer h3index(nmax)
      integer h3npipe, h3jpmax
      real * 8 h3eps2(IPMAX)
      real * 8 h3h2(IPMAX)
      real * 8 h3xi(3,IPMAX)
      real * 8 h3vi(3,IPMAX)
      real * 8 h3acc(3,IPMAX)
      real * 8 h3jerk(3,IPMAX)
      real * 8 h3pot(IPMAX)
      real * 8 h3t0(nmax)
      integer g6ncluster
      real * 8 g6fwork(3,ipmax, g6ncmax)
      real * 8 g6jwork(3,ipmax, g6ncmax)
      real * 8 g6pwork(ipmax, g6ncmax)
      integer  g6flagwork(ipmax, g6ncmax)
      integer  g6nj(g6ncmax)
      integer  g6jcid(nmax)
      integer  g6jcloc(nmax), g6error(g6ncmax)
      common/h3int/h3index
      common/h3real/h3eps2,h3h2,h3xi,h3vi,h3acc,h3jerk,h3pot,
     $     h3t0
      common/g6int/g6ncluster, g6nj, g6flagwork, g6jcid, g6jcloc,
     $             g6error
      common/g6real/g6fwork, g6jwork, g6pwork
