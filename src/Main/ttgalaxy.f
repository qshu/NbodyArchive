      SUBROUTINE TTGALAXY(X,Y,Z,T,MSC,RSC,TSC,VSC,TTPHIG,TTTDEP,IKEYPOT)
*
*       User-defined galactic potential
*       -------------------------------

*** FlorentR - new subroutine

* Compute the galactic potential TTPHIG as a function of the position
* X, Y, Z and time T.
*
*
************ IMPORTANT NOTES:    (please read carefully)
*
* common6 is not included to avoid the user to name its variables as
* those of the common6. Instead, the parameters are passed as arguments.
*
* To convert a value in physical units into nbody6 units, *divide* the
* physical quantity in Msun, pc, Myr, km/s
* by: MSC,RSC,TSC,VSC for the length, mass, velocity and time
* respectively.
*
* Keep in mind that this routine is executed 25 times for every
* star and tail member, at every timestep. That's a lot! So make sure
* you optimize your code!
* Avoid doing the same computation twice. Put as much as possible
* in the "FIRST" block (see below).
*
* Friendly advice: do not convert X, Y, Z, T
* into physical units (kpc, Myr). Convert the parameters instead, in the
* "FIRST" block if possible.
*
* Second friendly advice (because I am very friendly!): in subroutines,
* do not overwrite the potential. Do this instead: ttphig = ttphig + ...
* This way, you can easily sum up multiple components!
*
      INTEGER TTTDEP,IKEYPOT(11)
      REAL*8 X,Y,Z,T, MSC, RSC, TSC, VSC, TTPHIG

      TTPHIG = 0.D0 ! init (= no galaxy)

      IF(IKEYPOT(1).EQ.1) CALL pointmass(ttphig,x,y,z,t,msc,rsc,tsc,vsc)
      IF(IKEYPOT(2).EQ.1) CALL plummer(ttphig,x,y,z,t,msc,rsc,tsc,vsc)
      IF(IKEYPOT(3).EQ.1) CALL hernquist(ttphig,x,y,z,t,msc,rsc,tsc,vsc)
      IF(IKEYPOT(4).EQ.1) CALL nfw(ttphig,x,y,z,t,msc,rsc,tsc,vsc)
      IF(IKEYPOT(5).EQ.1) CALL nfw2(ttphig,x,y,z,t,msc,rsc,tsc,vsc)
      IF(IKEYPOT(6).EQ.1) CALL isothermal(ttphig,x,y,z,t,msc,
     &     rsc,tsc,vsc)
      IF(IKEYPOT(7).EQ.1) CALL nakedmndisk(ttphig,x,y,z,t,msc,
     &     rsc,tsc,vsc)
      IF(IKEYPOT(8).EQ.1) CALL pseudoexpdisk(ttphig,x,y,z,t,msc,
     &     rsc,tsc,vsc)
      IF(IKEYPOT(9).EQ.1) CALL bulgediskhalo(ttphig,x,y,z,t,msc,
     &     rsc,tsc,vsc)
      IF(IKEYPOT(10).EQ.1) CALL spiralarms(ttphig,x,y,z,t,msc,
     &     rsc,tsc,vsc)
      IF(IKEYPOT(11).EQ.1) CALL nfwcosmo(ttphig,x,y,z,t,msc,rsc,tsc,vsc)


* set TTTDEP to 0 if the potential is time-independent or adiabatically
* (slowly compared to the motion of the cluster) changing.
* Set TTTDEP to 1 if the potential is strongly time-dependent.

      TTTDEP = 0

      RETURN
      END

************************************************************************
************************************************************************
************************************************************************
      SUBROUTINE pointmass(ttphig,x,y,z,t,msc,rsc,tsc,vsc)

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, msc, rsc, tsc, vsc

      real*8 mass
      save mass

      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here
      
        mass = 1.d10 / msc ! Msun -> Nbody Units

        first = .false.
      ENDIF

* compute variables here
      ttphig = ttphig - mass / sqrt(x**2+y**2+z**2)

      RETURN
      END

************************************************************************
      SUBROUTINE plummer(ttphig,x,y,z,t,msc,rsc,tsc,vsc)

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, msc, rsc, tsc, vsc

      real*8 m, r02
      save m, r02

      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here

        m = 1.d13 / msc  ! Msun -> Nbody Units
        r02 = (3.d3 / rsc)**2 ! pc -> Nbody Units

        first = .false.
      ENDIF

* compute variables here
      ttphig = ttphig - m / sqrt(r02 + x**2+y**2+z**2)
      
      RETURN
      END
      
************************************************************************
      SUBROUTINE hernquist(ttphig,x,y,z,t,msc,rsc,tsc,vsc)

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, msc, rsc, tsc, vsc

      real*8 m, r0
      save m, r0

      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here

        m = 1.d13 / msc  ! Msun -> Nbody Units
        r0 = 3.d3 / rsc ! pc -> Nbody Units

        first = .false.
      ENDIF

* compute variables here

      ttphig = ttphig - m / (r0 + sqrt(x**2+y**2+z**2))
      
      RETURN
      END

************************************************************************
      SUBROUTINE nfw(ttphig,x,y,z,t,msc,rsc,tsc,vsc)

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, msc, rsc, tsc, vsc

      real*8 m, r0, r
      save m, r0
 
      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here
      
        m = 1.d13 / msc  ! Msun -> Nbody Units
        r0 = 30.d3 / rsc ! pc -> Nbody Units

        first = .false.
      ENDIF

* compute variables here
      r = sqrt(x**2+y**2+z**2)
      ttphig = ttphig - m / r * log(1.D0 + r/r0)

      RETURN
      END

************************************************************************
      SUBROUTINE nfw2(ttphig,x,y,z,t,msc,rsc,tsc,vsc)

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, msc, rsc, tsc, vsc

      real*8 mvir, rs, H0, c, r200, pi, r
      save mvir, rs, c

      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here
      
        mvir = 1.23805d+12 / msc  ! Msun -> Nbody units
        H0 = 71.d0 * 1.d-6 / vsc   ! km/s/Mpc -> km/s/pc ->  Nbody
        c = 15.8163d0 ! dimensionless
        pi = atan(1.d0)*4.d0

        r200 = (mvir / ( (4.d0/3.d0)*pi*200.d0 * 
     &     3.d0*H0**2/(8.d0*pi) ))**(1.d0/3.d0)
        rs = r200/c ! already in Nbody units

        first = .false.
      ENDIF

* compute variables here
      r = sqrt(x**2+y**2+z**2)
      ttphig = ttphig - mvir / r * log(1.d0 + r/rs) 
     &  / (log(1.d0 + c) - c/ (1.d0 + c))

      
      RETURN
      END

************************************************************************
      SUBROUTINE isothermal(ttphig,x,y,z,t,msc,rsc,tsc,vsc)

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, msc, rsc, tsc, vsc

      real*8 v02H
      save v02H

      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here
      
        v02H = 0.5d0 * (220.d0 / vsc)**2   ! Km/s -> Nbody units

        first = .false.
      ENDIF

* compute variables here
      ttphig = ttphig + v02H * log(x**2+y**2+z**2)
      
      RETURN
      END


************************************************************************
      SUBROUTINE nakedmndisk(ttphig,x,y,z,t,msc,rsc,tsc,vsc)

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, msc, rsc, tsc, vsc

      real*8 md, a, b2
      save md, a, b2

      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here

        md = 1.d11 / msc   ! Msun -> Nbody units
        a = 5.d3 / rsc       ! pc -> Nbody units
        b2 = (300.d0 / rsc)**2  ! pc -> Nbody units

        first = .false.
      ENDIF
      
      ttphig = ttphig - md/sqrt(x**2+y**2+(a+sqrt(z**2+b2))**2)
      
      RETURN
      END

************************************************************************
      SUBROUTINE pseudoexpdisk(ttphig,x,y,z,t,msc,rsc,tsc,vsc)
! mimic an exponential disk with 3 Miyamoto-Nagai disks
! Details in Smith et al. (2015)  

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, msc, rsc, tsc, vsc

      real*8 md, rd, hz, hzrd, brd, m1, m2, m3, a1, a2, a3, b2
      save m1, m2, m3, a1, a2, a3, b2

      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here
*      exponential disk parameters
        md = 1.d11 / msc ! Msun -> Nbody units
        rd = 5.d3 / rsc ! pc -> Nbody units
        hz = 1.d2 / rsc! pc -> Nbody units
        
* MN b (from Smith et al. 2015, Fig 5)
        hzrd = hz / rd
        brd = -0.269d0* hzrd**3 + 1.080d0* hzrd**2 + 1.092d0 * hzrd
        
* MN parameters (from Smith et al. 2015, Table 2 and Eq 7)
        m1= 0.0036d0*brd**4-0.0330d0*brd**3
     &      +0.1117d0*brd**2-0.1335d0*brd+0.1749d0
        m2=-0.0131d0*brd**4+0.1090d0*brd**3
     &      -0.3035d0*brd**2+0.2921d0*brd-5.7976d0
        m3=-0.0048d0*brd**4+0.0454d0*brd**3
     &      -0.1425d0*brd**2+0.1012d0*brd+6.7120d0
        a1=-0.0158d0*brd**4+0.0993d0*brd**3
     &      -0.2070d0*brd**2-0.7089d0*brd+0.6445d0
        a2=-0.0319d0*brd**4+0.1514d0*brd**3
     &      -0.1279d0*brd**2-0.9325d0*brd+2.6836d0
        a3=-0.0326d0*brd**4+0.1816d0*brd**3
     &      -0.2943d0*brd**2-0.6329d0*brd+2.3193d0
        m1 = m1 * md
        m2 = m2 * md
        m3 = m3 * md
        a1 = a1 * rd
        a2 = a2 * rd
        a3 = a3 * rd
        b2 = (brd * rd)**2
        
        first = .false.
      ENDIF

* compute variables here

      ttphig = ttphig - m1/sqrt(x**2+y**2+(a1+sqrt(z**2+b2))**2)
     &                - m2/sqrt(x**2+y**2+(a2+sqrt(z**2+b2))**2)
     &                - m3/sqrt(x**2+y**2+(a3+sqrt(z**2+b2))**2)
      
      RETURN
      END

************************************************************************
      SUBROUTINE bulgediskhalo(ttphig,x,y,z,t,msc,rsc,tsc,vsc)

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, msc, rsc, tsc, vsc

      integer i
      real*8 mb1, rb1, mb2, rb2, phib
      real*8 md(3), a(3), b2, phid
      real*8 vh2h, r02, phih
      real*8 r2
      save mb1, rb1, mb2, rb2, md, a, b2, vh2h, r02

      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here

*     Bulge parameters
        mb1 = 3.0d9 / msc   ! Msun -> Nbody units
        rb1 = 2700.d0 / rsc  ! pc -> Nbody units
      
        mb2 = 1.6d10 / msc  ! Msun -> Nbody units
        rb2 = 420.d0 / rsc   ! pc -> Nbody units

*     Disk parameters
        md(1) = 8.9d10 / msc   ! Msun -> Nbody units
        md(2) = -6.9d10 / msc  ! Msun -> Nbody units
        md(3) = 2.8d10 / msc   ! Msun -> Nbody units
        a(1) = 5.d3 / rsc       ! pc -> Nbody units
        a(2) = 15.8d3 / rsc    ! pc -> Nbody units
        a(3) = 33.d3 / rsc      ! pc -> Nbody units
        b2 = (300.d0 / rsc)**2  ! pc -> Nbody units

*     Halo parameters
        vh2h = 0.5d0 * (225.d0 / vsc)**2   ! km/s -> Nbody units
        r02 = (8.4d3 / rsc)**2          ! pc -> Nbody units

        first = .false.
      ENDIF

* compute variables here
      r2 = x**2 + y**2 + z**2
      
      phib = -mb1/sqrt(r2 + rb1**2) - mb2/sqrt(r2+rb2**2)
      
      phid = 0.d0
      do i=1,3
        phid = phid - md(i)/sqrt(x**2+y**2+(a(i)+sqrt(z**2+b2))**2)
      enddo

      phih = vh2h * log(r2 + r02)

      ttphig = ttphig + phib + phid + phih
      
      RETURN
      END

************************************************************************
      SUBROUTINE spiralarms(ttphig,x,y,z,t,msc,rsc,tsc,vsc)

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, msc, rsc, tsc, vsc


      integer N
      real*8 pi, alpha, H, Rs, r0, A, omega
      real*8 theta, xp, yp, r, KH, beta, D

      save N, pi, alpha, H, Rs, r0, A, omega


      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here

        pi = 4.0d0*atan(1.d0)
*     Spiral arm parameters for Cox & Gomez (2002), A&A
        N = 2
        alpha = 15.5d0*pi/180.d0
        H = 300.d0 / rsc   ! pc -> Nbody units
        Rs = 2600.d0 / rsc ! pc -> Nbody units
        r0 = 5600.d0 /rsc  ! pc -> Nbody units
        A = 0.0336d0 / msc * rsc**3  ! Msun / pc^3 -> Nbody units
        omega = 20.d0 * 0.001023d0 * tsc ! km/s/kpc -> Myr^-1 -> nbody

        first = .false.
      ENDIF

* compute variables here
      theta = -omega*t
      xp = cos(theta)*x - sin(theta)*y
      yp = sin(theta)*x + cos(theta)*y
      r = sqrt(xp**2 + yp**2 + z**2)

      KH = H*DBLE(N)/(r*sin(alpha))
      beta= KH*(1.d0 + 0.4d0*KH)
      D = (1.d0 + KH + 0.3d0*(KH)**2)/(1.d0+0.3d0*KH)
      
      ttphig = ttphig + ( -4.d0*pi*H*A*exp(-(r-r0)/Rs) *
     &    cos(2.d0*(atan2(yp,xp) - log(r/r0)/tan(alpha))) ) /
     &    (cosh(KH*z/beta)**beta*(KH*D))
      
      RETURN
      END


************************************************************************
      SUBROUTINE nfwcosmo(ttphig,x,y,z,t,msc,rsc,tsc,vsc)
! Buist & Helmi
      implicit none
      logical first
      real*8 x, y, z, t, ttphig, msc, rsc, tsc, vsc
      
      real*8 p1, p2, mgm, rs, expmz, zz, r, h0, omega_m, t0

      save p1, p2, mgm, rs, h0, omega_m, t0

      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here

! age of the universe at t = 0   (Myr -> Nbody units)
        t0 = 1169.16d0 / tsc 

! Omega_matter
        omega_m = 0.31d0
        
! (1 - Omega_m / Omega_m )^(1/3)
        p1 = ((1.d0-omega_m)/omega_m)**(1.d0/3.d0)

! Hubble constant km/s/Mpc -> pc / Myr / Mpc -> 1/Myr -> nbody units
        h0 = 68.d0 * 1.023d0 * 1.0d-6 * tsc

! 2 / 3  * 1/ ( H0 * sqrt(1-Omega_m)  )  : unit of time in nbody units
        p2 = 2.d0/(3.d0 * h0 * sqrt(1.d0-omega_m))

        mgm = 1.5d11 / msc ! 1.5e11 Msun
        rs = 16.0d3 / rsc ! 16 kpc

        first = .false.
      ENDIF

* compute variables here
      r = sqrt(x**2 + y**2 + z**2)


      zz = p1 * (sinh((t+t0)/p2))**(-2.d0/3.d0) -1.d0
!     zz = 0 ! static version

      ! exp(-0.1*z)
      expmz = EXP(-0.1d0*zz)
      
      ttphig = ttphig - mgm * expmz**2 / r * log(1.d0 + r/(rs*expmz))
      
      RETURN
      END


************************************************************************
************************************************************************
************************************************************************
**** Subroutine template
      SUBROUTINE mynewpot(ttphig,x,y,z,t,msc,rsc,tsc,vsc)

      implicit none
      logical first
      real*8 x, y, z, t, ttphig, msc, rsc, tsc, vsc


      save first
      data first /.TRUE./
      IF(first) THEN
* compute constants here


        first = .false.
      ENDIF

* compute variables here

      ttphig = ttphig - 0.d0
      
      RETURN
      END

************************************************************************
************************************************************************
************************************************************************
