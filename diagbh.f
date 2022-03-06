c     -----------------------------------------------------------
c     diagbh.F : diagnostics for heavy particles
c     -----------------------------------------------------------
      subroutine printbinary(x1,v1,m1,x2,v2,m2)
      REAL*8 x1(3),v1(3),m1,x2(3),v2(3),m2
      REAL*8 xr(3),vr(3),xc(3),vc(3),
     $     eb,e,a,mt,mu,l2,pot,rp, ke,lx,ly,lz
      integer k
      mt = m1 + m2
      mu = m1*m2/mt
      do 10 k = 1, 3
         xr(k)=x2(k)-x1(k)
         vr(k)=v2(k)-v1(k)
         xc(k)=(m1*x1(k)+m2*x2(k))/mt
         vc(k)=(m1*v1(k)+m2*v2(k))/mt
 10   continue
      pot = -m1*m2/sqrt(xr(1)**2+xr(2)**2+xr(3)**2)
      ke = 0.5*mu*(vr(1)**2+vr(2)**2+vr(3)**2)
      eb = pot + ke
      lx = mu*(xr(3)*vr(2) - xr(2)*vr(3))
      ly = mu*(xr(1)*vr(3) - xr(3)*vr(1))
      lz = mu*(xr(2)*vr(1) - xr(1)*vr(2))
      l2 = lx**2+ly**2+lz**2
c      write(6,*) pot,ke,l2
c     eb = -2mu/a ----
      a = - 0.5*m1*m2/eb
c      write(6,*) a, l2, mu, mt
      if(1-l2/(mu*mu*mt*a) .gt. 0) then
          e = sqrt(1-l2/(mu*mu*mt*a))
      else
          e = 0
      endif
      rp = a*(1-e)
      write(6,600)'B1 :',m1,x1,v1
      write(6,600)'B2 :',m2,x2,v2
 600  format(a5,1p7e16.7)
      write(6,601)eb,a,e,rp
 601  format(' EB,A,E,RP:', 1p4e16.7)
      end
      subroutine outbinary
      include 'common1.h'
      integer k,i
      if(nbh .eq. 2) then
      call printbinary(x(1,1),xdot(1,1),body(1),
     $        x(1,2),xdot(1,2),body(2))
      else
         do  i = 1, nbh
            write(6,600)body(i), (x(k,i),k=1,3),(xdot(k,i),k=1,3)
 600        format(' B :', 1p7e11.3)
         enddo
      endif
      end
      
      
