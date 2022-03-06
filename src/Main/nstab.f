c     general three-body stability algorithm
c     
c     system is unstable if nstab=1 is returned
c     system is stable if nstab=0 is returned
c     
c     Rosemary Mardling
c     School of Mathematical Sciences, Monash University
c     
c     version as of 16-11-07
c     email mardling@sci.monash.edu.au to be added to updates email list
c     preprint on astro-ph by New Year :-)
c     
c     sigma=period ratio (outer/inner) **should be > 1**
c     ai=inner semi-major axis
c     a0=outer semi-major axis
c     ei0=initial inner eccentricity
c     eo=outer eccentricity
c     relinc=relative inclination (radians)
c     m1, m2, m3=masses (any units; m3=outer body)
c     
c     valid for all inclinations
c     
c     MASS RATIO CONDITIONS
c     valid for systems with at least one of  m_2/m_1>0.05  OR  m_3/m_1>0.05
c     (so that one could have, for example,  m_2/m_1=0  and  m_3/m_1=0.1) 
c     OR BOTH m2/m1>0.01  AND  m3/m1>0.01
c     **future version will include other resonances to cover smaller mass ratios
c     
c     assumes resonance angle phi=0 because resonance overlap criterion doesn't recognize
c     instability outside separatrix.
c     
c     system is unstable if nstab=1 is returned
c     system is stable if nstab=0 is returned

      integer function nstab(ai,a0,ei0,eo,relinc,m1,m2,m3)

      implicit real*8 (a-h,m,o-z)
      common/params2/mm1,mm2,mm3
      common/savepi/pi
      save itime
      data itime/0/

c     reject outer pericentre inside inner apocentre
      if(a0*(1.d0-eo).lt.ai*(1.d0+ei0))then
         nstab=1
         return
      endif

      if(itime.eq.0)then
         pi=4.d0*datan(1.0d0)
         itime=1
      endif

      mm1=m1
      mm2=m2
      mm3=m3
      
      m12=m1+m2
      m123=m12+m3

c     set period ratio (outer/inner)
      a0ai=a0/ai
      sigma=sqrt(a0ai**3*m12/m123)
c     do not allow period ratio < 1
      if(sigma.lt.1.d0)then
         nstab=1
         return
      endif
      
      Mi2=m3/m123
      Mo2=(m1*m2/m12**2)*(m12/m123)**(2.d0/3.d0)
      Mi3=(m3/m12)*(m12/m123)**(4.d0/3.d0)*(m1-m2)/m12
      Mo3=(m1*m2/m12**2)*(m12/m123)*(m1-m2)/m12
      
      c22=3.d0/8.d0
      c20=0.25d0
      c31=sqrt(3.d0)/4.d0
      c33=-sqrt(5.d0)/4.d0
      
      e=eo
      
c     inclination coefficients

      win=0
      
      A=sqrt(1.d0-ei0**2)*cos(relinc)
      Z=(1.d0-ei0**2)*(1.d0+sin(relinc)**2)+5.d0*ei0**2*
     &     (sin(win)*sin(relinc))**2
      Del=z**2+25.d0+16.d0*A**4-10.d0*Z-20.d0*A**2-8.d0*A**2*Z
      
      eK=sqrt(abs((Z+1.d0-4.d0*A**2+sqrt(Del))/6.d0))
      cosIK=A/sqrt(1.d0-eK**2)
      sinIK=sqrt(1.d0-cosIK**2)
      
      gam222=0.25d0*(1.d0+cosIK)**2
      gam22m2=0.25d0*(1.d0-cosIK)**2
      gam220=0.5d0*sqrt(1.5d0)*sinIK**2
      gam200=0.5d0*(3.d0*cosIK**2-1.d0)
           
c     induced inner eccentricity
           ei=ein_induced(sigma,ei0,e,relinc)
           
c     octopole emax	   	   
           if(ABS(m1-m2)>1.D-14)then
              eoctmax=eoct(sigma,ei0,e)
              ei=max(eoctmax,ei)
           endif
           
           ei=max(eK,ei)
           ei=min(ei,1.0d0)
           
           n=sigma
           nstab=0
           
c     [n:1](222) resonance
           s221=-3.d0*ei+(13.d0/8.d0)*ei**3+(5.d0/192.d0)*ei**5         

           f22n=flmn(2,2,n,e)/(1.d0-e)**3

           An=abs(6.d0*c22*s221*f22n*(Mi2+Mo2*sigma**0.666d0)*gam222)
           phi=0.d0
           En=0.5d0*(sigma-n)**2-An*(1.d0+cos(phi))

c     [n+1:1](222) resonance	   
           f22n=flmn(2,2,n+1,e)/(1.d0-e)**3

           An=abs(6.d0*c22*s221*f22n*(Mi2+Mo2*sigma**0.666d0)*gam222)
           
           Enp1=0.5d0*(sigma-(n+1))**2-An*(1.d0+cos(phi))
           if((En.lt.0.d0).and.(Enp1.lt.0.d0)) nstab=1
           
c     [n:1](22-2) resonance
           s22m1=-(ei**3*(4480.d0+1880.d0*ei**2+1091.d0*ei**4))/15360.d0
           f22n=flmn(2,2,n,e)/(1.d0-e)**3

           An=abs(6.d0*c22*s22m1*f22n*(Mi2+Mo2*sigma**0.666d0)*gam22m2)
           phi=0.d0
           En=0.5d0*(sigma-n)**2-An*(1.d0+cos(phi))

c     [n+1:1](22-2) resonance	   
           f22n=flmn(2,2,n+1,e)/(1.d0-e)**3

           An=abs(6.d0*c22*s22m1*f22n*(Mi2+Mo2*sigma**0.666d0)*gam22m2)
           
           Enp1=0.5d0*(sigma-(n+1))**2-An*(1.d0+cos(phi))
           if((En.lt.0.d0).and.(Enp1.lt.0.d0)) nstab=1

c     [n:1](202) resonance
           s201=(ei*(-9216.d0+1152.d0*ei**2-48.d0*ei**4+ei**6))/9216.d0
           f22n=flmn(2,2,n,e)/(1.d0-e)**3

           An=abs(6.d0*sqrt(c20*c22)*s201*f22n*
     &          (Mi2+Mo2*sigma**0.666d0)*gam220)
           
           phi=0.d0
           En=0.5d0*(sigma-n)**2-An*(1.d0+cos(phi))

c     [n+1:1](202) resonance	   
           f22n=flmn(2,2,n+1,e)/(1.d0-e)**3

           An=abs(6.d0*sqrt(c20*c22)*s201*f22n
     &          *(Mi2+Mo2*sigma**0.666d0)*gam220)
           Enp1=0.5d0*(sigma-(n+1))**2-An*(1.d0+cos(phi))
           if((En.lt.0.d0).and.(Enp1.lt.0.d0)) nstab=1

c     [n:1](002) resonance
           s201=(ei*(-9216.d0+1152.d0*ei**2-48.d0*ei**4+ei**6))/9216.d0
           f20n=flmn(2,0,n,e)/(1.d0-e)**3

           An=abs(3.d0*c20*s201*f20n*(Mi2+Mo2*sigma**0.666d0)*gam200)
           
           phi=0.d0
           En=0.5d0*(sigma-n)**2-An*(1.d0+cos(phi))

c     [n+1:1](002) resonance	   
           f20n=flmn(2,0,n+1,e)/(1.d0-e)**3

           An=abs(3.d0*c20*s201*f20n*(Mi2+Mo2*sigma**0.666d0)*gam200)
           
           Enp1=0.5d0*(sigma-(n+1))**2-An*(1.d0+cos(phi))
           if((En.lt.0.d0).and.(Enp1.lt.0.d0)) nstab=1
           
           end

c     -----------------------------------------
c     Asymptotic expression for f^(lm)_n(e) for all e<1 and n.
c     
      real*8 function flmn(l,m,n,e)

      implicit real*8 (a-h,o-z)
      common/savepi/pi

      if(e.lt.5.d-3)then
         if(m.eq.n)then
            flmn=1
         else
            flmn=0
         endif
         return
      endif
      
      rho=n*(1.d0-e)**1.5d0
      
      xi=(a_cosh(1.d0/e)-sqrt(1.d0-e**2))/(1.d0-e)**1.5d0
      
      flmn=(1.d0/(2.d0*pi*n))*2.d0**m*(sqrt(2.d0*pi)/facfac(l,m))*
     .     ((1.d0+e)**(DBLE(3.d0*m-l-1.d0)/4.d0)/e**m)*
     .     (rho**(DBLE(l+m+1.d0)/2.d0))*
     .     exp(-rho*xi)

      end


c     -----------------------------------------
      real*8 function ein_induced(sigma,ei0,e,relinc)

      implicit real*8 (a-h,m,o-z)
      common/params2/m1,m2,m3
      common/savepi/pi

      m123=m1+m2+m3
      n=sigma

      gam222=0.25d0*(1.d0+cos(relinc))**2
      gam220=0.5d0*sqrt(1.5d0)*sin(relinc)**2
      gam200=0.5d0*(3.d0*cos(relinc)**2-1.d0)
                
      f22n=flmn(2,2,n,e)/(1.d0-e)**3
      f20n=flmn(2,0,n,e)/(1.d0-e)**3
                
      prod222=f22n*gam222
      prod220=f22n*gam220
      prod200=f20n*gam200
                          
      prod=max(prod222,prod220,prod200)
                               
      a=4.5d0*(m3/m123)*(2.d0*pi*n)*prod/sigma**2
                               
      ein_induced=sqrt(ei0**2+a**2)
                               
      end
c     -----------------------------------------
c     eoct.f
c     
c     calculates maximum eccentricity for arbitrary coplanar system
c     using Mardling (2007) MNRAS in press
c     
      real*8 function eoct(sigma,ei0,eo)
      implicit real*8 (a-h,m,o-z)
      common/params2/m1,m2,m3
      common/savepi/pi

      m12=m1+m2
      m123=m12+m3
      aoai=((m123/m12)*sigma**2)**0.3333d0
      al=1.d0/aoai
      
      epso=sqrt(1.d0-eo**2)
      
      eeq=1.25d0*al*eo/epso**2/abs(1.d0-sqrt(al)*(m2/m3)/epso)

      AA=abs(1.d0-ei0/eeq)  
      
      if(AA.lt.1.d0)then
         eoct=(1.d0+AA)*eeq
      else
         eoct=ei0+2.d0*eeq
      endif
      
      end

c     -----------------------------------------
      real*8 function a_cosh(x)
      real*8 x

      a_cosh=log(x+dsqrt(x**2-1.d0))

      end
c     -----------------------------------------
      real*8 function sgn(x)
      real*8 x

      if(x.lt.0.d0)then
         sgn=-1.d0
      else
         sgn=1.d0
      endif

      end
c     -----------------------------------------
      real*8 function facfac(l,m)
      implicit real*8 (a-h,o-z)

      prod=1.d0

      n=l+m-1

      do i=1,n,2
         prod=prod*i
      enddo

      facfac=prod

      end
