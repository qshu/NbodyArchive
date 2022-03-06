*       Three-body stability test 
*
*       Rosemary Mardling, School of Mathematical Sciences, Monash University
*       8 September 2003

	real*8 function stability(mm1,mm2,mm3,ein,eout,inc)

        parameter (ne=100)
        implicit real*8 (a-h,m,o-z)
	real*8 inc
        real*8 f22_10(ne),f22_11(ne)
        real*8 f20_10(ne),f20_11(ne)
        real*8 curlyE_22(ne),curlyE_20(ne)

	INCLUDE 'stability.dat'

        if(eout.ge.1.0)then
           xfac=40.0*(1.0+mm3/(mm1+mm2))
           stability=2.8*xfac**0.4
           return
        endif

        e=0.0
        eress=0.0
        Ttau=0.0
	emax=0.995
	pi=4.*atan(1.)

	m1=1
	m2=mm2/mm1
	m3=mm3/mm1

	m12=m1+m2
	m123=m1+m2+m3
	a1=m1/m12
	a2=m2/m12

	Cm=2*pi*2*(3*m1*m2)**0.5/(m12**2*m123)**(1./3.)

	C22=3./8.
	C20=1./4.

	alpha22=(1+cos(inc))**2/4
	alpha20=(2-3*sin(inc)**2)/2

	C2a=sqrt(2*C22*alpha22)*Cm
	C2b=sqrt(2*C20*ABS(alpha20))*Cm

	CE=(m3/m12)**2*(m12/m123)**2/(m1*m2)**2

	CEa=CE*(2*C22)**2*alpha22
	CEb=CE*(2*C20)**2*alpha20

*       l=2, m=2 
*       --------
	n=10
	gg=1
	Ttau1=0

	do j=1,ne
	   eold=e
	   e=j*emax/real(ne)

	   energy=CEa*curlyE_22(j)/(1+ein**2/2)**3+ein**2

	   DeltaTtau1=C2a*energy**0.25*sqrt(abs(f22_10(j)))*n**(1./3.)/
     .         (1+ein**2/2)

	   DeltaTtau2=C2a*energy**0.25*sqrt(abs(f22_11(j)))*(n+1)**(1./3.)/
     .         (1+ein**2/2)

	   Ttauold=Ttau1
	   Ttau1=n+DeltaTtau1
	   Ttau2=(n+1)-DeltaTtau2

	   gold=gg
	   gg=Ttau2-Ttau1

	   if(gg*gold.lt.0)then
	      eress=eold-gold*(e-eold)/(gg-gold)
	      Ttau=Ttauold+(eress-eold)*(Ttau1-Ttauold)/(e-eold)
	      go to 1
	   endif	 
	enddo

1	continue

	eres_a=eress
	Ttau_a=Ttau

	if(eress.lt.0)then
	   eres_a=0
	   Ttau_a=0
	endif

        eres_b=0
        Ttau_b=0
	if(inc.lt.pi/2)go to 5

*       l=2, m=0 
*       --------
	n=10
	gg=1
	Ttau1=0

	do j=1,ne
	   eold=e
	   e=j*emax/real(ne)

	   energy=CEb*curlyE_20(j)/(1+ein**2/2)**3+ein**2

	   DeltaTtau1=C2b*energy**0.25*sqrt(abs(f20_10(j)))*n**(1./3.)/
     .         (1+ein**2/2)

	   DeltaTtau2=C2b*energy**0.25*sqrt(abs(f20_11(j)))*(n+1)**(1./3.)/
     .         (1+ein**2/2)

	   Ttauold=Ttau1
	   Ttau1=n+DeltaTtau1
	   Ttau2=(n+1)-DeltaTtau2

	   gold=gg
	   gg=Ttau2-Ttau1

	   if(gg*gold.lt.0)then
	      eress=eold-gold*(e-eold)/(gg-gold)
	      Ttau=Ttauold+(eress-eold)*(Ttau1-Ttauold)/(e-eold)
	      go to 11
	   endif	 
	enddo

11	continue

	eres_b=eress
	Ttau_b=Ttau

 5	alpha=(1+eres_a)/(1-eres_a)**3
	C_a=Ttau_a/alpha**0.6

        if(inc.lt.pi/2)then
           Ttau_b=0
           C_b=0.0
        else
	   alpha=(1+eres_b)/(1-eres_b)**3
           C_b=Ttau_b/alpha**0.6
        endif

	C=C_a

	eress=eres_a
	Ttau=Ttau_a
	if(C_b.gt.C_a)then
           C=C_b
*          eress=eres_b
*          Ttau=Ttau_b
	endif

	alpha=(1+eout)/(1-eout)**3
	rf=C*alpha**0.6
	Rpcrit=((m123/m12)*rf**2)**0.33333*(1-eout)
*          preliminary care for small ein, retrograde         
        if(inc.gt.pi/2)then
           if(ein.lt.0.5)Rpcrit=2.d0*Rpcrit
        end if
*
*       stability=0.9*Rpcrit
        stability=1.0*Rpcrit

	end
