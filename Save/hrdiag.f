***
      SUBROUTINE hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                  r,lum,kw,mc,rc)
*
*
*       H-R diagram for stars of different metallicity.
*       -----------------------------------------------
*
*       Computes the new mass, luminosity, radius & stellar type.
*       Input (MASS, AJ, TM, TN, LUMS & TSCLS) supplied by routine STAR.
*       Ref: P.P. Eggleton, M.J. Fitchett & C.A. Tout (1989) Ap.J. 347, 998.
*
*       Revised 27th March 1995 by C. A. Tout;
*       24th October 1995 to include metallicity;
*       14th November 1996 to include naked helium stars;
*       28th February 1997 to allow accretion induced supernovae.
*
*       Revised 5th April 1997 by J. R. Hurley
*       to include Z=0.001 as well as Z=0.02, convective overshooting,
*       MS hook and more elaborate CHeB.
*       Final version by JRH 8/98 for Z in range [0.0001,0.03].
*
      implicit none
*
      integer kw,kwp
*
      real*8 mass,aj,mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20)
      real*8 r,lum,mc,rc
      real*8 mch,mxns,mlp
      parameter(mch=1.44d0,mxns=1.8d0,mlp=12.d0)
      real*8 mtc
* 
      real*8 thook,thg,tbagb,tau,tloop,taul,tauh,tau1,tau2,dtau,texp
      real*8 lx,ly,dell,alpha,beta,neta
      real*8 rx,ry,delr,rzams,rtms,gamma,rmin,taumin
      parameter(taumin=5.0d-08)
      real*8 mcmax,mcx,mcy,mcbagb,lambda
      real*8 am,xx,rdgen,mew,lum0,kap,zeta,ahe,aco
      parameter(lum0=7.0d+04,kap=-0.5d0,ahe=6.964d0,aco=44.313d0)
*
      real*8 thookf,tblf
      real*8 lalphf,lbetaf,lnetaf,lhookf,lgbtf,lmcgbf,lzhef,lpertf
      real*8 rzamsf,rtmsf,ralphf,rbetaf,rgammf,rhookf
      real*8 rgbf,rminf,ragbf,rzahbf,rzhef,rhelmf,rpertf
      real*8 mcgbtf,mcgbf,mcheif,mcagbf
      external thookf,tblf
      external lalphf,lbetaf,lnetaf,lhookf,lgbtf,lmcgbf,lzhef,lpertf
      external rzamsf,rtmsf,ralphf,rbetaf,rgammf,rhookf
      external rgbf,rminf,ragbf,rzahbf,rzhef,rhelmf,rpertf
      external mcgbtf,mcgbf,mcheif,mcagbf
*
*
*       ---------------------------------------------------------------------
*       MASS    Stellar mass in solar units (input: old; output: new value).
*       AJ      Current age in Myr.
*       MT      Current mass in solar units (used for R).
*       TM      Main sequence time.
*       TN      Nuclear burning time.
*       TSCLS   Time scale for different stages.
*       LUMS    Characteristic luminosity.
*       GB      Giant Branch parameters
*       ZPARS   Parameters for distinguishing various mass intervals.
*       R       Stellar radius in solar units.
*       TE      Effective temperature (suppressed).
*       KW      Classification type (1 - 15).
*       MC      Core mass.
*       ---------------------------------------------------------------------
*
*
      mc = 0.d0
      tbagb = tscls(2) + tscls(3)
      thg = tscls(1) - tm
*
* Make evolutionary changes to stars that have not reached KW > 5.
*
      if(kw.gt.6) goto 90
*
      if(aj.lt.tscls(1))then
*
*        Either on MS or HG
*
         rzams = rzamsf(mass)
         rtms = rtmsf(mass)
*
         if(aj.lt.tm)then
*
*           Main sequence star.
*
            tau = aj/tm
            thook = thookf(mass)*tscls(1)
            zeta = 0.01d0
            tau1 = MIN(1.d0,aj/thook)
            tau2 = MAX(0.d0,
     &             MIN(1.d0,(aj-(1.d0-zeta)*thook)/(zeta*thook)))
*
            dell = lhookf(mass,zpars(1))
            dtau = tau1**2 - tau2**2
            alpha = lalphf(mass)
            beta = lbetaf(mass)
            neta = lnetaf(mass)
            lx = LOG10(lums(2)/lums(1))
            if(tau.gt.taumin)then
               xx = alpha*tau + beta*tau**neta +
     &              (lx - alpha - beta)*tau**2 - dell*dtau
            else
               xx = alpha*tau + (lx - alpha)*tau**2 - dell*dtau
            endif
            lum = lums(1)*10.d0**xx
*
            delr = rhookf(mass,zpars(1))
            dtau = tau1**3 - tau2**3
            alpha = ralphf(mass)
            beta = rbetaf(mass)
            gamma = rgammf(mass)
            rx = LOG10(rtms/rzams)
* Note that the use of taumin is a slightly pedantic attempt to
* avoid floating point underflow. It IS overkill!
            if(tau.gt.taumin)then
               xx = alpha*tau + beta*tau**10 + gamma*tau**40 +
     &              (rx - alpha - beta - gamma)*tau**3 - delr*dtau
            else
               xx = alpha*tau + (rx - alpha)*tau**3 - delr*dtau
            endif
            r = rzams*10.d0**xx
*
            if(mass.lt.(zpars(1)-0.3))then
               kw = 0
* This following is given by Chris for low mass MS stars which will be 
* substantially degenerate. We need the Hydrogen abundance, X, which we 
* calculate from Z assuming that the helium abundance, Y, is calculated
* according to Y = 0.24 + 2*Z
               rdgen = 0.0258d0*((1.d0+zpars(11))**(5.d0/3.d0))*
     &                          (mass**(-1.d0/3.d0))
               r = MAX(rdgen,r)
            else
               kw = 1
            endif
*
         else 
*
*           Star is on the HG
*
            if(mass.le.zpars(2))then
               mc = mcgbf(lums(3),GB,lums(6))
            elseif(mass.le.zpars(3))then
               mc = mcheif(mass,zpars(2),zpars(9))
            else
               mc = mcheif(mass,zpars(2),zpars(10))
            endif
*
* Test whether core mass has reached total mass.
*
            if(mc.ge.mt)then
               aj = 0.d0
               if(mass.gt.zpars(2))then
*
* Zero-age helium star
*
                  mc = mt
                  mass = mt
                  kw = 7
                  CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               else
*
* Zero-age helium white dwarf.
*
                  mc = mt
                  kw = 10
               endif
            else
               tau = (aj - tm)/thg
               lum = lums(2)*(lums(3)/lums(2))**tau
               if(mass.le.zpars(3))then
                  rx = rgbf(mt,lums(3))
               else
* He-ignition and end of HG occur at Rmin
                  rmin = rminf(mass)
                  ry = ragbf(mt,lums(4),zpars(2))
                  rx = MIN(rmin,ry)
                  if(mass.le.mlp)then
                     texp = log(mass/mlp)/log(zpars(3)/mlp)
                     rx = rgbf(mt,lums(4))
                     rx = rmin*(rx/rmin)**texp
                  endif
                  tau2 = tblf(mass,zpars(2),zpars(3))
                  if(tau2.eq.0.0) rx = ry
               endif
               r = rtms*(rx/rtms)**tau
               kw = 2
            endif
*
         endif
*
* Now the GB, CHeB and AGB evolution.
*
      elseif(aj.lt.tscls(2))then
*
*        Red Giant.
*
         kw = 3
         lum = lgbtf(aj,GB(1),GB,tscls(4),tscls(5),tscls(6))
         if(mass.le.zpars(2))then
* Star has a degenerate He core which grows on the GB
            mc = mcgbf(lum,GB,lums(6))
         else
* Star has a non-degenarate He core which may grow, but
* only slightly, on the GB
            tau = (aj - tscls(1))/(tscls(2) - tscls(1))
            mcx = mcheif(mass,zpars(2),zpars(9))
            mcy = mcheif(mass,zpars(2),zpars(10))
            mc = mcx + (mcy - mcx)*tau
         endif
         r = rgbf(mt,lum)
         if(mc.ge.mt)then
            aj = 0.d0
            if(mass.gt.zpars(2))then
*
* Zero-age helium star
*
               mc = mt
               mass = mt
               kw = 7
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
            else
*
* Zero-age helium white dwarf.
*
               mc = mt
               kw = 10
            endif
         endif
*
      elseif(aj.lt.tbagb)then
*
*       Core helium burning star.
*
         if(kw.eq.3.and.mass.le.zpars(2))then
            mass = mt
            CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
            aj = tscls(2)
         endif
         kw = 4
         if(mass.le.zpars(2))then
            mcx = mcgbf(lums(4),GB,lums(6))
         else
            mcx = mcheif(mass,zpars(2),zpars(10))
         end if
         tau = (aj - tscls(2))/tscls(3)
         mc = mcx + (mcagbf(mass) - mcx)*tau
         if(mass.le.zpars(2))then
            lx = lums(5)
            ly = lums(7)
            rx = rzahbf(mt,mc,zpars(2))
            rmin = rgbf(mt,lx)*zpars(13)**(mass/zpars(2))
            texp = MIN(MAX(0.4d0,rmin/rx),2.5d0)
            ry = ragbf(mt,ly,zpars(2))
            if(rmin.lt.rx)then
               taul = (log(rx/rmin))**(1.d0/3.d0)
            else
               rmin = rx
               taul = 0.d0
            endif
            tauh = (log(ry/rmin))**(1.d0/3.d0)
            tau2 = taul*(tau - 1.0d0) + tauh*tau
            r = rmin*exp(abs(tau2)**3)
            lum = lx*(ly/lx)**(tau**texp)
         elseif(mass.gt.zpars(3))then
*
* For HM stars He-ignition takes place at Rmin in the HG, and CHeB
* consists of a blue phase (before tloop) and a RG phase (after tloop).
*
            tau2 = tblf(mass,zpars(2),zpars(3))
            tloop = tscls(2) + tau2*tscls(3)
            rmin = rminf(mass)
            rx = ragbf(mt,lums(4),zpars(2))
            rmin = MIN(rmin, rx)
            if(mass.le.mlp) then
               texp = log(mass/mlp)/log(zpars(3)/mlp)
               rx = rgbf(mt,lums(4))
               rx = rmin*(rx/rmin)**texp
            else
               rx = rmin
            end if
            texp = MIN(MAX(0.4d0,rmin/rx),2.5d0)
            lum = lums(4)*(lums(7)/lums(4))**(tau**texp)
            if(aj.lt.tloop)then
               ly = lums(4)*(lums(7)/lums(4))**(tau2**texp)
               ry = ragbf(mt,ly,zpars(2))
               taul = 0.d0
               if(rmin.ne.rx) taul = (log(rx/rmin))**(1.d0/3.d0)
               tauh = (log(ry/rmin))**(1.d0/3.d0)
               tau = (aj - tscls(2))/(tau2*tscls(3))
               tau2 = taul*(tau - 1.0d0) + tauh*tau
               r = rmin*exp(abs(tau2)**3)
            else
               r = ragbf(mt,lum,zpars(2))
            end if
         else
*
* For IM stars CHeB consists of a RG phase (before tloop) and a blue
* loop (after tloop).
*
            tau2 = 1.d0 - tblf(mass,zpars(2),zpars(3))
            tloop = tscls(2) + tau2*tscls(3)
            if(aj.lt.tloop)then
               tau = (tloop - aj)/(tau2*tscls(3))
               lum = lums(5)*(lums(4)/lums(5))**(tau**3)
               r = rgbf(mt,lum)
            else
               lx = lums(5)
               ly = lums(7)
               rx = rgbf(mt,lx)
               rmin = rminf(mt)
               texp = MIN(MAX(0.4d0,rmin/rx),2.5d0)
               ry = ragbf(mt,ly,zpars(2))
               if(rmin.lt.rx)then
                  taul = (log(rx/rmin))**(1.d0/3.d0)
               else
                  rmin = rx
                  taul = 0.d0
               endif
               tauh = (log(ry/rmin))**(1.d0/3.d0)
               tau = (aj - tloop)/(tscls(3) - (tloop - tscls(2)))
               tau2 = taul*(tau - 1.d0) + tauh*tau
               r = rmin*exp(abs(tau2)**3)
               lum = lx*(ly/lx)**(tau**texp)
            endif
         endif
* 
* Test whether core mass exceeds total mass.
*
         if(mc.ge.mt)then
*
* Evolved MS naked helium star.
*
            kw = 7
            xx = (aj - tscls(2))/tscls(3)
            mass = mt
            CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
            aj = xx*tm
         endif
      else
*
*        Asymptotic Red Giant.
*
* Set the core mass and then test whether it exceeds either
* the total mass or the maximum allowed core mass.
*
* On the AGB the He core mass remains constant until at Ltp it
* is caught by the C core mass and they grow together.
*
         mcbagb = mcagbf(mass)
         mcx = mcgbtf(tbagb,GB(8),GB,tscls(7),tscls(8),tscls(9))
         mcmax = MAX(MAX(mch,0.773d0*mcbagb-0.35d0),1.05d0*mcx)
         if(aj.lt.tscls(13))then
            mcx = mcgbtf(aj,GB(8),GB,tscls(7),tscls(8),tscls(9))
            mc = mcbagb
            lum = lmcgbf(mcx,GB)
            if(mt.le.mc)then
*
* Evolved naked helium star as the envelope is lost but the
* star has not completed its interior burning. The star becomes
* a post-HeMS star.
*
               kw = 9
               mt = mc
               mass = mt
               mc = mcx
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               if(mc.le.GB(7))then
                  aj = tscls(4) - (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &                            (mc**(1.d0-GB(5)))
               else
                  aj = tscls(5) - (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &                            (mc**(1.d0-GB(6)))
               endif
               aj = MAX(aj,tm)
               goto 90
            else
               kw = 5
            endif
         else
            kw = 6
            mc = mcgbtf(aj,GB(2),GB,tscls(10),tscls(11),tscls(12))
            lum = lmcgbf(mc,GB)
*
* Approximate 3rd Dredge-up on AGB by limiting Mc.
*
            lambda = MIN(0.9d0,0.3d0+0.001d0*mass**5)
            tau = tscls(13)
            mcx = mcgbtf(tau,GB(2),GB,tscls(10),tscls(11),tscls(12))
            mcy = mc
            mc = mc - lambda*(mcy-mcx)
            mcx = mc
            mcmax = MIN(mt,mcmax)   
         endif
         r = ragbf(mt,lum,zpars(2))
*
* Mc,x represents the C core mass
*
         if(mcx.ge.mcmax)then
            aj = 0.d0
            mc = mcmax
            if(mc.lt.mch)then
               if(mcbagb.lt.1.6)then
*     
* Zero-age Carbon/Oxygen White Dwarf
*
                  mt = mc
                  kw = 11
               else
*     
* Zero-age Oxygen/Neon White Dwarf
*
                  mt = mc
                  kw = 12
               endif
            else
               if(mcbagb.lt.1.6)then
*
* Star is not massive enough to ignite C burning.
* so no remnant is left after the SN
*
                  mt = 0.d0
                  r = 0.d0
                  kw = 15
               else
                  mt = 1.17d0 + mc/11.12d0
                  mc = mt
                  if(mt.le.mxns)then
*
* Zero-age Neutron star
*
                     kw = 13
                  else
*
* Zero-age Black hole
*
                     kw = 14
                  endif  
               endif
            endif
         endif
*
      endif
*
 90   continue
*
      if(kw.ge.7.and.kw.le.9)then
*
* Naked Helium Star
*
         rx = rzhef(mt)
         if(aj.lt.tm)then
*
* Main Sequence
*
            kw = 7
            tau = aj/tm
            am = MAX(0.d0,0.85d0-0.08d0*mass)
            lum = lums(1)*(1.d0+0.45d0*tau+am*tau**2)
            am = MAX(0.d0,0.4d0-0.22d0*LOG10(mt))
            r = rx*(1.d0+am*(tau-tau**6))
* Star has no core mass and hence no memory of its past
* which is why we subject mass and mt to mass loss for
* this phase.
            mc = 0.d0
            if(mt.lt.zpars(10)) kw = 10
         else
*
* Helium Shell Burning
*
            kw = 9
            lum = lgbtf(aj,GB(8),GB,tscls(4),tscls(5),tscls(6))
            r = rhelmf(kw,mt,lum,rx,lums(2))
            mc = mcgbf(lum,GB,lums(6))
            mtc = MIN(mt,1.45d0*mt-0.31d0)
            mcmax = MIN(mtc,MAX(mch,0.773d0*mass-0.35d0))
            if(mc.ge.mcmax)then
               aj = 0.d0
               mc = mcmax
               if(mc.lt.mch)then
                  if(mass.lt.1.6)then
*     
* Zero-age Carbon/Oxygen White Dwarf
*
                     mt = MAX(mc,(mc+0.31d0)/1.45d0)
                     kw = 11
                  else
*     
* Zero-age Oxygen/Neon White Dwarf
*
                     mt = mc
                     kw = 12
                  endif
               else
                  if(mass.lt.1.6)then
*
* Star is not massive enough to ignite C burning.
* so no remnant is left after the SN
*
                     mt = 0.d0
                     r = 0.d0
                     kw = 15
                  else
                     mt = 1.17d0 + mc/11.12d0
                     mc = mt
                     if(mt.le.mxns)then
*
* Zero-age Neutron star
*
                        kw = 13
                     else
*
* Zero-age Black hole
*
                        kw = 14
                     endif
                  endif  
               endif
            endif
         endif
      endif
*
      if(kw.ge.10.and.kw.le.12)then
*
*        White dwarf.
*
         mc = mt
         if(mc.ge.mch)then
*
* Accretion induced supernova with no remnant
*
            kw = 15
            aj = 0.d0
            mt = 0.d0
            r = 0.d0
         else
            if(kw.eq.10)then
               xx = zpars(14)/ahe
            else
               xx = zpars(14)/aco
            endif
            lum = 635.d0*mt*xx/(aj+0.1d0)**1.4
            r = 0.0115d0*SQRT(MAX(7.5614367d-07,(mch/mt)**(2.d0/3.d0)
     &                                      - (mt/mch)**(2.d0/3.d0)))
         endif
      endif
*
      if(kw.eq.13)then
*
*        Neutron Star.
*
         mc = mt
         if(mc.gt.mxns)then
*
* Accretion induced Black Hole?
*
            kw = 14
            aj = 0.d0
         else
            lum = 40.d0/(MAX(aj,1.d0))**3.d0
            r = 1.0d-05
         endif
      endif
*
      if(kw.eq.14)then
*
*        Black hole
*
         mc = mt
         lum = 1.0d-06/mt
         r = 4.24d-06*mt
      endif
*
* Calculate the core radius and the luminosity and radius of the
* remnant that the star will become.
*
      tau = 0.d0
      if(kw.le.1.or.kw.eq.7)then
         rc = 0.d0
      elseif(kw.le.3)then
         if(mass.gt.zpars(2))then
            lx = lzhef(mc)
            rx = rzhef(mc)
            rc = rx
         else
            lx = 635.d0*mc*zpars(14)/(ahe*0.1d0**1.4)
            rx = 0.0115d0*SQRT(MAX(7.5614367d-07,
     &           (mch/mc)**(2.d0/3.d0)-(mc/mch)**(2.d0/3.d0)))
            rc = 5.d0*rx
         endif
      elseif(kw.eq.4)then
         tau = (aj - tscls(2))/tscls(3)
         kwp = 7
         CALL star(kwp,mc,mc,tm,tn,tscls,lums,GB,zpars)
         am = MAX(0.d0,0.85d0-0.08d0*mc)
         lx = lums(1)*(1.d0+0.45d0*tau+am*tau**2)
         rx = rzhef(mc)
         am = MAX(0.d0,0.4d0-0.22d0*LOG10(mc))
         rx = rx*(1.d0+am*(tau-tau**6))
         CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
         rc = rx
      elseif(kw.eq.5)then
         kwp = 9
         tau = 3.d0*(aj-tbagb)/(tn-tbagb)
         CALL star(kwp,mc,mc,tm,tn,tscls,lums,GB,zpars)
         lx = lmcgbf(mcx,GB)
         if(tau.lt.1.0) lx = lums(2)*(lx/lums(2))**tau
         rx = rzhef(mc)
         rx = rhelmf(kwp,mc,lx,rx,lums(2))
         CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
         rc = rx
      elseif(kw.le.9)then
         lx = 635.d0*mc*zpars(14)/(aco*0.1d0**1.4)
         rx = 0.0115d0*SQRT(MAX(7.5614367d-07,
     &        (mch/mc)**(2.d0/3.d0) - (mc/mch)**(2.d0/3.d0)))
         rc = 5.d0*rx
      else
         rc = r
      endif
*
* Perturb the luminosity and radius due to small envelope mass.
*
      if(kw.ge.2.and.kw.le.9.and.kw.ne.7)then
         mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
         if(kw.ge.8) mew = ((mtc-mc)/mtc)*5.d0
         if(mew.lt.1.0)then
            xx = lpertf(mt,mew)
            lum = lx*(lum/lx)**xx
            if(r.le.rx)then
               xx = 0.d0
            else
               xx = rpertf(mt,mew,r,rx)
            endif
            r = rx*(r/rx)**xx
         endif
         rc = MIN(rc,r)
      endif
*
      return
      end
***
