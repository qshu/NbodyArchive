        IMPLICIT REAL*8 (A-H,M,O-Z)                                            
        PARAMETER (NMX=200,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,                    !
     &  NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)                                         !
         COMMON/DataForRoutines1/X(NMX3),V(NMX3),WTTL,M(NMX),                  !
     &   XC(NMX3),WC(NMX3),MC(NMX)                                             !
     &  ,XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),N                   !
      COMMON/DataForChainRoutinesTwo/MMIJ,CMX(3),CMV(3)
     & ,ENERGY,Energr,CHTIME          !
      common/softening/ee,cmethod(3),Clight,NofBH                              !
      common/TIMECOMMON/Taika,timecomparison,steppi              
      common/spincommon/spin(3)! the relative spin of M(1) !Spin=spin*G*M^2/c
      common/tolerancecommon/EPS
        common/slowi/aat(nmx),rix(nmx),rsch(nmx),OK(nmx,-1:1)
     &,VC(-2:NMX3+3)
        common/rajahairio/tiny_pertu
        common/Hidastukset/NumberOfSlowdowns
c       The following publications contain the info needed for understanding this code


c \bibitem[Hellstr{\"o}m and Mikkola(2010)]{2010CeMDA.106..143H} Hellstr{\"o}m, C., Mikkola, S.\ 2010.\ Explicit algorithmic regularization in the few-body problem for velocity-dependent perturbations.\ Celestial Mechanics and Dynamical Astronomy 106, 143-156. 

c \bibitem[Mikkola and Merritt(2008)]{2008AJ....135.2398M} Mikkola, S., Merritt, D.\ 2008.\ Implementing Few-Body Algorithmic Regularization with Post-Newtonian Terms.\ The Astronomical Journal 135, 2398-2405. 


c \bibitem[Mikkola and Merritt(2006)]{2006MNRAS.372..219M} Mikkola, S., Merritt, D.\ 2006.\ Algorithmic regularization with velocity-dependent forces.\ Monthly Notices of the Royal Astronomical Society 372, 219-223. 


c \bibitem[Mikkola and Aarseth(2002)]{2002CeMDA..84..343M} Mikkola, S., Aarseth, S.\ 2002.\ A Time-Transformed Leapfrog Scheme.\ Celestial Mechanics and Dynamical Astronomy 84, 343-354. 


c \bibitem[Mikkola and Aarseth(1996)]{1996CeMDA..64..197M} Mikkola, S., Aarseth, S.~J.\ 1996.\ A Slow-down Treatment for Close Binaries.\ Celestial Mechanics and Dynamical Astronomy 64, 197-208. 


c \bibitem[Mikkola and Aarseth(1993)]{1993CeMDA..57..439M} Mikkola, S., Aarseth, S.~J.\ 1993.\ An implementation of N-body chain regularization.\ Celestial Mechanics and Dynamical Astronomy 57, 439-459. 


c \bibitem[Mikkola and Tanikawa(2013)]{2013NewA...20...38M} Mikkola, S., Tanikawa, K.\ 2013.\ Implementation of an efficient logarithmic-Hamiltonian three-body code.\ New Astronomy 20, 38-41. 


c \bibitem[Mikkola and Tanikawa(2013)]{2013MNRAS.430.2822M} Mikkola, S., Tanikawa, K.\ 2013.\ Regularizing dynamical problems with the symplectic logarithmic Hamiltonian leapfrog.\ Monthly Notices of the Royal Astronomical Society 430, 2822-2827. 


c \bibitem[Mikkola and Tanikawa(1999)]{1999MNRAS.310..745M} Mikkola, S., Tanikawa, K.\ 1999.\ Algorithmic regularization of the few-body problem.\ Monthly Notices of the Royal Astronomical Society 310, 745-749. 
