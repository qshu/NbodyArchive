      SUBROUTINE BHOUT2
*
*
*     Analysis of the BHs influence
*
*
      INCLUDE 'common6.h'
*
      COMMON/BHole/   BHx(3), BHv(3),BHTIME, CURJ
*
      REAL*8 dist2, dist(3,NMAX), vr, vrarr(NMAX), rtidal,
     &      hv(3), rrel(3), mu, rreldot(3), rdot2, r2, Av(3),
     &     ecc, TBHout, Energy, ecosphi
      INTEGER CURJ,Nold, Ntotold
      LOGICAL ISSET
      DATA isset/.FALSE./
      DATA TBHout/0.0D0/
      SAVE
*
      jpair = 0
*
 2    CONTINUE
      IF (CURJ.ge.IFIRST)then
         CALL XVPRED(CURJ,0)
      ELSE
         JPAIR = KVEC(CURJ)
         CALL XVPRED(N+JPAIR,0)
         CALL RESOLV(JPAIR,1)
      end if
      BHx(1)=X(1,CURJ)
      BHx(2)=X(2,CURJ) 
      BHx(3)=X(3,CURJ)
      BHv(1)=XDOT(1,CURJ)
      BHv(2)=XDOT(2,CURJ)
      BHv(3)=XDOT(3,CURJ)
      BHBODY=BODY(CURJ)
      BHTIME=TIME+TOFF
      rn    = BHBODY/(vc**2)
*
*
*      roi=0.02D0
      rbhmin=10*rn
      noi=0
      NNB=LIST(1,CURJ)
      if (jpair.gt.0)  nnb=list(1,2*jpair+1)
*     search for passing stars
*     only search in the neighbor-list:
      DO 10 L=2, NNB
         I=LIST(L,CURJ)
         if (jpair.gt.0) I=LIST(L,2*jpair+1)
         CALL XVPRED(I,0)
         dist2=(BHx(1)-X(1,I))**2+(BHx(2)-X(2,I))**2+(BHx(3)-X(3,I))**2
         vr=0.0D0
            if (sqrt(dist2).lt.rbhmin) rbhmin=sqrt(dist2)
            noi=noi+1
         DO 11 K=1,3
            dist(K,I)=(X(K,I)-BHx(K))/sqrt(dist2)
 11      CONTINUE
         vr=dist(1,I)*XDOT(1,I)+dist(2,I)*XDOT(2,I)+dist(3,I)*XDOT(3,I)
c     Omit check, if index has changed:
         if (name(I).ne.BK(5)
     &        .and.I.ne.KVEC(CURJ)+N.and.Ntot.eq.Ntotold) then
         IF (vr*vrarr(I).lt.0.and.vr.gt.0.0D0)then 
*     define dummy rt for binaries:
            IF (I.lt.N)then
            rtidal = RADIUS(I)*(2*BHBODY/BODY(I))**(1/3.0D0)
            ELSE
               rtidal =RADIUS(2*(I-N)-1)
     &                       +RADIUS(2*(I-N))
            END IF
*
*     Calculate the reduced mass:
      mu=BODY(CURJ)*BODY(I)/(BODY(CURJ)+BODY(I))
*
      rreldot(1)=XDOT(1,I)-BHV(1)
      rreldot(2)=XDOT(2,I)-BHV(2)
      rreldot(3)=XDOT(3,I)-BHV(3)
*
      rdot2=rreldot(1)**2+rreldot(2)**2+rreldot(3)**2
*
      rrel(1)=X(1,I)-BHx(1)
      rrel(2)=X(2,I)-BHx(2)
      rrel(3)=X(3,I)-BHx(3)
*
      vcm1=(bhbody*bhv(1)+body(I)*XDOT(1,I))/(BHBODY+BODY(I))
      vcm2=(bhbody*bhv(2)+body(I)*XDOT(2,I))/(BHBODY+BODY(I))
      vcm3=(bhbody*bhv(3)+body(I)*XDOT(3,I))/(BHBODY+BODY(I))
*
      r2=rrel(1)**2+rrel(2)**2+rrel(3)**2
*
            hv(1)=rrel(2)*rreldot(3)-rrel(3)*rreldot(2)
            hv(2)=rrel(3)*rreldot(1)-rrel(1)*rreldot(3)
            hv(3)=rrel(1)*rreldot(2)-rrel(2)*rreldot(1)
*
*     Calculate the runge-lenz vector:
            Av(1)=(rreldot(2)*hv(3)-rreldot(3)*hv(2))/(BHBODY+BODY(I))-
     &           rrel(1)/sqrt(r2)
            Av(2)=(rreldot(3)*hv(1)-rreldot(1)*hv(3))/(BHBODY+BODY(I))-
     &           rrel(2)/sqrt(r2)
            Av(3)=(rreldot(1)*hv(2)-rreldot(2)*hv(1))/(BHBODY+BODY(I))-
     &           rrel(3)/sqrt(r2)
*
*     Energy in the Kepler-problem:
            Energy=0.5D0*BHBODY*(
     &           (BHV(1)-VCM1)**2+
     &           (BHV(2)-VCM2)**2+
     &           (BHV(3)-VCM3)**2) + 
     &           0.5D0*BODY(I)*(
     &           (XDOT(1,I)-VCM1)**2+
     &           (XDOT(2,I)-VCM2)**2+
     &           (XDOT(3,I)-VCM3)**2)-
     &           BODY(I)*BHBODY/SQRT(R2)
*            Energy=0.5D0*mu*vr**2+mu*(hv(1)**2+hv(2)**2+hv(3)**2)/
*    &           (2.0D0*r2)-mu*(BHBODY+BODY(I))/sqrt(r2)
            ecc=sqrt(Av(1)**2+Av(2)**2+Av(3)**2)
            ecosphi=(Av(1)*rrel(1)+Av(2)*rrel(2)+Av(3)*rrel(3))/sqrt(r2)
*     Pericenter according to Runge-Lenz approach:
            qperi=(hv(1)**2+hv(2)**2+hv(3)**2)/
     &           ((BHBODY+BODY(I))*(1.0D0+ecc))
            ahalf=qperi/(1.0D0-ecc)
            if(rank.eq.0) then
            write(6,301), BHTIME, name(I),
     &           sqrt(dist2), vr, rtidal, ecc, qperi, ahalf, 
     &           ecosphi/ecc, Energy
 301        format (' Pericenter:',1pe15.6,' name=',1i7,' rmin=',
     &           1pe15.6,' vr=',1pe15.6,' rt=',1pe15.6,' ecc=',1pe15.6,
     &           ' qmin=',1pe15.6,' a=',1pe15.6,' cos(phi)=',1d15.6,
     &           ' en=',1pe15.6)
*               CALL FILE_INIT(2)
               WRITE (20) name(I), BHTIME,
     &              sqrt(dist2), rtidal, RADIUS(I), ecc, qperi, ahalf, 
     &              mu**2*(hv(1)**2+hv(2)**2+hv(3)**2), 
     &              Energy, BHBODY, BODY(I), ecosphi/ecc
               CALL FLUSH(20)
            end if
*     Always compare to the magnified rtidal:
            if (dist2.lt.rtidal**2*FLOAT(BK(6))**2)then
               if (rank.eq.0) print*, 'BHout calls disrupt:',
     &              BHTIME,name(I), qperi, 
     &              sqrt(dist2), rtidal
               CALL DISRUPT(I, rtidal)
            end if
         end if
      end if
*     Remember the radial velocity:
      vrarr(I)=vr
 10   CONTINUE


*     Write the BH-diagnostics:
      if (BHTIME.gt.TBHout) then
         if(rank.eq.0) then
*            CALL FILE_INIT(2)
            WRITE (18) BHTIME, (BHx(j), BHv(j),j=1,3),
     &           BHBODY,
     &           rn, rbhmin, noi, rc, (RDENS(j),j=1,3)
            CALL FLUSH(18)
         end if
         TBHout = BHTIME+0.01D0
      end if


      Nold=N
      Ntotold=Ntot
      if (.not.isset)then
         KSTAR(BK(5))=14
         isset=.TRUE.  
      end if
      RETURN
      END
