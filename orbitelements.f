       SUBROUTINE ORBITELEMENTS(I,XI,XIDOT)
*
*
*       Drag force & first derivative.
*       -------------------------
*
      INCLUDE 'common6.h'
      REAL*8  RSTARI,RSTARJ,EXCENTR1(NMAX),RVX,RVY,RVZ,VRT1,VRT2,VRT1T2,
     &        INCLIN1(NMAX),ANGLE1,INCLIN_L(NMAX),XOLDN(ID,NMAX),
     &        EXCENTR2(NMAX),TINMASS,CENER,STARNUMBER1,CEXC,PROVR(NMAX),
     &        SENERGY(NMAX),POINT1(NMAX),POINT2(NMAX),ECC_SUPER2(NMAX),
     &        POINT3(NMAX),CONST1,CONST3,R_MAX(NMAX),ECC_SUPER1(NMAX),
     &        R_MIN(NMAX),XJ(NMAX),XJDOT(NMAX),XJN(NMAX),XCD(NMAX),
     &        P_ORBIT,PR1,PR3,R1R3,REAL_MIN(NMAX),REAL_MAX(NMAX),
     &        ANGLE2(NMAX),R_CURR(NMAX),RTEK,XI(3),XIDOT(3),RDENS1,
     &        A_ORBIT(NMAX)
      INTEGER STARNUMBER
      LOGICAL LPR,MAX_CANDIDATE(NMAX),MIN_CANDIDATE(NMAX),APOGEI(NMAX),
     &        PERIGEI(NMAX),FIRST1
      SAVE FIRST1
      DATA FIRST1 /.TRUE./
*
      LPR=NAME(I).EQ.1927
*
      TINMASS=CMBLHOLE+1.0E-10
*       Recalling the Name of the Star
      J=NAME(I)
*       Memorizing previous positions
      IF(ITIMER(J).EQ.0) THEN
         MAX_CANDIDATE(J)=.FALSE.
         MIN_CANDIDATE(J)=.FALSE.
         APOGEI(J)=.FALSE.
         PERIGEI(J)=.FALSE.
*       Ficticious values of Rmax & Rmin
         R_MIN(J)=1000.0D0
         R_MAX(J)=0.0D0
         GOTO 555
      END IF
*
      DO K=1,3
         XJ(K)=XOLD(K,J)
         XJDOT(K)=XDOTOLD(K,J)
         XJN(K)=XOLDN(K,J)-RDENS(K)
         XCD(K)=XI(K)-RDENS(K)
      END DO
*
*       Turning-Points of Orbits
*      VRT1=XCD(1)*XIDOT(1)+XCD(2)*XIDOT(2)+XCD(3)*XIDOT(3)
*      VRT2=XJ(1)*XJDOT(1)+XJ(2)*XJDOT(2)+XJ(3)*XJDOT(3)
      VRT1=XI(1)*XIDOT(1)+XI(2)*XIDOT(2)+XI(3)*XIDOT(3)
      VRT2=XJ(1)*XJDOT(1)+XJ(2)*XJDOT(2)+XJ(3)*XJDOT(3)
      VRT1T2=VRT1*VRT2
*
      R_CURR(J)=SQRT(XCD(1)**2.0+XCD(2)**2.0+XCD(3)**2.0)
      RTEK=SQRT(XI(1)**2.0+XI(2)**2.0+XI(3)**2.0)
      RDENS1=SQRT(RDENS(1)**2.0+RDENS(2)**2.0+RDENS(3)**2.0)
      IF(VRT1T2.LE.0.0D0) THEN
      IF(VRT1.GT.0.0D0 .AND. RTEK.LT.R_MIN(J)) THEN    
          MIN_CANDIDATE(J)=.TRUE.
          R_MIN(J)=RTEK
          REAL_MIN(J)=R_CURR(J)
      IF(LPR) PRINT*,'MIIIIIIIIIIIIIIIIIIIIIIIIIIIIII',R_MIN(J)
          DO K=1,3
          XOLD2(K,J)=XI(K)
          END DO
      ELSEIF(VRT1.LT.0.0D0 .AND. RTEK.GT.R_MAX(J)) THEN
          MAX_CANDIDATE(J)=.TRUE.
          R_MAX(J)=RTEK
          REAL_MAX(J)=R_CURR(J)
      IF(LPR) PRINT*,'MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',R_MAX(J)
          DO K=1,3
          XOLD1(K,J)=XI(K)
          END  DO
      END IF
      ENDIF
*
      IF(MIN_CANDIDATE(J)) THEN
         ANGMIN(J)=ABS(XOLD2(1,J)*XI(1)+XOLD2(2,J)*XI(2)+
     &                 XOLD2(3,J)*XI(3))/(R_MIN(J)*RTEK) 
        IF(ANGMIN(J).LT.0.85) THEN
          IF(R_MIN(J).LT.RTEK) THEN
       IF(LPR) PRINT*,'MREALIIIIIIIIIII',R_MIN(J),ANGMIN(J)
          PERIGEI(J)=.TRUE.
          RSTARMIN(J)=REAL_MIN(J)
          END IF
          MIN_CANDIDATE(J)=.FALSE.
          MAX_CANDIDATE(J)=.FALSE.
          R_MIN(J)=10.0D0
          R_MAX(J)=0.0D0
        END IF
      ENDIF
*
      IF(MAX_CANDIDATE(J)) THEN
         ANGMAX(J)=ABS(XOLD1(1,J)*XI(1)+XOLD1(2,J)*XI(2)+
     &                 XOLD1(3,J)*XI(3))/(R_MAX(J)*RTEK) 
        IF(ANGMAX(J).LT.0.85) THEN
          IF(R_MAX(J).GT.RTEK) THEN
        IF(LPR) PRINT*,'MIIREALMAXXXXXXXXXXXX',R_MAX(J)
          APOGEI(J)=.TRUE.
          RSTARMAX(J)=REAL_MAX(J)
          END IF
          MAX_CANDIDATE(J)=.FALSE.
          MIN_CANDIDATE(J)=.FALSE.
          R_MAX(J)=0.0
          R_MIN(J)=10.0
        END IF
      ENDIF
*
*           Calculating of Mass inside of Sphere of R_MAX(I)      
*      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*      IF(LPR) THEN
*      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(TIME.GT.0.0) THEN
*      IF(MAX_CANDIDATE(J)) THEN
         DO 519 K=IFIRST,NTOT
         IF(K.EQ.J) GO TO 519
         RSTARJ=SQRT((X(1,K)-RDENS(1))**2.0+(X(2,K)-RDENS(2))**2.0+
     &               (X(3,K)-RDENS(3)**2.0))
*        IF(REAL_MAX(J).GE.RSTARJ) TINMASS=TINMASS+BODY(K)
         IF(R_CURR(J).GE.RSTARJ) TINMASS=TINMASS+BODY(K)
 519     CONTINUE
*
*       Laplace Eccentricity 'ECC_L(I)' and Inclination (instant mean) 
      RVX=XCD(1)*XIDOT(3)-XIDOT(1)*XCD(3) 
      RVY=XCD(1)*XIDOT(2)-XIDOT(1)*XCD(2)
      RVZ=XCD(2)*XIDOT(3)-XIDOT(2)*XCD(3)    
      ECC_L(J)=SQRT((((-RVX*XIDOT(3)-RVY*XIDOT(2)))+
     &            XCD(1)*TINMASS/R_CURR(J))**2+
     &               (((RVY*XIDOT(1)-RVZ*XIDOT(3)))+
     &            XCD(2)*TINMASS/R_CURR(J))**2+
     &               (((RVZ*XIDOT(2)+RVX*XIDOT(1)))+
     &            XCD(3)*TINMASS/R_CURR(J))**2)/TINMASS
      ANGLE1= (XCD(1)*XIDOT(2)-XCD(2)*XIDOT(1))/
     &        (SQRT((XCD(2)*XIDOT(3)-XCD(3)*XIDOT(2))**2.0+
     &              (XCD(3)*XIDOT(1)-XCD(1)*XIDOT(3))**2.0+
     &              (XCD(1)*XIDOT(2)-XCD(2)*XIDOT(1))**2.0))
      INCLIN_L(J)=(360/TWOPI)*ACOS(ANGLE1)
*
*       Three-point's method eccentricity 'ECC_3P'
*     IF(APOGEI(J)) THEN
      DO K=1,3
      XJ(K)=XJ(K)-RDENS(K)
      END DO
      CONST1=SQRT((XJ(2)*XCD(3)-XJ(3)*XCD(2))**2.0+
     &            (XJ(3)*XCD(1)-XJ(1)*XCD(3))**2.0+
     &            (XJ(1)*XCD(2)-XJ(2)*XCD(1))**2.0)/
     &       SQRT((XJN(2)*XCD(3)-XJN(3)*XCD(2))**2.0+
     &            (XJN(3)*XCD(1)-XJN(1)*XCD(3))**2.0+
     &            (XJN(1)*XCD(2)-XJN(2)*XCD(1))**2.0)
      CONST3=SQRT((XJ(2)*XJN(3)-XJ(3)*XJN(2))**2.0+
     &            (XJ(3)*XJN(1)-XJ(1)*XJN(3))**2.0+
     &            (XJ(1)*XJN(2)-XJ(2)*XJN(1))**2.0)/
     &       SQRT((XCD(2)*XJN(3)-XCD(3)*XJN(2))**2.0+
     &            (XCD(3)*XJN(1)-XCD(1)*XJN(3))**2.0+
     &            (XCD(1)*XJN(2)-XCD(2)*XJN(1))**2.0)
      POINT1(J)=SQRT(XJN(1)**2.0+XJN(2)**2.0+XJN(3)**2.0) 
      POINT2(J)=SQRT(XJ(1)**2.0+XJ(2)**2.0+XJ(3)**2.0) 
      POINT3(J)=SQRT(XCD(1)**2.0+XCD(2)**2.0+XCD(3)**2.0) 
      P_ORBIT=(CONST1*POINT1(J)+CONST3*POINT3(J)-POINT2(J))/
     &        (CONST1+CONST3-1.0)
      PR1=P_ORBIT-POINT1(J)
      PR3=P_ORBIT-POINT3(J)
      R1R3=SQRT((XJN(2)*XCD(3)-XJN(3)*XCD(2))**2.0+
     &          (XJN(3)*XCD(1)-XJN(1)*XCD(3))**2.0+
     &          (XJN(1)*XCD(2)-XJN(2)*XCD(1))**2.0)
      ECC_3P(J)=(SQRT((PR1*XCD(1)-PR3*XJN(1))**2.0+
     &                (PR1*XCD(2)-PR3*XJN(2))**2.0+
     &                (PR1*XCD(3)-PR3*XJN(3))**2.0))/R1R3
      A_ORBIT(J)=ABS(P_ORBIT)/(1-ECC_3P(J)**2.0)
      INCLIN3(J)=(360.0/TWOPI)*ACOS((XJ(1)*XCD(2)-
     &                             XJ(2)*XCD(1))/
     &        SQRT((XJ(2)*XCD(3)-XJ(3)*XCD(2))**2.0+
     &             (XJ(3)*XCD(1)-XJ(1)*XCD(3))**2.0+
     &             (XJ(1)*XCD(2)-XJ(2)*XCD(1))**2.0))
      INCLIN5(J)=(360.0/TWOPI)*ACOS((XJN(1)*XJ(2)-
     &                               XJN(2)*XJ(1))/
     &        SQRT((XJN(2)*XJ(3)-XJN(3)*XJ(2))**2.0+
     &             (XJN(3)*XJ(1)-XJN(1)*XJ(3))**2.0+
     &             (XJN(1)*XJ(2)-XJN(2)*XJ(1))**2.0))
*      END IF 
      END IF
*
*     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*      END IF
*     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*      Excentricity & Inclination of Orbits
*     IF(MAX_CANDIDATE(J)) THEN
      IF(TIME.GT.0.0) THEN
*     IF(APOGEI(J).AND.PERIGEI(J)) THEN
      APOGEI(J)=.FALSE.
      PERIGEI(J)=.FALSE.
      ECC_PA(J)=ABS(RSTARMAX(J)-RSTARMIN(J))/(RSTARMAX(J)+RSTARMIN(J))
*
      K11=K11+1
      I11=0
      NUMBER1(K11)=J
      IF(K11.GT.1) THEN
      DO K=1,K11-1
      IF(NUMBER1(K).EQ.NUMBER1(K11)) THEN
      I11=1
      END IF
      END DO   
      K11=K11-I11
      END IF
*
*      IF(LPR) THEN
      IF(I11.EQ.0) THEN
      PRINT*,'***************************** a new star:'
      PRINT*,'N,NUMBER(I):',K11,NUMBER1(K11),NAME(I),NTOT
      PRINT*,'INCLIN BY PLANE & ANG. MOMENTUM:',INCLIN_L(J)
      PRINT*,'INCLIN BY 3 POINTS             :',INCLIN3(J),INCLIN5(J)
      WRITE(22,518) K11,NUMBER1(K11),INCLIN3(J),INCLIN_L(J)
 518  FORMAT(2(I5),D15.5)
      INCLDISTR(K11)=INCLIN_L(J)
      INCLDISTR1(K11)=INCLIN3(J)
      EXCEDISTR(K11)=ECC_PA(I)
      EXCEDISTR1(K11)=ECC_L(J)
      EXCEDISTR2(K11)=ECC_3P(J)
      EXCEDISTR3(K11)=EXCENTR222(J)
*
      RSTARI=SQRT(XCD(1)*XCD(1)+XCD(2)*XCD(2)+XCD(3)*XCD(3))        
      EXCENTR222(J)=SQRT(1.0-((RSTARMIN(J)**2.0)/(RSTARMAX(J)**2.0)))
      PRINT*,'ECCENTRCITY BY Rmax & Rmin:',ECC_PA(J)
      PRINT*,'MOMENTARY ECCENTRCITY     :',ECC_L(J)
      PRINT*,'ECCENTR BY 3 POINTS!      :',ECC_3P(J)
      PRINT*,'ECCENTRCITY IN HARMONIC P :',EXCENTR222(J)
      PRINT*,'ECC_SUPER1,2!      :',ECC_SUPER1(J),ECC_SUPER2(J)
      PRINT*,'X_POSITION R1,R2,R3:',XCD(1),XJ(1),XJN(1)
      PRINT*,'Y_POSITION R1,R2,R3:',XCD(2),XJ(2),XJN(2)
      PRINT*,'Z_POSITION R1,R2,R3:',XCD(3),XJ(3),XJN(3),UGOL
      PRINT*,'r1,r2,r3,R1R3        :',POINT1(J),POINT2(J),POINT3(J),R1R3
      PRINT*,'RMN, RMX,RCRR,RDNS:',RSTARMIN(J),RSTARMAX(J),RSTARI,RDENS1
*             SENERGY(I)  
*      END IF
      END IF
*
*
      IF(K11.EQ.5000) THEN
      PRINT*,'######################################'
      PRINT*,'######################################'
      PRINT*,'######################################'
      PRINT*,'######################################'
      PRINT*,'######################################'
      PRINT*,'######################################'
      PRINT*,'######################################'
      PRINT*,'######################################'
      PRINT*,'######################################'
      K11=0
      NFILE=NFILE+1
      I2=0
      DO IANGLE=0,189,1
      STARNUMBER=0
      DO I1=1,N
      IF(INCLDISTR(I1).GT.IANGLE.AND.INCLDISTR(I1).LE.IANGLE+1) THEN
      STARNUMBER=STARNUMBER+1
      I2=I2+1
      END IF
      END DO
      WRITE(3,513) IANGLE,STARNUMBER,I2
 513  FORMAT(3(I5))
      END DO
*
      I2=0
      DO IANGLE=0,180,1
      STARNUMBER=0
      DO I1=1,N
      IF(INCLDISTR1(I1).GT.IANGLE.AND.INCLDISTR1(I1).LE.IANGLE+1) THEN
      STARNUMBER=STARNUMBER+1
      I2=I2+1
      END IF
      END DO
      WRITE(4,524) IANGLE,STARNUMBER,I2
 524  FORMAT(3(I5))
      END DO
*
      DO CEXC=0.0,3.99,0.01
      STARNUMBER=0
      DO I1=1,N
      IF(EXCEDISTR(I1).GT.CEXC.AND.EXCEDISTR(I1).LE.CEXC+0.01) THEN
      STARNUMBER=STARNUMBER+1
      END IF
      END DO
      WRITE(35,520) CEXC,STARNUMBER
 520  FORMAT(D15.5,I5)
      END DO
*
      DO CEXC=0.0,3.99,0.01
      STARNUMBER=0
      DO I1=1,N
      IF(EXCEDISTR1(I1).GT.CEXC.AND.EXCEDISTR1(I1).LE.CEXC+0.01) THEN
      STARNUMBER=STARNUMBER+1
      END IF
      END DO
      WRITE(36,523) CEXC,STARNUMBER
 523  FORMAT(D15.5,I5)
      END DO
*
      DO CEXC=0.0,3.99,0.01
      STARNUMBER=0
      DO I1=1,N
      IF(EXCEDISTR2(I1).GT.CEXC.AND.EXCEDISTR2(I1).LE.CEXC+0.01) THEN
      STARNUMBER=STARNUMBER+1
      END IF
      END DO
      WRITE(39,531) CEXC,STARNUMBER
 531  FORMAT(D15.5,I5)
      END DO
*
      DO CEXC=-2.0,10.99,0.01
      STARNUMBER=0
      DO I1=1,N
      IF(A_ORBIT(I1).GT.CEXC.AND.A_ORBIT(I1).LE.CEXC+0.01) THEN
      STARNUMBER=STARNUMBER+1
      END IF
      END DO
      WRITE(40,532) CEXC,STARNUMBER
 532  FORMAT(D15.5,I5)
      END DO
*
      DO CEXC=0.0,5.0,0.05
      STARNUMBER=0
      DO I1=1,N
      IF(AXISDISTR(I1).GT.CEXC.AND.AXISDISTR(I1).LE.CEXC+0.01) THEN
      STARNUMBER=STARNUMBER+1
      END IF
      END DO
      WRITE(37,526) CEXC,STARNUMBER
 526  FORMAT(D15.5,I5)
      END DO
*
      CALL ENERGY
      DO I1=1,N
      SENERGY(I1)=0.5*BODY(I1)*(XDOT(1,I1)**2+XDOT(2,I1)**2+
     &                          XDOT(3,I1)**2)+PHIDBL(I1)
      PARTENERGY(I1)=SENERGY(I1)
      END DO
*
      DO CENER=-10.0,1.0,0.01
      STARNUMBER=0
      DO I1 =1,N
      IF(PARTENERGY(I1).GT.CENER.AND.PARTENERGY(I1).LE.CENER+0.01) THEN
      STARNUMBER=STARNUMBER+1
      END IF
      END DO
      WRITE(38,522) CENER,STARNUMBER
 522  FORMAT(D15.5,I5)
      END DO
*
      END IF
*
      END IF
*
*//////////////////////////////////////////////////////////////////////////
      IF(LPR) THEN
*      PRINT*,'BBBBBBBBBBBBBBBBBBBBBBBBBBBBB',BODY(J)
      I1=24
      IF(FIRST1) THEN
*      OPEN(UNIT=1,STATUS='NEW',FORM='UNFORMATTED',FILE='ARI1.bin')
*      OPEN(UNIT=2,STATUS='NEW',FORM='UNFORMATTED',FILE='ARI2.bin')
      FIRST1=.FALSE.
      END IF
*
      WRITE (1) ,TIME,(XI(K),K=1,3)
      WRITE (2) ,TIME,INCLIN_L(J),ECC_3P(J),ECC_L(J)
      WRITE(54,511) TIME,RRDENS,(XI(K),XIDOT(K),K=1,3)
      WRITE(56,512) TIME,INCLIN_L(J),ECC_3P(J),ECC_L(J)
 511  FORMAT(1X,1P,13(D15.5,1X))
 512  FORMAT(1X,1P,4(D15.5,1X))
      END IF
 555  ITIMER(J)=1
      DO K=1,3
      XOLDN(K,J)=XOLD(K,J)
      XOLD(K,J)=XI(K)
      XDOTOLD(K,J)=XIDOT(K)
      END DO
*
      RETURN
*
      END
      





























