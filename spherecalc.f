      SUBROUTINE CUBECALC
*
*     Calculate effects of BH motion
*
      INCLUDE 'common6.h'
*
*     Definig Variables needed
      REAL*8 TFINI,OBTIME,OUTSTB,ZAVZB,TIMESB
      REAL*8 CUX,CUY,CUZ,POX,POY,POZ
      REAL*8 DENS
      INTEGER COUNT,INDIC
*
*     How long should the cube be observed?
      TFINI=2.0
*
*     Time when a cube should be observed
      OBSRVA=1.5
      OBSRVB=4.5
      OBSRVC=7.5
      OBSRVD=10.5
      OBSRVE=13.5
      OBSRVF=16.5
*
      IF(TIME.LT.OBSRVA.OR.TIME.GT.(OBSRVA+TFINI))THEN
        IF(TIME.LT.OBSRVB.OR.TIME.GT.(OBSRVB+TFINI))THEN
          IF(TIME.LT.OBSRVC.OR.TIME.GT.(OBSRVC+TFINI))THEN
            IF(TIME.LT.OBSRVD.OR.TIME.GT.(OBSRVD+TFINI))THEN
              IF(TIME.LT.OBSRVE.OR.TIME.GT.(OBSRVE+TFINI))THEN
                IF(TIME.LT.OBSRVF.OR.TIME.GT.(OBSRVF+TFINI))THEN
                  TIMESB=0.0D0
                  RETURN
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
*
*     What output step?
      OUTSTB=0.001 
*
*     Define sphere size
      
*
*     Save cube position
      IF(TIMESB.EQ.0.D0)THEN
         DO 10 J=1,NTOT
         IF(NAME(J).EQ.(NZERO-1))THEN
           POX=X(1,J)
           POY=X(2,J)
           POZ=X(3,J)
           PRINT*,''
           WRITE(*,1) CUX,CUY,CUZ
    1      FORMAT(' Cube size',F15.7,F15.7,F15.7)
         ENDIF
   10 CONTINUE
      ENDIF
*
*     Look whether new output is necessary
      ZAVZB=TIME-TIMESB
      IF(ZAVZB.LT.OUTSTB)THEN
         RETURN
      ENDIF
      TIMESB=TIME
*
*
*     Is first BH still in cube?
      DO 20 J=1,NTOT
        IF(NAME(J).EQ.(NZERO-1))THEN
          IF(X(1,J).GT.(POX-CUX/2) .AND. X(1,J).LT.(POX+CUX/2) .AND.
     &       X(2,J).GT.(POY-CUY/2) .AND. X(2,J).LT.(POY+CUY/2) .AND.
     &       X(3,J).GT.(POZ-CUZ/2) .AND. X(3,J).LT.(POZ+CUZ/2))THEN
            INDIC=1
          ELSE
            INDIC=0
          ENDIF
        ENDIF
   20 CONTINUE
*
*     Count number of particles within cube
      DO 30 J=1,NTOT
        IF(X(1,J).GT.(POX-CUX/2) .AND. X(1,J).LT.(POS+CUX/2) .AND.
     &     X(2,J).GT.(POY-CUY/2) .AND. X(2,J).LT.(POY+CUY/2) .AND.
     &     X(3,J).GT.(POZ-CUZ/2) .AND. X(3,J).LT.(POZ+CUZ/2))THEN
           COUNT=COUNT+1
        ENDIF
   30 CONTINUE
*
*     Caluculate particle density
      DENS=COUNT/(CUX*CUY*CUZ)
*
*     Output data
      WRITE(*,40) TIME,POX,POY,POZ,COUNT,DENS,INDIC
   40 FORMAT('  DENSaroundFBH',F20.10,E20.7,E20.7,E20.7,I8,E15.4,I5)
*
*     Set particle counter to zero
      COUNT=0
*
      RETURN
*
      END


