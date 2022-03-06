      SUBROUTINE SMBH_ADJUST
*
*
*       Positions and velocity adjustment due to the SMBH.
*       --------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8 RI(NMAX),R_SHELL,R_SHELL_OLD,RHO_SHELL,
     &       M_SHELL,M_TOT,F_TOT,POT_I,F_BNDR,POT_SH
      SAVE R_SHELL_OLD
      DATA R_SHELL_OLD/-3/
      SAVE R_MID_OLD
      DATA R_MID_OLD/-3.5/
*
      if(rank.eq.0) then
*        DO 85 I = 1,N
         DO 85 I = IFIRST,NTOT
             WRITE (10,100) BODY(I),(X(K,I),K=1,3),(XDOT(K,I),K=1,3)
   85     CONTINUE
  100 FORMAT(1X,1P,7(1X,D18.10))
      end if
*      
      M_TOT = 0.0
      POT_I = 0.0
      F_TOT = 0.0
      POT_CH = 0.0
      DO 96 R_SHELL = -2, 1, 0.05
      M_SHELL=0.0
*      
      DO 95 I = IFIRST,NTOT
      RI(I)=(X(1,I)**2.0+X(2,I)**2.0+X(3,I)**2.0)**0.5
      IF(log10(RI(I)).LE.R_SHELL.AND.log10(RI(I)).GT.R_SHELL_OLD) THEN
      M_SHELL=M_SHELL+BODY(I)   
      POT_SH = POT_SH + 
     &BODY(I)*(PHIDBL(I)+SMBH/DSQRT(X(1,I)**2+X(2,I)**2+X(3,I)**2))/2.0
      END IF
   95 CONTINUE
*   
      R_MID = (10**R_SHELL + 10**R_SHELL_OLD)/2.0
      M_TOT = M_TOT + M_SHELL
      F_TOT = (F_BNDR + M_TOT/(10**R_SHELL)**2.0)/2.0
      POT_I = POT_I - F_BNDR*(R_MID - R_MID_OLD) 
      RHO_SHELL=M_SHELL/(2.0/3.0*TWOPI*
     &             (10**(R_SHELL*3.0)-10**(R_SHELL_OLD*3.0)))
      R_SHELL_OLD = R_SHELL
      R_MID_OLD = R_MID
      WRITE(36,8) R_MID,RHO_SHELL,M_TOT,F_TOT,POT_I 
      F_BNDR = M_TOT/(10**R_SHELL)**2.0
   8  FORMAT(5(D15.5))
      PRINT*,'POT_I & POT_SH', POT_I,POT_SH, POT-EBLCKHL
   96 CONTINUE    
      STOP
*        
      RETURN
      END
