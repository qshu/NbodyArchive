      SUBROUTINE GW_DECISION(SEP,M1,M2,VSTAR,ECC,DT,RES)              
*
*     Decide if you want to integrate or not  
      INTEGER RES
      REAL*8 SEP,M1,M2,VSTAR,ECC,DT,DE(5)        

      CALL TIDES3(SEP,M1,M2,ECC,VSTAR,DT,DE)
*
      RES=0

*      PRINT *,'GW DECISION b', DT, DE(4)/ECC, DE(3)/SEP, M1, M2, DE   
      IF ((DE(4)/ECC).GT.5.0D-3.OR.(DE(3)/SEP).GT.5.0D-3) THEN
*        do not bypass the GW emission effect
        RES=1
*        PRINT *, 'ALLOW INTEGRATION',DE(4)/ECC,DE(3)/SEP, DT
      END IF   


      RETURN 
   
      END

