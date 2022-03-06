PRO example

zeit = DBLARR(3250)
xpos = DBLARR(3250)
ypos = DBLARR(3250)
zpos = DBLARR(3250)
rpos = DBLARR(3250)
 
;OPENU, lun, '/u/And/omarov/Nbody_4/ARI1.bin, /GET_LUN, /F77_UNFORMATTED
OPENR, lun, '/u/And/omarov/Nbody_4/ARI1.bin, /GET_LUN, /F77_UNFORMATTED

FOR I=0,100 DO BEGIN
;WHILE NOT EOF(inunit) DO BEGIN
;READU, lun,vrem,xi

;zeit(i)= time
;xpos(i)=xi(o)
;ypos(i)=xi(1)
;zpos(i)=xi(2)
READU, inunit,TIME,xpos,ypos,zpos
;ENDWHILE
ENDFOR
	    
;FREE_LUN, lun
;FREE_LUN, inunit

PRINT, TIME,xpos,ypos,zpos
;SET_PLOT, 'PS'
PLOT_IO, TIME, ypos 

END
