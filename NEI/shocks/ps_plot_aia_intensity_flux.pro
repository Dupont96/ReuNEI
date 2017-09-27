pro ps_plot_aia_intensity_flux


set_plot,'ps'
!p.font=0
device,xsize=16,ysize=12,xoff=2,$
yoff=2,bit=8,filename='ps_aia_intenflux.eps',/color,Encapsulated=1

restore,'AIA_intensity_0193.sav',/ver         
time1=time
utbase1=utbase
int1=intensity
uttime1=uttime
restore,'AIA_intensity_0211.sav',/ver 
time2=time
utbase2=utbase
int2=intensity
uttime2=uttime
restore,'AIA_intensity_0171.sav',/ver
time3=time
utbase3=utbase
int3=intensity
uttime3=uttime

restore,'AIA_intensity_0335.sav',/ver
time4=time
utbase4=utbase
int4=intensity
uttime4=uttime


!P.background=255
!P.color=0    
linecolors
utplot,time1,int1,utbase1,yrange=[10,55],/xstyle,/ystyle,position=[0.12,0.15,0.9,0.9],$
charsize=1.,thick=1.5,/norm,/nodata,ytitle='Intensity (DN s!E-1!Npix!E-1!N)',$
xtitle='Time' 
oplot,time1,int1,color=3,thick=2
oplot,time2,int2*5,color=11,thick=2
oplot,time3,int3,color=8,thick=2
oplot,time4,int4*250,color=6,thick=2
;oplot,time4,smooth(int4*250,3),color=6,thick=2

oplot,[534,534],[10,55],linestyle=2,color=0,thick=1.5
oplot,[755,755],[10,55],linestyle=3,color=0,thick=1.5

xyouts,0.42,0.55,'193',color=3,charsize=1,/norm  
xyouts,0.28,0.75,'211x5',color=11,charsize=1,/norm
xyouts,0.7,0.44,'171',color=8,charsize=1,/norm
xyouts,0.5,0.7,'335x250',color=6,charsize=1,/norm

xyouts,0.31,0.2,'Shock',color=0,charsize=1,/norm
xyouts,0.42,0.2,'CME',color=0,charsize=1,/norm

device,/close
set_plot,'X'

end

