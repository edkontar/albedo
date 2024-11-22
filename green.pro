pro green,out,mu,y0,dy

x0=1./y0
green=0.

IF (x0 LT 1e-2) THEN begin
fmu=1.844*mu^0.844

;Green_WLZ88
;angle averaged cross-section
;White,Lightman,Zdziarski,1988,ApJ,331,939
;green=fmu*green_wlz88

;wlz88,value,y0,dy

wlz88_vol2,value,y0,dy

green=value*fmu

end

IF ((x0 GT 1e-2) AND (x0 LT 31.6)) then begin
    IF (dy LT 2.) THEN BEGIN
    ;*********************************************
    c=fltarr(4,5,8)
    c[0,0,*]=[ 0.8827,-0.5885, 0.7988,-0.4117,0,0,0,0]
    c[0,1,*]=[ 0.0517, 0.1076,-0.1691, 0.0678,0,0,0,0]
    c[0,2,*]=[ 0.0014, 0.0043,      0,      0,0,0,0,0]
    c[0,3,*]=[-0.0003, -0.0027,0.0021,      0,0,0,0,0]

    c[1,0,*]=[-1.3259, 1.6214,-3.7007, 2.1722, 0, 1, 20.3,1.063]
    c[1,1,*]=[ 0.0790,-0.4029, 0.6619, 0.2210, 0,-1, 17.7,1.041]
    c[1,2,*]=[ 0.0751,-0.1848, 0.4068,-0.4126, 0, 1, 09.8,1.185]
    c[1,3,*]=[-0.0020,-0.0394, 0.1004,-0.0597, 0, 0,    0,    0]

    c[2,0,*]=[ 3.2953,-3.6996, 7.9837,-4.6855, 0, 0,    0,    0]
    c[2,1,*]=[-0.5278, 0.9857,-1.7454,-0.3464, 0, 1, 27.6,0.959]
    c[2,2,*]=[-0.1919, 0.5798,-1.2879, 1.4885, 0,-1, 18.3,0.986]
    c[2,3,*]=[ 0.0200, 0.0832,-0.0333,-0.2370, 0, 1, 16.6,1.086]

    c[3,0,*]=[-2.2779, 2.4031,-4.0733, 1.9499, 0,     -1, 19.6,0.968]
    c[3,1,*]=[ 0.4790,-1.0166, 3.1727,-4.0108, 3.0545,-1, 30.4,0.957]
    c[3,2,*]=[ 0.1122,-0.3580, 0.4985,-0.3750,-0.5349, 1, 31.2,0.972]
    c[3,3,*]=[-0.0160,-0.0471,-0.0536, 0.2006, 0.2929,-1, 30.6,1.009]
    c[3,4,*]=[ 0.0005, 0.0021,-0.0447, 0.1749,-0.2303, 1, 16.9,1.130]
    ;*********************************************

    b=fltarr(4,5)
    for i=0,3 do begin
    for j=0,4 do begin
    b[i,j]=total(c[i,j,0:4]*[1,mu^1,mu^2,mu^3,mu^4])+c[i,j,5]*exp(c[i,j,6]*(mu-c[i,j,7]))
    end
    end

    a0=fltarr(4)
    for i=0,3 do a0[i]=alog(y0)*(total(b[i,*]*[1.,alog(y0),alog(y0)^2,alog(y0)^3]))

    log_log=alog((2.+y0)/(dy+y0))/alog((2.+y0)/y0)

    sum=total(a0*[1,log_log,log_log^2,log_log^3])

    green1,f,mu,dy

    green= f*((dy+1)/(dy+y0))*exp(sum)

   IF ((x0 GT 1e-2) AND (x0 LT 1e-1)) THEN BEGIN
     b0=0.006+0.089*mu-0.102*mu^2+0.056*mu^3
     a4=1.+b0*(0.326*((y0/10.)^1.692-1.)+alog(y0/10.))
     green,green10,mu,9.999,dy

     green= green10*a4
     ;needs a value green at y0=10
   END
   END


   IF (dy GT 2.) THEN BEGIN

   cc=fltarr(8,4)

   cc[1,*]=[ 0.618,-0.830,  0,     0]
   cc[2,*]=[ 0.128,-0.132,  0,     0]
   cc[3,*]=[ 0.632,-0.875,  0,      0]
   cc[4,*]=[-0.672,-0.286,  0.717, -0.481]
   cc[5,*]=[0.0126,-0.0160, 0.0077, 0]
   cc[6,*]=[0.0111, 0.0030,-0.0014, 0]
   cc[7,*]=[-2.437,-0.3280,-0.2600, 0.279]

   bb=fltarr(8)
   for j=1,7 do bb[j]=total(cc[j,*]*[1,mu,mu^2,mu^3])

   a5=0.763*y0^(-0.507)+0.007*y0^(-1.553)+0.033*alog(y0)+0.230
   a6=75.180/y0
   a7=0.100*y0^1.623+0.115*alog(y0)-0.100
   a8=(y0 LT 1.)*(bb[1]+bb[2]*alog(y0))+(y0 GT 1.)*(bb[3]+bb[4]*alog(y0)+bb[5]*(bb[6]+1./y0)^bb[7])

   green1,f,mu,dy

   green=f*((dy+1.)/(dy+y0))*a5*exp(a7*(1.+a8/(dy+y0))*alog(1.+a6/(dy+y0)))

   END



END

dycc=1.-sqrt(1.-(mu-0.05)^2)
IF (dy LT dycc) then green=0.

;end function
;***********************************************************************************
;return, green

;IF (511./y0 GT 12. ) THEN BEGIN
absorption,coef,mu,y0,dy
;print, 'Green function  =',out
out=green*coef

;ENDIF  ELSE BEGIN
;non_rel_limit,result,mu,y0,dy
; solution of transfer equations
;out =result

;ENDELSE

;out=green

end
