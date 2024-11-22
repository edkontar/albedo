;+
; PROJECT:
;   RHESSI/albedo
;
; PURPOSE:
;
;   Calculates Green matrix for X-ray analysis following the paper:
;   https://ui.adsabs.harvard.edu/abs/2006A%26A...446.1157K/abstract

; CATEGORY:
;   HESSI, Spectra
;
; CALLS: No externals
;
; INPUTS: none
;
;
; OUTPUTS: IDL save files for different angles:
; eg green_compton_mu005.dat
; where mu=0.05 is the cos of the heliocentric angle, so cos(theta)=0.05 
;
; OPTIONAL OUTPUTS:
;   plots reference albedo correction to compare with the data
;   https://ui.adsabs.harvard.edu/abs/2007A%26A...466..705K/abstract
;
; KEYWORDS:
;   none
;
; COMMON BLOCKS:
;   none
;
; SIDE EFFECTS:
;
;
; RESTRICTIONS:
;
;
; MODIFICATION HISTORY:
;   Version 1, eduard(at)astro.gla.ac.uk, 23 November, 2005
;   returns a structure with resolution and errors
;- Eduard@glasgow updated names to upload to github, November 2024


Function Fin,x
return,exp(-x*511./300.)*x^(-.9)
end

pro make_all_green_1kev

RESOLVE_ROUTINE,'absorption', /COMPILE_FULL_FILE
RESOLVE_ROUTINE,'green', /COMPILE_FULL_FILE

mc2=511.
N=597
Njana=N

;M=303

mu_arr=fltarr(18)/20.+0.5

;mu0=cos(45.*!PI/180.)
;print, '!4m!3 =',mu0


EBinNumber=N
e2N=fltarr(2,EbinNumber)
base=(1000.-2.-float(EBinNumber)*.3)^(1./EBinNumber)
Ebin=2.+findgen(EBinNumber+1)*.3+base^findgen(EBinNumber+1)
;e2N(0,*)=EBin(0:N_elements(EBin)-2)
;e2N(1,*)=EBin(1:N_elements(EBin)-1)

e2N(0,*)=findgen(N)+3.
e2N(1,*)=findgen(N)+4.



for s=0,18 do begin

;s=18
mu0=0.05+s*0.05

;mu0=0.77

print,'current mu =',mu0

save_file='green_compton_mu'+string(format='(I3.3)',5*FLOOR(mu0*20.),/print)+'.dat'



e2jana=e2N

de=e2jana(1,*)-e2jana(0,*)
de=transpose(de)

e=(e2N(0,*)+e2N(1,*))/2.

;E =transpose(900./float(findgen(N)+1.) )
;E=findgen(N)*4.+1.
x =E/mc2
x2N=e2N/mc2
x0=x

Fref=fltarr(N)
A=fltarr(N,N)

;stop

y =1./x
y0=1./x0
yy=511./e2N
;dy=transpose(yy(0,*)-yy(1,*))
GreenF=fltarr(100)

y0_max=511./4.
deltay=y0_max*1.2^(-findgen(100))

for k=0,99 do begin
green,aa,mu0,y0_max,deltaY(k)
greenF(k)=aa
end

plot_oo,deltaY,greenF,yrange=[1e-4,1],psym=1,xrange=[1e-2,1e3]

;stop

For i=1,N-1 do begin
; output energy domain

   for j=i,N-1 do begin
   ; input energy domain

   DeltaY_x0x=y(i)-y0(j)
   DeltaY    =y0(j-1)-y0(j)

   IF (DeltaY GT 2.) THEN BEGIN
   Ndy=90
   N2=40
   dy2N=fltarr(2,Ndy)

   dy2N(0,0:N2-1)=(findgen(N2)+0.)*2./(N2)
   dy2N(1,0:N2-1)=(findgen(N2)+1.)*2./(N2)

   dy2N(0,N2:Ndy-1)=1.+(deltaY-1.)^((findgen(Ndy-N2)+0.)/(Ndy-n2))
   dy2N(1,N2:Ndy-1)=1.+(deltaY-1.)^((findgen(Ndy-N2)+1.)/(Ndy-n2))


   ENDIF ELSE BEGIN
   Ndy=1+floor(DeltaY*60.)
   dy2N=fltarr(2,Ndy)

   dy2N(0,0:Ndy-1)=(findgen(Ndy)+0.)*deltaY/(Ndy)
   dy2N(1,0:Ndy-1)=(findgen(Ndy)+1.)*deltaY/(Ndy)
   ENDELSE

   dy_m=(dy2N(1,*)+dy2N(0,*))/2.
   dy  = dy2N(1,*)-dy2N(0,*)

   sum=0.
   for k=0, Ndy-1 do begin
   green,aa,mu0,y0(j)+dy_m(k),DeltaY_x0x+dy_m(k)
   A(i,j)=A(i,j)+aa*dy(k)
   END


    ;print,'Energy =',e(j),' keV ...................................... OK'
   END

;IF (i EQ 10 ) THEN STOP
print,'Energy =',e(i),' keV ...................................... OK'


IF (e(i) LT 50.) THEN BEGIN
non_rel_limit,result,mu0,511./e(i),dy
IF (Total(A(i,*)) LT result) THEN  A(i,*)=A(i,*)*result/Total(A(i,*))
END


END

Ajana=A
;Ajana=fltarr(Njana,Njana)
;test=fltarr(Njana,Njana)

;for i=0, Njana-1 do begin
;irange=round([e2jana(0,i)-3.,e2jana(1,i)-3.])
;for j=0, Njana-1 do begin
;jrange=round([e2jana(0,j)-3.,e2jana(1,j)-3.])

;Ajana(i,j)=total($
;total(A(irange(0):irange(1),jrange(0):jrange(1)),1)$
;)/(de(i))

;end
;end


b=transpose(Ajana)

x=transpose((e2jana(1,*)+e2jana(0,*)))/(2.*511.)

;b=smooth(b,5)

f0=fin(x)
fa=b##(f0)

e=x*511.

window,0
plot_oo,e,f0*x, yrange=[1e-2,1],ystyle=1
oplot,e,Fa*x,line=1
oplot,e,(F0+Fa)*x,line=2


p = {cos_theta:mu0,edges:e2jana,albedo:b}
SAVE,FILENAME = save_file,p

ff=exp(-e/2.)+e^(-2)
one=fltarr(Njana,Njana)
for k=0, Njana-1 do one(k,k)=1.
b1=invert(one+b)

fb=b1##ff
window,1
!P.Multi=[0,1,2]
plot_oo,e,fb
oplot,e,ff,line=1

plot,e,-smooth(deriv(alog(e),alog(fb)),5),/xlog
!P.Multi=0


F17=(3./e)^(2)
a17=b##f17

window,4
plot_oo,e,F17,xrange=[3,500],yrange=[1e-5,1.],xstyle=1,ytitle='Photon Flux',xtitle='Energy, keV'
oplot,e,a17,line=1
oplot,e,a17+F17,line=2

;figure for

end


f0=1e8*(exp(-e/2.)+1./e^(2.5))
a0=b##f0

window,5
plot_oo,e,f0,xrange=[3,500],xstyle=1,ytitle='Photon Flux',xtitle='Energy, keV'
oplot,e,a0,line=2
oplot,e,f0+a0,line=1

Set_plot, 'PS'
Device,filename='albedo_mu077_pow2_reference.ps', xsize=12,ysize=12,xoffset=2,yoffset=2,ENCAPSULATED=0

!P.Multi=[0,1,2]

plot_oo,e,F17,xrange=[3,500],yrange=[1e-5,1.],xstyle=1,ytitle='Photon Flux',xtitle='Energy, keV'
oplot,e,a17,line=1
oplot,e,a17+F17,line=2

plot_oo, e, a17/f17, xrange=[3,500],yrange=[1e-2,1],xstyle=1,ytitle='Reflectivity',xtitle='Energy, keV'

!P.Multi=0

DEVICE, /CLOSE
SET_PLOT, 'win'


print,'All files calculated ...................... OK'

 stop

end