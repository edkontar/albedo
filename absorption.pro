function Compton,x

z=1.2 ; number of e- per H-atom

re=2.82d-13
s_T=8.*!PI*re^2/(3.*1d-24)

a1=(1.d0+x)/x^2
a2=2.*(1.d0+x)/(1.+2.*x)-alog(1.d0+2.*x)/x
a3=alog(1.d0+2.*x)/(2.*x)-(1.d0+3.*x)/(1.d0+2.*x)^2

result=2.*!PI*(re*1e+12)^2*(a1*a2+a3)/s_T
return,result

end

function my_Lambda,x
;const=3e-5
;c=3e-5
;c=720./(1.*1.*511.^3)
;c=1e-6 ;7.2e-2/((2.82)^2*(5.11)^3*8.*!PI/3.)

z=1.2 ; number of e- per H-atom
re=2.82e-13
s_T=8.*!PI*re^2/(3.*1e-24)
;Thomson cross-section in cm^-24 units

s_a=fltarr(3)
IF (511./x GT 8.33)                      THEN s_a=[701.2*x^3,0.,0.]
IF (511./x GT 7.11) AND (511./x LT 8.33) THEN s_a=[629.0*x^3,0.,0.]
IF (511./x GT 4.04) AND (511./x LT 7.11) THEN s_a=[433.9*x^3,-2.4*511.*x^2,0.75*x*511.^2]
IF (511./x GT 1.00) AND (511./x LT 4.04) THEN s_a=[352.9*x^3,18.7*511.*x^2,0.]
s_a=s_a/(511.^3)

;s_a=s_a/(511.^3)

return,z*S_T/(z*S_T+total(s_a))
;return,total(s_a)*511.^3/x^3

end

function LnLambda,x
;const=3e-5
;c=3e-5
;c=720./(1.*1.*511.^3)
;c=701.2/(0.665*511.^3)
z=1.2 ; number of e- per H-atom
re=2.82e-13
s_T=8.*!PI*re^2/(3.*1e-24)
;Thomson cross-section in cm^-24 units

s_a=fltarr(3)

IF (511./x GT 8.33)                      THEN s_a=[701.2*x^3,0.,0.]
IF (511./x GT 7.11) AND (511./x LT 8.33) THEN s_a=[629.0*x^3,0.,0.]
IF (511./x GT 4.04) AND (511./x LT 7.11) THEN s_a=[433.9*x^3,-2.4*511.*x^2,0.75*x*511.^2]
IF (511./x GT 1.00) AND (511./x LT 4.04) THEN s_a=[352.9*x^3,18.7*511.*x^2,0.]

s_a=s_a/(511.^3)

return,alog(z*S_T/(z*S_T+total(s_a)))

end

function ln_int,x
;c=3e-5
;c=720./(1.*1.*511.^3)
;c=701.2/(0.665*511.^3)

s_a=fltarr(3)
IF (511./x GT 8.33)                      THEN s_a=[701.2*x^3,0.,0.]
IF (511./x GT 7.11) AND (511./x LT 8.33) THEN s_a=[629.0*x^3,0.,0.]
IF (511./x GT 4.04) AND (511./x LT 7.11) THEN s_a=[433.9*x^3,-2.4*511.*x^2,0.75*x*511.^2]
IF (511./x GT 1.00) AND (511./x LT 4.04) THEN s_a=[352.9*x^3,18.7*511.*x^2,0.]
s_a=s_a/(0.665*511.^3)

part1=x*(3.-alog(1.+c*x^3))+alog((1.-c^(1./3.)*x+c^(2./3.)*x^2)/(1.+c^(1./3)*x)^2)/(2.*c^(1./3.))
part2=-sqrt(3)*atan((2.*c^(1./3.)*x-1.)/sqrt(3))

;return, part1+part2
return,0

end

pro absorption,out,mu,y0,dy
; absorption

out=1.

IF (dy LT 2.) THEN BEGIN
a_mu=0.802-1.019*mu+2.528*mu^2-3.198*mu^3+1.457*mu^4+8.1e-7*mu^(-4)
h_dy=exp(0.381*dy)-2.
int1=QROMO('lnlambda', y0,y0+dy)
;out =my_Lambda(y0)*(a_mu+(1.-a_mu))*(1.+h_dy*sqrt(1.-my_Lambda(y0)))*exp(ln_int(y0+dy)-ln_int(y0))
out =my_Lambda(y0)*(a_mu+(1.-a_mu))*(1.+h_dy*sqrt(1.-my_Lambda(y0)))*exp(int1)
END

IF (dy GT 2) THEN BEGIN
int2= QROMO('lnlambda',y0,y0+dy-0.15)
;print,'INt=',int2
;out =exp(ln_int(y0+dy-0.15)-ln_int(y0))
out =exp(int2)
END

IF (OUT LT 1e-6) THEN out=0.
IF (OUT GT 1) THEN out=1.

;print, 'absorption =',out

end

pro non_rel_limit,result,mu,y0,dy
; green function in non-relativistic limit

L0=my_Lambda(y0)

IF (L0 LT 0.9 ) THEN BEGIN

c1=3.-2.*mu^2+3.*mu^4
c2=alog(1.+1./mu)
c3=(3.*mu^2-1.)*(0.5-mu)
F1=3.*L0*mu*(c1*c2+c3)/16.


k1=1./sqrt(1.-L0)-1
k2=1.-alog(1.+sqrt(3.*(1.-L0))/sqrt(3.*(1-L0)))
k3=(1.+mu*sqrt(3))/(1.+mu*sqrt(3.*(1-L0)))
br=(1.+k1*k2)*k3-1.
F2=L0*mu*alog(1.+1./mu)*br/2.
result=F1+F2

END


end
