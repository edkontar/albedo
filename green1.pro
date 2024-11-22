pro green1, f, mu,dy
; green function for x0=1
; e.g. for 511 keV


dyc=1.-sqrt(1.-mu^2)
dyd=1.+sqrt(1.-mu^2)

IF ((dy LT 2.) AND (dy GT dyd)) then begin
 b0=0.911-0.549*(1.-mu)^(1.471)
 b1=0.254-0.041*mu^(-3.798)
green1=(b0+b1*(dy-2.))/(1.+dy)
end

IF (dy LT dyd) then begin
b2=0.161+0.439*mu^0.189-0.791*mu^7.789
b3=-0.871+1.740*mu^0.176-1.088*mu^0.907
b3=(0. LT b3)*b3 +0.*(0. GT b3)
b4=0.934+0.054*mu^(-0.666)
b5=(mu LT 0.524)*(4.647-11.574*mu+11.046*mu^2)+(mu GT 0.524)*1.615
b6=0.012+0.199*mu^.939+0.861*mu^8.999
b7=2.150
b8=-1.281+2.374*mu-4.332*mu^2+2.630*mu^3
b9=2.882+0.035*mu^(-1.065)
b10=30.28+69.29*mu
b11=1.037*(1.-sqrt(1.-0.874*mu^2))^(0.123)
b12=0.123
b_max=(b4-dy)
b_max=(0. LT b_max)*b_max+0.*(0. GT b_max)
b_max=b_max^b5
;print,b_max
frac1=(b2+b3*b_max +b6*(b7-dy)^b8)^b9
frac2=(1.+dy)*(1.+exp(b10*(b11-dy^b12)))
green1=frac1/frac2
end

IF (dy GT 2.) THEN BEGIN
b13=56.50-27.54*mu+671.8*mu^2-1245.*mu^3+708.4*mu^4
b14=-0.897-0.826*mu+2.273*mu^2-0.984*mu^3
b15=1.240-1.297*mu
b16=-0.490+1.181*mu-1.038*mu^2
fas=1.448*(mu^2.033+0.745*mu^1.065)
gas=1.283*fas*(dy)^(-1.5)
green1=GAS*exp((b14+b15*(dy+1.)^b16)*alog(1.+b13/(dy+1.)) )
END

;green1=green1/x^2
;green1=green1*(1.+dy)

f=green1
;return, f

end
