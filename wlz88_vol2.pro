pro wlz88_vol2,greenf,y0,dy
; angle averaged green function
; from

C=0.10
D=0.56

IF (dy LT 2.) THEN greenf=C
IF (dy GT 2.) THEN greenf=D*dy^(-1.5)

;print,'the value is:',greenf
end