function convert_func, p
common amoebainfo, x_iraf, wl_soln_iraf

trial_wl=poly(x_iraf, p)

rmse=total(sqrt((trial_wl-wl_soln_iraf)^2), /double)

print, p
print, rmse

plot, wl_soln_iraf
oplot, trial_wl, color=200
wait, 0.2

return, rmse

end 


pro convert_IRAF;, func, ord, xmin, xmax, c
;Quick script to convert IRAF wl soln to ordinary polynomial
;coefficients for my IDL IRCS_cal code
;FUNC: function id (2 for legendre)
;ORD: order
;XMIN: db entry after order
;XMAX: db entry after xmin
;C: coefficients (array of remaining entries)

func=2
ord=5
xmin=9.63d0
xmax=987.53d0
c=dblarr(ord)
c[0]=2.322905043791022d0
c[1]=-0.02873394489300029d0
c[2]=-7.423267763775305d-4
c[3]=7.686501323716287d-5
c[4]=-6.055156262998915d-6

common amoebainfo, x_iraf, wl_soln_iraf

;Define stuff
npix=1024L
wl_iraf=dblarr(npix, ord)
n_iraf=dblarr(npix)

z=dblarr(npix, ord)

;Make independent variable vector
x_iraf=dindgen(npix)

n_iraf=(2d0*x_iraf-(xmax+xmin))/(xmax-xmin)

;Make z vectors
z[*,0]=replicate(1d0, npix)
z[*,1]=n_iraf

;Compute wavelength
for i=0, ord-1 do begin
    if i ge 2 then z[*,i]=((2*(i+1)-3) * n_iraf * z[*,i-1]-(i-1)*z[i-2])/(i)
    wl_iraf[*,i]=c[i]*z[*,i]
endfor

wl_soln_iraf=total(wl_iraf, 2, /double)

wl_soln_iraf=reverse(wl_soln_iraf)

plot, wl_soln_iraf

stop



guess=[2.2912196d0,6.2602232d-05,-2.4059731d-09,-1.0215283d-12,-3.4594357d-16]
scale=[5d-05,5d-08,1d-11,1d-13,5d-18]
ftol=1d-10

r=amoeba3(ftol, scale=scale, p0=guess, function_name='convert_func', function_value=fval)
stop
print, r[0:4]
print, fval[0]

end
