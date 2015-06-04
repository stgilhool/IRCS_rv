function f_xy, x, y, n

fxy=(16*x*(1-x)*y*(1-y)*sin(n*!pi*x)*sin(n*!pi*y))^2

return, fxy

end


function f_xy_chi, p, np, funa=funa

n=replicate(18L, np)

;x=funa.x
;y=funa.y



funca=f_xy(p[0,*],p[1,*],n)

chi=1-funca
chi=reform(chi, np)
help, chi
return, chi

end



pro ga_test

;MAKE SURFACE
nx=101L
ny=101L
x=dindgen(nx)/(nx-1)
y=dindgen(ny)/(ny-1)
n=18L
xx=rebin(x,nx,ny)
yy=transpose(rebin(y,nx,ny))






funa={x:xx, y:yy}
ndim=2
lim=[[0,1],[0,1]]
npop=300L

r_val=dblarr(ndim, npop)
save_gfit=dblarr(npop)
;r_fit=dblarr(npop)

;FIND MINIMUM of function
sol=solber('f_xy_chi',ndim, lim=lim, funa=funa, /plot_flag, /print_flag, ngen_max=1000L, save_gen=r_val, save_gfit=save_gfit)

print, sol

;Sort the saved results by fitness
sort_ord=sort(save_gfit)

r_val=r_val[*, sort_ord]
save_gfit=save_gfit[sort_ord]

r_val_uniq=r_val[*,0]
save_gfit_uniq=save_gfit[0]

for i=1,n_elements(save_gfit)-1 do begin
    if r_val[0,i] ne r_val[0,i-1] or r_val[1,i] ne r_val[1,i-1] then begin
        r_val_uniq=[[r_val_uniq], [r_val[*,i]]]
        save_gfit_uniq=[save_gfit_uniq, save_gfit[i]]
    endif
endfor


set_plot, 'x'
load_ct=12
window, 1
fxy=f_xy(xx,yy,n)
contour, 1-fxy, x, y
oplot, r_val_uniq[0,0:100], r_val_uniq[1,0:100], ps=6, color=99999

print, 'x', 'y', 'fitness'
print, r_val_uniq[*, 0:100], reform(save_gfit_uniq[0:100], 1, 101)



stop


end
