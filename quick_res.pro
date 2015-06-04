pro quick_res

file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_127_t7newpen_r2free_spec.fits'

model=dblarr(128L, 14)
obs=dblarr(128L, 14)
err=dblarr(128L, 14)
wl=dblarr(128L, 14)
tell=dblarr(128L, 14)
nh3=dblarr(128L, 14)
star=dblarr(128L, 14)
restest=dblarr(128L, 14)

for vis=0,13 do begin
    
    a=mrdfits(file, vis+1)
    model[*,vis]=a.model
    obs[*,vis]=a.obs
    err[*,vis]=a.error
    wl[*,vis]=a.wl_grid
    tell[*,vis]=a.telluric
    nh3[*,vis]=a.nh3
    star[*,vis]=a.stellar
    restest[*,vis]=a.residuals

endfor

res=model-obs


!p.multi=0

plot, wl, res, ps=3, title="All residuals"



stop

end
