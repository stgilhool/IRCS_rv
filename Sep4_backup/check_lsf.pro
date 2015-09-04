pro check_lsf


;rfile='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_127_smooth30_r2.fits'
rfile='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_290_417_Mar23_mpfit2_wl_gh_r3.fits'
;rfile='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_127_smooth.fits'


;r2file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_127_four2.fits'
;r2file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_127_run2sign.fits'
r2file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_290_417_Mar23_mpfit2_wl_gh_r1.fits'
;sfile='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_127_four2_spec.fits'

sfile='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_127_smooth30_r2_spec.fits'

write_results=0
;fits_info, rfile, /silent, n_ext=n_ext
n_ext=14L

fiducial=mrdfits(rfile, 1)
oversamp=fiducial.oversamp
npix_select=fiducial.npix_select
first_pix=fiducial.first_pix
npix_lsf=oversamp*10L + 1
npix_model=npix_select*oversamp
n_gh0=n_elements(fiducial.gh0_coeff_index)
n_gh1=n_elements(fiducial.gh1_coeff_index)
n_wl=n_elements(fiducial.delta_wl_index)
n_other=n_elements(fiducial.other_index)

delta_rv=dblarr(n_ext)
rv_guess=dblarr(n_ext)
final_rv=dblarr(n_ext)
bcv=dblarr(n_ext)
d_bcv=dblarr(n_ext)
mjd=dblarr(n_ext)
gh0=dblarr(n_ext, n_gh0)
gh1=dblarr(n_ext, n_gh1)
h2o=dblarr(n_ext)
co2ch4=dblarr(n_ext)
d_wl=dblarr(n_ext, n_wl)
other=dblarr(n_ext, n_other)
lsf=dblarr(n_ext, npix_lsf, npix_model)
lsf_penalty=dblarr(n_ext)
model=dblarr(n_ext, npix_select)
obs=dblarr(n_ext, npix_select)
res=dblarr(n_ext, npix_select)
status=dblarr(n_ext)
wl_grid=dblarr(n_ext, npix_select)
err=dblarr(n_ext, npix_select)
stellar=dblarr(n_ext, npix_select)
telluric=dblarr(n_ext, npix_select)
nh3=dblarr(n_ext, npix_select)
chi2=dblarr(n_ext)
chi22=dblarr(n_ext)

delta_rv2=dblarr(n_ext)
for visit=0,n_ext-1 do begin
    ss=mrdfits(r2file, visit+1)
    delta_rv2[visit]=ss.delta_rv
    chi22[visit]=ss.chi2
endfor


for visit=0,n_ext-1 do begin

    str=mrdfits(rfile, visit+1)
    delta_rv[visit]=str.delta_rv
    rv_guess[visit]=str.guess[str.delta_rv_index]
    mjd[visit]=str.mjd
    final_rv[visit]=str.final_rv
    gh0[visit,*]=str.result[str.gh0_coeff_index]
    gh1[visit,*]=str.result[str.gh1_coeff_index]
    h2o[visit]=str.result[str.h2o_depth_index]
    co2ch4[visit]=str.result[str.co2ch4_depth_index]
    d_wl[visit, *]=str.result[str.delta_wl_index]
    other[visit, *]=str.result[str.other_index]
    bcv[visit]=str.bcv
    d_bcv[visit]=str.delta_bcv
    chi2[visit]=str.chi2
    status[visit]=str.status

    
    lsf[visit,*,*]=ircsrv_lsf(gh0[visit,*], gh1_vec=gh1[visit,*], first_pix=first_pix, npix_select=npix_select, neg_penalty=lsf_pen)
    lsf_penalty[visit]=lsf_pen
if write_results eq 1 then begin
     ocs=mrdfits(sfile, visit+1)
     obs[visit,*]=ocs.obs
     model[visit,*]=ocs.model
     res[visit,*]=ocs.residuals
     err[visit,*]=ocs.error
     stellar[visit,*]=ocs.stellar
     telluric[visit,*]=ocs.telluric
     nh3[visit,*]=ocs.nh3
     chi2[visit]=ocs.chi2
     wl_grid[visit,*]=ocs.wl_grid
endif
    

 endfor


if write_results eq 1 then begin
    o=replicate({bcv:0d0, delta_bcv:0d0, mjd:0d0, final_rv:0d0, $
                 delta_rv:0d0, $
                 h2o_depth:0d0, $
                 co2ch4_depth:0d0, $
                 delta_wl:replicate(0d0,n_wl), $
                 gh0_coeff:replicate(0d0, n_gh0), $
                 gh1_coeff:replicate(0d0, n_gh1), $
                 other:replicate(0d0, n_other), $
                 wl_grid:replicate(0d0, npix_select), $
                 obs:replicate(0d0, npix_select), $
                 model:replicate(0d0, npix_select), $
                 stellar:replicate(0d0, npix_select), $
                 telluric:replicate(0d0, npix_select), $
                 nh3:replicate(0d0, npix_select), $
                 sigma:replicate(0d0, npix_select), $
                 chi2:0d0 $
                }, n_ext)
    
    for i=0, n_ext-1 do begin
        o[i].bcv=bcv[i]
        o[i].delta_bcv=d_bcv[i]
        o[i].mjd=mjd[i]
        o[i].final_rv=final_rv[i]
        o[i].delta_rv=delta_rv[i]
        o[i].h2o_depth=h2o[i]
        o[i].co2ch4_depth=co2ch4[i]
        o[i].delta_wl=d_wl[i,*]
        o[i].gh0_coeff=gh0[i,*]
        o[i].gh1_coeff=gh1[i,*]
        o[i].other=other[i,*]
        o[i].wl_grid=wl_grid[i,*]
        o[i].obs=obs[i,*]
        o[i].model=model[i,*]
        o[i].stellar=stellar[i,*]
        o[i].telluric=telluric[i,*]
        o[i].nh3=nh3[i,*]
        o[i].sigma=err[i,*]
        o[i].chi2=chi2[i]
    endfor
    

 outfile= '/home/stgilhool/RV_projects/IRCS_rv/smoothingresults.fits'            
stop
mwrfits, o, outfile                     
endif



window, 1, xsize=1800L, ysize=1000L

!p.multi=[0,3,1]
x1=npix_select*oversamp/3
x2=npix_select*oversamp/2
x3=2*npix_select*oversamp/3

for visit=0,n_ext-1 do begin
;stop
fincheck=where(finite(lsf[visit,*,*]) eq 0, nancount)
if nancount eq 0 then begin
titlestr="Exposure no: "+ strtrim(visit,2) + $
  " | chi2: " + strtrim(chi2[visit],2) + $
  " | lsf_penalty: " + strtrim(lsf_penalty[visit],2)
    plot, lsf[visit,*,x1], title=titlestr
    plot, lsf[visit,*,x2]
    plot, lsf[visit,*,x3]
stop
endif
endfor


!p.multi=0
plot, mjd, final_rv, ps=6, symsize=1.5
oplot, mjd, delta_rv2+d_bcv, ps=6, color=200






stop

end
