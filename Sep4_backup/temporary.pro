pro check_lsf

rfile='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Feb10run_lsfprobs/GJ273_18Jan2011_AB_0_127.fits'

n_ext=14L

fiducial=mrdfits(rfile, 1)
oversamp=fiducial.oversamp
npix_select=fiducial.npix_select
npix_lsf=oversamp*10L + 1
npix_model=npix_select*oversamp
n_gh0=n_elements(fiducial.gh0_coeff_index)
n_gh1=n_elements(fiducial.gh1_coeff_index)
n_wl=n_elements(fiducial.delta_wl_index)
n_other=n_elements(fiducial.other_index)

delta_rv=dblarr(n_ext)
final_rv=dblarr(n_ext)
mjd=dblarr(n_ext)
gh0=dblarr(n_ext, n_gh0)
gh1=dblarr(n_ext, n_gh1)
h2o=dblarr(n_ext)
co2ch4=dblarr(n_ext)
d_wl=dblarr(n_ext, n_wl)
other=dblarr(n_ext, n_other)
lsf=dblarr(n_ext, npix_lsf, npix_model)

for visit=0,n_ext-1 do begin

    str=mrdfits(rfile, visit+1)
    delta_rv[visit]=str.delta_rv
    mjd[visit]=str.mjd
    final_rv[visit]=str.final_rv
    gh0[visit,*]=str.result[str.gh0_coeff_index]
    gh1[visit,*]=str.result[str.gh1_coeff_index]
    
    lsf[visit,*,*]=ircsrv_lsf(gh0[visit,*], gh1_vec=gh1[visit,*], npix_select=npix_select)

endfor

window, 1, xsize=1800L, ysize=1000L

!p.multi=[0,3,1]
x1=npix_select*oversamp/3
x2=npix_select*oversamp/2
x3=2*npix_select*oversamp/3

for visit=0,13 do begin

    plot, lsf[visit,*,x1], title="Exposure no: "+ strtrim(visit,2)
    plot, lsf[visit,*,x2]
    plot, lsf[visit,*,x3]
stop
endfor

end
