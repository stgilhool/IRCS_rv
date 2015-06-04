pro find_window


rootpath='/home/stgilhool/RV_projects/IRCS_rv/data/'
modelpath=rootpath+'supplemental/'

;Const
oversamp=7L
npix_lsf=(oversamp*10L)+1L
bigc=299792458.D
npix=1024L

;Read in calibration results
calib_file='GJ273_18Jan2011_AB1_5_7_7_13.fits'
calib_ext=13
;model_file=calibpath+calib_file
model_file=rootpath+'epoch/18Jan2011/calib_results/'+calib_file
model_par=mrdfits(model_file, calib_ext)

;;;Parameters from earlier calibration
;wl
wl_coeff=model_par.wl_result
wl_scale=model_par.wl_scale


;READ IN PARAMETERS FROM 2-basis 2-run FITS
    ;rfile='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_127_run2sign.fits'
rfile='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_128_255_Mar16bigpen_r1.fits'
f=mrdfits(rfile, 7)
re=f.result

;delta_rv[index]=re[f.delta_rv_index]
delta_wl_coeff=re[f.delta_wl_index]
h2o_depth=re[f.h2o_depth_index]
co2ch4_depth=re[f.co2ch4_depth_index]
delta_wl_coeff=re[f.delta_wl_index]
gh0_coeff=re[f.gh0_coeff_index]
gh1_coeff=re[f.gh1_coeff_index]
other=re[f.other_index]

tau_scale=other[0]	
;params[index]=re



;;;Contstruct wl array and oversampled wl array from coeffs
x=dindgen(npix)
xx=(dindgen(npix*oversamp)-(oversamp/2))/oversamp

wl_soln=poly(x, wl_coeff)
wl_soln_over=poly(xx, wl_coeff)

;READ IN TELLURIC
h2ostr=mrdfits(modelpath+'sky_h2o.fits',1)
co2ch4str=mrdfits(modelpath+'sky_co2andch4.fits', 1)

wl_telluric_long=h2ostr.wave
h2o_long=h2ostr.trans
co2ch4_long=co2ch4str.trans

m24index=where(wl_telluric_long gt 2.275 and wl_telluric_long lt 3.365)
wl_telluric=wl_telluric_long[m24index]
h2o=h2o_long[m24index]
co2ch4=co2ch4_long[m24index]



;Telluric: Adjust optical depths and construct telluric spectrum
h2o_scaled=h2o^h2o_depth
co2ch4_scaled=co2ch4^co2ch4_depth
telluric_long=h2o_scaled*co2ch4_scaled


;put telluric spectrum onto grid
    
telluric=interpol(telluric_long, wl_telluric, wl_soln)


;Find a good window
npix_window=51L
n_best=100L
n_windows=npix-npix_window

best_pixels=dblarr(npix_window, n_best)
best_telluric=dblarr(npix_window, n_best)
best_totals=dblarr(n_best)
best_wl=dblarr(npix_window, n_best)

 for win=0, n_windows-1 do begin
     pixels=lindgen(npix_window)+win
    
     win_vec=telluric[pixels]
    
     value_total=total(win_vec, /double)
    
     lowest_total=min(best_totals, lowest_index)
    
     if value_total gt lowest_total then begin
         best_totals[lowest_index]=value_total
         best_pixels[*, lowest_index]=pixels
         best_telluric[*, lowest_index]=win_vec
         best_wl[*, lowest_index]=wl_soln[pixels]
     endif

 endfor
value_total_vec=dblarr(n_windows)
win_num=lindgen(n_windows)

for win=0, n_windows-1 do begin
    
    pixels=lindgen(npix_window)+win
    win_vec=telluric[pixels]
    value_total_vec[win]=total(win_vec, /double)
    
endfor

plot, win_num, value_total_vec, /xs

stop
    
    

!p.multi=0
window, 1, xsize=1200, ysize=650


for bestwin=0, n_best-1 do begin
    
    plot, best_pixels[*,bestwin], best_telluric[*,bestwin], /xs
    stop
endfor

end
