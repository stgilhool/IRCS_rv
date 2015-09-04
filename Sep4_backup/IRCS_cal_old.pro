function amoebafunction, p

common amoeba_info, wl_lab, int_lab, int_obs

c=299792458.D
npix=1024L
;oversamp=7D

;;;Parameters from amoeba
;don't need rad_vel for initial wl calibration
; rad_vel=p[0]
; wl_0=p[1]
; wl_1=p[2]
; wl_2=p[3]
; wl_3=p[4]


wl_0=p[0]
wl_1=p[1]
wl_2=p[2]
wl_3=p[3]
tau_scale=p[5]

;;;vectors for oversampling
;xx=lindgen(npix)
;xx_oversamp=findgen(npix*oversamp)/oversamp
;vector for rebinning
;undersamp_index=(lindgen(npix)-1)*oversamp

;;;Reconstruct wl array including perturbations
x=dindgen(npix)
freeparam_0=replicate(wl_0, npix)
freeparam_1=wl_1*x
freeparam_2=wl_2*(x^2)
freeparam_3=wl_3*(x^3)

trial_soln=freeparam_0+freeparam_1+freeparam_2+freeparam_3


;Apply a trial lsf to the lab spectrum
int_lab_conv=gauss_smooth(int_lab, p[4])

;Adjust the optical depth of the lab spectrum
int_lab_depth=int_lab_conv^tau_scale


;;;Doppler shift
;Shift template wl scale (in observer frame) by trial rv amount
;wl_model=wl_temp_obs_frame*((1.d)+(rad_vel/c))


;Interpolate model onto trial wl grid
int_lab_final=interpol(int_lab_depth, wl_lab, trial_soln)

; try comparing the max depth features
; features_obs=sort(int_trial)
; features_lab=sort(int_lab)
; features_obs_posn=wl_lab[where(features_obs lt 175)]
; features_lab_posn=wl_lab[where(features_lab lt 175)]

; err=0
; chi2=0
; for feature=0, n_elements(features_obs_posn)-1 do begin
;     err=(features_lab_posn[feature]-features_obs_posn[feature])^2
;     chi2=chi2+err
; endfor


;;;Get chi^2
err=0
chi2=0
 for pixel=69, npix-70 do begin  ;cut off 50 pixels on either end
   
     err=((int_obs[pixel]-int_lab_final[pixel])^2);/(error[pixel]^2)
    
     chi2=chi2+err
    
 endfor
        
;;;Do animation of a fitting spectrum, zoomed in on a line at low wl and
;;;lines at high wl
;index1=where(wl_ref gt 1.588d4 and wl_ref lt 1.591d4)
;index2=where(wl_ref gt 1.637d4 and wl_ref lt 1.641d4)
;!p.multi=[0,2,1]
;plot, wl_ref[index1], observation[index1], /xs, yr=[0,1.1] 
;oplot, wl_ref[index1], model_obs[index1], color=200
;plot, wl_ref[index2], observation[index2], /xs, yr=[0,1.1] 
;oplot, wl_ref[index2], model_obs[index2], color=200

;!p.multi=[0,2,2]
; index1=where(wl_ref gt 1.588d4 and wl_ref lt 1.591d4)
; index2=where(wl_ref gt 1.60d4 and wl_ref lt 1.603d4)
; index3=where(wl_ref gt 1.607d4 and wl_ref lt 1.61d4)
; index4=where(wl_ref gt 1.638d4 and wl_ref lt 1.641d4)
; plot, wl_ref[index1], observation[index1], /xs, yr=[0.6,1.1] 
; oplot, wl_ref[index1], model_obs[index1], color=200
; plot, wl_ref[index2], observation[index2], /xs, yr=[0.6,1.1] 
; oplot, wl_ref[index2], model_obs[index2], color=200
; plot, wl_ref[index3], observation[index3], /xs, yr=[0.6,1.1] 
; oplot, wl_ref[index3], model_obs[index3], color=200
; plot, wl_ref[index4], observation[index4], /xs, yr=[0.6,1.1] 
; oplot, wl_ref[index4], model_obs[index4], color=200


;;;Do animation of fitting spectrum for full spectrum
 !p.multi=0
 plot, trial_soln, int_obs , /xs, yr=[0,1.5];, xr=[2.31,2.32]
 oplot, trial_soln, int_lab_final, color=200
 wait, 0.05

print, chi2
return, chi2

end




pro IRCS_cal_old

common amoeba_info, wl_lab, int_lab, int_obs

;read in both spectra
readcol, 'NH3_lab.dat', wl_lab, int_lab, format='D,D'
readcol, 'NH3_obs.dat', trash, int_obs_backwards, format='I,D'

;flip the spectrum to read low to high
int_obs=reverse(int_obs_backwards)

;directly compare them using AMOEBA
;;;;;;;;;;;;;;;;;;;;;;
;;;  Run amoeba   ;;;;
;;;;;;;;;;;;;;;;;;;;;;

;;;Define amoeba inputs
ftol=1e-10
;;;Define entries in guess vector
;define delta wl
del_wl=(wl_lab[1023]-wl_lab[0])/1023d0 

;guess0=0.d                   ;guess for rv
guess1=wl_lab[0]-(25d0*del_wl)                     ;guess for wl_param_0
guess2=0.9d0*del_wl                     ;guess for wl_param_1
guess3=0d                     ;guess for wl_param_2
guess4=0d                     ;guess for wl_param_3
guess5=1d			;guess for gaussian sigma
guess6=1.7d	  		;guess for tau
;guess3=5d-6                     ;guess for spec_shift_0
;guess4=5d-8                     ;guess for spec_shift_1
;;;Define entries in scale vector
;scale0=1000.d
scale1=7*del_wl
scale2=0.1d0*del_wl
scale3=1d-8
scale4=1d-10
scale5=0.5d
scale6=0.3d
;;;Define guess and scale vectors
guess=[guess1,guess2,guess3,guess4,guess5,guess6]
scale=[scale1,scale2,scale3,scale4,scale5,scale6]

!p.multi=[0,2,2]
;;;Run amoeba
r=amoeba3(ftol, scale=scale, p0=guess,function_name='amoebafunction', $
         function_value=fval, nmax=2000)
 if n_elements(r) ne 0 then begin ;if amoeba returns a result
     ;velocity[file_num]=r[0]
     wl_param0=r[0]
     wl_param1=r[1]
     wl_param2=r[2]
     wl_param3=r[3]
     chi2_all=fval[0]
     
 endif else begin ;if amoeba doesn't return a result
     ;velocity=!Values.F_NaN
     wl_param0=!Values.F_NaN
     wl_param1=!Values.F_NaN
     wl_param2=!Values.F_NaN
     wl_param3=!Values.F_NaN
     chi2_all=!Values.F_NaN
     
 endelse

!p.multi=0

;;;Print final results for a given file
;print, "File: ", file_num
;print, "RV: ", velocity[file_num]
print, "WL_0: ", wl_param0
print, "WL_1: ", wl_param1
print, "WL_2: ", wl_param2
print, "WL_3: ", wl_param3
print, "LSF_sigma: ", r[4]
print, "Optical_depth", r[5]
print, "Chi^2: ", fval[0]
stop

;;;Write to file
;output_file='result_'+result_label+'_epoch.fits'
;if epoch_num eq initial_epoch then mwrfits, mjd, output_file, /create $
;else mwrfits, mjd, output_file
;mwrfits, velocity, output_file
;mwrfits, mjd_epoch, output_file
;mwrfits, avg_rv, output_file
;mwrfits, err_bar, output_file
;mwrfits, chi2_all, output_file
;mwrfits, wl_param0, output_file
;mwrfits, wl_param1, output_file 
end
