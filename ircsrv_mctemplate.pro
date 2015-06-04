

pro ircsrv_mctemplate, epoch=epoch, object=object, trace=trace, visualize=visualize, first_pix=first_pix, npix_select=npix_select, mode=mode

npix=1024L
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;MAIN BODY ;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if n_elements(epoch) eq 0 then epoch='18Jan2011'
if n_elements(object) eq 0 then object='GJ273'
if n_elements(trace) eq 0 then trace='AB'


if n_elements(visualize) eq 0 then visualize=1
;Define keywords for fitting just a range
if n_elements(first_pix) eq 0 then first_pix=0
if n_elements(npix_select) eq 0 then npix_select=512L ;npix_select=npix-first_pix
if n_elements(mode) eq 0 then mode='mpfit'


;set file paths
rootpath='/home/stgilhool/RV_projects/IRCS_rv/data/'
epochpath=rootpath+'epoch/'+epoch+'/'
objectpath=epochpath+'final_spectra/'
flatpath=epochpath+'final_spectra/'
calibpath=epochpath+'calib_results/'
outpath=epochpath+'temp_results/'
modelpath=rootpath+'supplemental/'
outputpath=epochpath+'rv_results/'
rvshift1outpath=rootpath+'rvshift1_results/'
templatepath=rootpath+object+'/'

;p=paths(epoch=epoch, object=object, trace=trace)
;rootpath=p.dataroot
;epochpath=p.epochfold
;objectpath=p.specfold
;flatpath=p.specfold
;outpath=p.tempfold
;modelpath=p.dataroot+'supplemental/'


;Define flat file names and full path file names
ABflat_filename=strjoin(['NH3', object, epoch, 'AB1'], '_')+'.fits'
BAflat_filename=strjoin(['NH3', object, epoch, 'BA1'], '_')+'.fits'
ABflat_file=flatpath+strjoin(['NH3', object, epoch, 'AB1'], '_')+'.fits'
BAflat_file=flatpath+strjoin(['NH3', object, epoch, 'BA1'], '_')+'.fits'


;Define object file names and full path file names
ABobj_listname=strjoin([object, epoch, 'AB'], '_')+'.list'
BAobj_listname=strjoin([object, epoch, 'BA'], '_')+'.list'
ABobj_list=objectpath+strjoin([object, epoch, 'AB'], '_')+'.list'
BAobj_list=objectpath+strjoin([object, epoch, 'BA'], '_')+'.list'


;Read in object filenames and then spectra
readcol, ABobj_list, ABobj_filename, format='A'
readcol, BAobj_list, BAobj_filename, format='A'

ABobj_file=objectpath+ABobj_filename
BAobj_file=objectpath+BAobj_filename

n_ABobj=n_elements(ABobj_file)
n_BAobj=n_elements(BAobj_file)


;Initialize arrays of spectra
temp_str=mrdfits(ABobj_file[0], 1)
temp_spectrum=temp_str.spectrum
npix=n_elements(temp_spectrum)
ABspec_arr=dblarr(n_ABobj, npix)
BAspec_arr=dblarr(n_BAobj, npix)
ABflat_arr=dblarr(n_ABobj, npix)
BAflat_arr=dblarr(n_BAobj, npix)

ABmjd_arr=dblarr(n_ABobj)
BAmjd_arr=dblarr(n_BAobj)

ABerr_arr=dblarr(n_ABobj, npix)
BAerr_arr=dblarr(n_BAobj, npix)

;Set continuum normalization params
low_reject=0.53
high_reject=3.0

;Loop through each object file and flip and normalize it
for f=0, n_ABobj-1 do begin
    ABobj_spec=mrdfits(ABobj_file[f], 1)
    BAobj_spec=mrdfits(BAobj_file[f], 1)

    ABspec_backwards=ABobj_spec.spectrum
    BAspec_backwards=BAobj_spec.spectrum

    ;Flip the spectra
    ABspec=reverse(ABspec_backwards)
    BAspec=reverse(BAspec_backwards)

    ;Read in and flip errors
    ABerr=reverse(ABobj_spec.sigma)
    BAerr=reverse(BAobj_spec.sigma)

    ;Get MJD
    ABhead=ABobj_spec.header
    BAhead=BAobj_spec.header
    
    ;Normalize the model spectrum and the observed spectrum
    ABnorm=continuum_fit(dindgen(npix), ABspec, low_rej=low_reject, high_rej=high_reject)
    BAnorm=continuum_fit(dindgen(npix), BAspec, low_rej=low_reject, high_rej=high_reject)
    ;Populate arrays
    ;spectra
    ABspec_arr[f,*]=ABspec/ABnorm
    BAspec_arr[f,*]=BAspec/BAnorm
    ;errors
    ABerr_arr[f,*]=ABerr/ABnorm
    BAerr_arr[f,*]=BAerr/BAnorm
    ;MJD
    ABmjd_arr[f]=sxpar(ABhead, 'MJD')
    BAmjd_arr[f]=sxpar(BAhead, 'MJD')

endfor

;;;;;;;

;Read in calibration results
calib_file='GJ273_18Jan2011_AB1_5_7_7_13.fits'
calib_ext=13
;model_file=calibpath+calib_file
model_file=rootpath+'epoch/18Jan2011/calib_results/'+calib_file
model_par=mrdfits(model_file, calib_ext)

common modelinfo, delta_rv_index, h2o_depth_index, co2ch4_depth_index, delta_wl_index, gh0_coeff_index, gh1_coeff_index, other_index, lin_switch, wl_telluric, h2o, co2ch4, npixels, int_lab, wl_lab, temp_spec, temp_wl, oversamp, npixselect, firstpix, int_obs, err, fmode, visit, last_guess, visual, wl_soln, wl_soln_select, wl_soln_over, wl_soln_over_select, x, xx, x_select, xx_select, bcv, delta_bcv

npixels=1024L
visual=visualize
npixselect=npix_select
firstpix=first_pix
fmode=mode

;Assume lin_switch is on
lin_switch=1

;;;Parameters from earlier calibration
;wl
wl_coeff=model_par.wl_result
wl_scale=model_par.wl_scale


;;;

;Other parameters and constants
oversamp=model_par.oversamp
npix_lsf=(oversamp*10L)+1L
bigc=299792458.D
npix_model=npix_select*oversamp
npix_over=npix*oversamp


;Read in LAB SPECTRUM
modelfile= modelpath+'NH3_model.dat'
readcol, modelfile, wl_lab, int_lab, format='D,D'

norm_lab=continuum_fit(wl_lab, int_lab, low_rej=low_reject, high_rej=high_reject)
int_lab_copy=int_lab
int_lab=int_lab_copy/norm_lab



;;;Contstruct wl array and oversampled wl array from coeffs
x=dindgen(npix)
xx=(dindgen(npix*oversamp)-(oversamp/2))/oversamp

wl_soln=poly(x, wl_coeff)
wl_soln_over=poly(xx, wl_coeff)

;Redefine those vectors to reflect only the specified range
x_select=dindgen(npix_select)+first_pix
xx_select=(dindgen(npix_model)-(oversamp/2))/oversamp+first_pix

wl_soln_select=poly(x_select, wl_coeff)
wl_soln_over_select=poly(xx_select, wl_coeff)


;Read in Template
templatefile=templatepath+object+'_template.dat'
readcol, templatefile, temp_wl, temp_spec, format='D,D'


;oversample the template and template_wl
if n_elements(temp_wl) eq npix then begin
    temp_wl_over=interpol(temp_wl, x, xx)
    template_over=interpol(temp_spec, temp_wl, temp_wl_over)
endif else begin
    temp_wl_over=temp_wl
    template_over=temp_spec
endelse


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


;FIX TO MAKE MORE ELEGANT: BCV of template
stru=mrdfits(ABobj_file[0], 1)
head1=stru.header
;get bcv correction
bcvcorr_ircs, head1, params1
bcv0=1000d0 * params1[0]


;Define some arrays
n_exp=n_ABobj-2
lsf=dblarr(npix_lsf, npix_model, n_exp)
tell_nh3_spectrum=dblarr(npix_over, n_exp)
stellar=dblarr(npix_model, n_exp)
wl_grid_over=dblarr(npix_over, n_exp)
err=dblarr(npix, n_exp)
final_spectrum=dblarr(npix_model)
delta_rv=dblarr(n_exp)
wl_grid_over=dblarr(npix_over, n_exp)
wl_grid_down=dblarr(npix, n_exp)
wl_template_shifted=dblarr(npix_over, n_exp)
norm_factor=dblarr(npix_select, n_exp)
obs=dblarr(npix, n_exp)

;Now fit each observation

;for visit=0, n_ABobj-1 do begin


;if visit eq 5 or visit eq 7 then continue

varray=[0,1,2,3,4,6,8,9,10,11,12,13]

foreach visit, varray, index do begin

;Make observation array
    obs[*, index]=ABerr_arr[visit, *]
;get error for this observation
    err[*, index]=ABerr_arr[visit, *]
;Mask entries
;Cut off pixels on either end of full 1024 pixel spectrum
    npix_trim_start=0L
    npix_trim_end=10L
    bigerr=1d10                 ;error value for masked pixels
    if n_elements(err) eq 1024L then begin
        ;err[0:npix_trim_start-1, index]=bigerr
        ;err[-1L*npix_trim_end:-1, index]=bigerr
    endif



;READ IN PARAMETERS FROM 2-basis 2-run FITS
    rfile='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_127_run2sign.fits'

    f=mrdfits(rfile, visit+1)
    re=f.result

    delta_rv[index]=re[f.delta_rv_index]
    delta_wl_coeff=re[f.delta_wl_index]
    h2o_depth=re[f.h2o_depth_index]
    co2ch4_depth=re[f.co2ch4_depth_index]
    delta_wl_coeff=re[f.delta_wl_index]
    gh0_coeff=re[f.gh0_coeff_index]
    gh1_coeff=re[f.gh1_coeff_index]
    other=re[f.other_index]

    tau_scale=other[0]	


;;;BEGIN MODELING

;Perturb wl soln
    delta_wl=poly(xx, delta_wl_coeff)
    wl_grid_over[*, index]=wl_soln_over+delta_wl[index]

    delta_wl_down=poly(x, delta_wl_coeff)
    wl_grid_down[*, index]=wl_soln+delta_wl_down


;NH3 Lab Spectrum: Interpolate onto the oversampled trial grid
;adjust optical depth first
    int_lab_depth=int_lab^tau_scale
    ammonia=interpol(int_lab_depth, wl_lab, wl_grid_over[*, index])




;Telluric: Adjust optical depths and construct telluric spectrum
    h2o_scaled=h2o^h2o_depth
    co2ch4_scaled=co2ch4^co2ch4_depth
    telluric_long=h2o_scaled*co2ch4_scaled
;put telluric spectrum onto grid
    
    telluric=interpol(telluric_long, wl_telluric, wl_grid_over[*, index])


;Product of the three
    tell_nh3_spectrum[*,index] = ammonia * telluric


;Make LSF

    lsf[*,*,index] = ircsrv_lsf(gh0_coeff, gh1_vec=gh1_coeff, oversamp=oversamp, first_pix=first_pix, npix_select=npix_select, neg_penalty=penalty)

;Shift stellar wl
    wl_template_shifted[*, index]=temp_wl_over*(1d0 + delta_rv[index]/bigc)


;Correct small normalization errors by multiplying smooth function
    if n_elements(other) gt 1 then begin
        nparam_other=n_elements(other)
        nparam_norm=nparam_other-1
        norm_pts_y=other[1:nparam_norm]
        
                                ;X vector of same dimensions as ircs-resolution spectrum
        norm_pts_xx=lindgen(npix_select)+first_pix
        
                                ;Create and populate vector with x-coord of node pts
        norm_pts_x=dblarr(nparam_norm)
        
        norm_pts_x[0]=first_pix
        norm_pts_x[nparam_norm-1]=first_pix+npix_select-1
        x_increment=npix_select/(nparam_norm-2)
        for i=1, nparam_norm-2 do norm_pts_x[i]=norm_pts_x[i-1]+x_increment
        
                                ;Interpolate between nodes with a spline
        norm_factor[*,index]=interpol(norm_pts_y, norm_pts_x, norm_pts_xx, /spline)    
        
    endif

endforeach

                                ;Truncate observations and errors to model length
                                ;Indices that we are using
select_x = lindgen(npix)
select_index = select_x[first_pix:first_pix+npix_select-1]
select_xx= lindgen(npix_over)
select_index_over = select_xx[first_pix*oversamp:(first_pix*oversamp)+npix_model-1]

temp_display=template_over[select_index_over]






;NOW DO THE BILLIONS OF TRIALS

n_trials=1e2
;maxiter=1e9 - 1


                                ;Make 3D arrays out of some of the
                                ;arrays
wl_template_shifted_arr=rebin(wl_template_shifted, npix_over, n_exp, n_trials)
wl_grid_over_arr=rebin(wl_grid_over, npix_over, n_exp, n_trials)
tell_nh3_spectrum_arr=rebin(tell_nh3_spectrum, npix_over, n_exp, n_trials)
norm_factor_arr=rebin(norm_factor, npix_select, n_exp, n_trials)    
obs_select=rebin(obs[select_index, *], npix_select, n_exp, n_trials)
err_select=rebin(err[select_index, *], npix_select, n_exp, n_trials)
;lsf_arr=rebin(lsf, npix_lsf, npix_model, n_exp, n_trials)
lsf_arr=temporary(lsf)
help, lsf
help, lsf_arr

;Prepare other arrays
final_template=dblarr(npix_model)
chi2_arr=dblarr(npix_select, n_trials)

win_width_down=3L
win_width_over=win_width_down * oversamp
n_windows=npix_select-win_width_down
min_arr=dblarr(2, n_windows+1)    

bestfile= '/home/stgilhool/RV_projects/IRCS_rv/data/GJ273/mctemplate.dat'

if file_test(bestfile) then begin
    templatestr=mrdfits(bestfile, 1)
    n_trials_completed=templatestr.n_trials
    final_template=templatestr.template
    min_arr=templatestr.min_arr
endif else n_trials_completed=0L





while n_trials_completed lt 1e9 do begin
    clock=tic()

    ;Graph Template
    plot, temp_display, /xs, yr=[0.2,1.3]
    wait,1


    

    ;Make big data arrays for template and shifts
    new_template=template_over
    new_template[0:49]=1d0
    template_arr=rebin(new_template, npix_over, n_exp, n_trials)
    temp_shift=randomn(seed, npix_over, n_trials)*0.02d0
    temp_shift_arr=rebin(reform(temp_shift, npix_over, 1, n_trials), npix_over, n_exp, n_trials)

    template_arr=temporary(template_arr)+temp_shift_arr ;this is our baby
    template_shift_arr=0
    template_trial=temporary(template_arr)

    ;stellar_trial=dblarr(npix_over, n_exp, n_trials) 
    model_over=dblarr(npix_model, n_exp, n_trials)
    ;model=dblarr(npix_select, n_exp, n_trials)
    
    
    
                                ;Stellar Template: RV shift the stellar template
    stellar_trial=interpol(template_trial, wl_template_shifted_arr, wl_grid_over_arr)
    
                                ;Now convolve to make full model
    product_spectrum=tell_nh3_spectrum_arr * stellar_trial
    
        
            
;            product=product_spectrum[*,index]
    for exp=0, n_exp-1 do begin
        for trial=0, n_trials-1 do begin
        ;convclockt=ic()
        model_over[*,exp,trial]= ircsrv_convolve(product_spectrum[*,exp,trial], lsf_arr[*,*,exp], oversamp=oversamp, first_pix=first_pix, npix_select=npix_select)
        ;toc, convclock
    endfor
    endfor
        
                                ;Downsample to IRCS resolution
    model = ircsrv_downsample(model_over, oversamp)
        
        
                                ;final IRCS sampled 2d array of models
                                ;for all exposures (n_exp x npix_select)
    model_final=temporary(model)*norm_factor_arr
    
        
                
                ;Calculate chi2
    chi2_allexp=((model_final-obs_select)/err_select)^2
    chi2_arr=total(chi2_allexp, 2, /double)

	;Select best template by windowing
	;Let's start with a three IRCS pixel window
  

    for l_pix=0, n_windows do begin

        window=lindgen(win_width_down)+l_pix
        window_over=lindgen(win_width_over)+(l_pix*oversamp)
        min_err=min(total(chi2_arr[window, *],1, /double), best_trial_index)
            
        if n_trials_completed eq 0 then begin
            min_arr[*,l_pix]=[min_err, best_trial_index]
            final_template[window_over]=template_trial[window_over, best_trial_index]
            oplot, final_template, color=200, ps=3
            
        endif else begin
                                ;check if new chunk is better than old chunk
            if min_err lt min_arr[0,l_pix] then begin
                min_arr[*,l_pix]=[min_err, best_trial_index]
                final_template[window_over]=template_trial[window_over, best_trial_index]
                oplot, final_template, color=200, ps=3
            endif
        endelse
    endfor
    wait, 1

n_trials_completed=n_trials+n_trials_completed

print, strtrim(n_trials_completed,2)+" trials completed"

outstr={n_trials:n_trials_completed , $
        template:final_template, $
        min_arr:min_arr, $
        bigchi:total(min_arr[0,*]) $
        }

mwrfits, outstr, bestfile, /create
toc, clock
endwhile



stop


end
