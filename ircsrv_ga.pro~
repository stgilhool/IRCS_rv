function rvmodel, p, fakekey=fakekey

common modelinfo, delta_rv_index, h2o_depth_index, co2ch4_depth_index, delta_wl_index, gh0_coeff_index, gh1_coeff_index, other_index, lin_switch, wl_telluric, h2o, co2ch4, npix, int_lab, wl_lab, template, wl_template, oversamp, npix_select, first_pix, observation, err, mode, visit, last_guess, visualize,wl_soln, wl_soln_select, wl_soln_over, wl_soln_over_select, x, xx, x_select, xx_select, bcv, delta_bcv, parscale_all, wl_start, wl_index, sig_start, rv_index, param_m, param_b,ghga_coeff_index, otherga_index

;Record parameters for posterity
last_guess=p

;rescale
p_scaled=p/parscale_all



;Read in parameters
delta_rv=p_scaled[delta_rv_index]
delta_wl_coeff=p_scaled[delta_wl_index]
gh0_coeff=p_scaled[gh0_coeff_index]
if lin_switch eq 1 then gh1_coeff=p_scaled[gh1_coeff_index]
h2o_depth=p_scaled[h2o_depth_index]
co2ch4_depth=p_scaled[co2ch4_depth_index]
other=p_scaled[other_index]

tau_scale=other[0]

;Constants
bigc=299792458d0 ;m/s
bigerr=1d10

;;;BEGIN MODELING

;Define wl grid
;x=dindgen(npix)
;xx=(dindgen(npix*oversamp)-(oversamp/2))/oversamp

;wl_grid=poly(x, wl_coeff)
;wl_grid_over=poly(xx, wl_coeff)

;Perturb wl soln
delta_wl=poly(xx, delta_wl_coeff)
wl_grid_over=wl_soln_over+delta_wl

delta_wl_down=poly(x, delta_wl_coeff)
wl_grid_down=wl_soln+delta_wl_down


;NH3 Lab Spectrum: Interpolate onto the oversampled trial grid
;adjust optical depth first
int_lab_depth=int_lab^tau_scale
ammonia=interpol(int_lab_depth, wl_lab, wl_grid_over)


;Stellar Template: RV shift the stellar template
wl_template_shifted=wl_template*(1d0 + (1000d0*delta_rv)/bigc)
stellar=interpol(template, wl_template_shifted, wl_grid_over)


;Telluric: Adjust optical depths and construct telluric spectrum
h2o_scaled=h2o^h2o_depth
co2ch4_scaled=co2ch4^co2ch4_depth
telluric_long=h2o_scaled*co2ch4_scaled
;put telluric spectrum onto grid
telluric=interpol(telluric_long, wl_telluric, wl_grid_over)


;Product of the three
product_spectrum = ammonia * telluric * stellar

;stop

;Take a detour to make LSF
if lin_switch eq 1 then begin
    lsf = ircsrv_lsf(gh0_coeff, gh1_vec=gh1_coeff, oversamp=oversamp, first_pix=first_pix, npix_select=npix_select, neg_penalty=lsf_penalty)
endif else $
  lsf = ircsrv_lsf(gh0_coeff, oversamp=oversamp, first_pix=first_pix, npix_select=npix_select, neg_penalty=lsf_penalty)




;Now convolve to make full model
model_over = ircs_convolve(product_spectrum, lsf, oversamp=oversamp, first_pix=first_pix, npix_select=npix_select)


;Downsample to IRCS resolution
model = downsample_tophat(model_over, oversamp)


;Correct small normalization errors by multiplying smooth function
;if n_elements(other) eq 2 then begin
;    norm_factor=

if n_elements(other) gt 1 then begin
    nparam_other=n_elements(other)
    nparam_norm=nparam_other-1
    norm_pts_y=other[1:nparam_norm]

    ;X vector of same dimensions as ircs-resolution spectrum
    norm_pts_xx=lindgen(npix_select)+first_pix

    ;Create and populate vector with x-coord of node pts
    norm_pts_x=dblarr(nparam_norm)

    if nparam_norm ge 2 then begin
        norm_pts_x[0]=first_pix
        norm_pts_x[nparam_norm-1]=first_pix+npix_select-1
        if nparam_norm ge 3 then begin
            x_increment=npix_select/(nparam_norm-2)
            for i=1, nparam_norm-2 do norm_pts_x[i]=norm_pts_x[i-1]+x_increment
        endif
    
    ;Interpolate between nodes with a spline
        if npix_select ge 128L and nparam_norm ge 5 then norm_factor=interpol(norm_pts_y, norm_pts_x, norm_pts_xx, /spline) $
        else norm_factor=interpol(norm_pts_y, norm_pts_x, norm_pts_xx)
    endif else if nparam_norm eq 1 then norm_factor=replicate(norm_pts_y, npix_select)

    ;Correct model
    model = model * norm_factor

endif

    

;Select appropriate parts of all vectors
;Indices that we are using
select_x = lindgen(npix)
select_index = select_x[first_pix:first_pix+npix_select-1]
;Check that this equals npix_model
if n_elements(select_index) ne n_elements(model) then message, "something's wrong with select"

;Select appropriate parts of all vectors
;downsample oversampled vectors
stellar_down=downsample_tophat(stellar, oversamp)
telluric_down=downsample_tophat(telluric, oversamp)
ammonia_down=downsample_tophat(ammonia, oversamp)
;select
stellar_select=stellar_down[select_index]
telluric_select=telluric_down[select_index]
ammonia_select=ammonia_down[select_index]
observation_select = observation[select_index]
err_select = err[select_index]
model_select=model
wl_grid_select=wl_grid_down[select_index]


;Fix model entries that are NaNs, if necessary
nan=where(finite(model) eq 0, nancount)
if nancount gt 0 then begin
    ;err_select[nan]=bigerr
    ;nan_penalty=10d8*nancount
    nan_penalty=0.1d0*nancount
    model_select[nan]=1d0
endif else nan_penalty=0


;;;
;Calculate CHI2
;;;

if lsf_penalty gt 1d-7 then lsf_penalty= 1d0+lsf_penalty

penalty=1d0-exp(-1d0*(lsf_penalty+nan_penalty))

residuals=(observation_select - model_select)

deviation = residuals/err_select
;deviation_squared = deviation^2


chi1_vec=deviation * exp(penalty)

chi1_vec_nopen=deviation



chi2=total(chi1_vec^2, /double, /nan)
chi2_nopen=total(chi1_vec_nopen^2, /double, /nan)

;;;;
;If visualize is set, do animation of fitting spectrum for full spectrum
;;;;



if (visualize eq 1 ) then begin ;or mode eq 'one_call') then begin

    ;Set up big dots
    phi=findgen(32)*(!PI*2/32.)
    phi = [ phi, phi(0) ]
    usersym, cos(phi), sin(phi), /fill


    title_str="RV shift: ["+strtrim(delta_rv,2) + ", " + strtrim(delta_rv+delta_bcv,2) +"] | visit: " + strtrim(visit,2)+ " | chi2: "+strtrim(chi2,2) + " | chi2_nopen: " + strtrim(chi2_nopen,2) + " | D_WL : [" + strtrim(wl_start,2) + ", " + strtrim(delta_wl_coeff,2) + "] | SIG : [" + strtrim(sig_start,2)

    plot, wl_grid_select, ammonia_select,yr=[0.3,1.1], /xs, title=title_str, ytitle='NH3', charsize=1.5
    plot, wl_grid_select, stellar_select,yr=[0.3,1.1], /xs, ytitle='Star'

    plot, wl_grid_select, telluric_select, yr=[0.3, 1.1], /xs, ytitle='Tell'
    
    plot, wl_grid_select, observation_select, yr=[0.3, 1.1], /xs
    oplot, wl_grid_select, model_select, ps=8, color=200, symsize=0.5

    plot, wl_grid_select, residuals, xtitle="Wavelength (microns)", yr=[-0.1, 0.1], ps=3, /xs
    
endif



;;;;
;Report results
;;;;


print, chi2, chi2_nopen
print, "D_WL: [" + strtrim(wl_start,2) + ", " + strtrim(delta_wl_coeff,2) + "]" 
print, "SIGMA: [" + strtrim(sig_start,2) + ", " + strtrim(gh0_coeff[0],2) + "]" 

print, "RV shift: ["+strtrim(delta_rv,2) + ", " + strtrim(delta_rv+delta_bcv,2) +"] "
print, "----------"

case mode of

    'ga': begin
        return, chi2
    end

    'mpfit': begin
        return, chi1_vec
    end

    
    'tnmin': begin
        return, chi2
    end


    'amoeba': begin
        return, chi2
    end


    'one_call': begin
        model_output={rv: delta_rv, $
                      bcv: bcv, $
                      delta_bcv:delta_bcv, $
                      params: p, $
                      h2o_depth:p[h2o_depth_index], $
                      co2ch4_depth:p[co2ch4_depth_index], $
                      gh0_coeff:p[gh0_coeff_index], $
                      delta_wl:p[delta_wl_index], $
                      other:p[other_index], $
                      wl_grid: wl_grid_select, $
                      nh3: ammonia_select, $
                      stellar: stellar_select, $
                      telluric: telluric_select, $
                      lsf: lsf, $
                      model: model_select, $                      
                      obs: observation_select, $
                      residuals: residuals, $
                      error:err_select, $
                      chi2: chi2}
        return, model_output
    end
       
    else: message, "appropriate mode not specified"
endcase

end


function ga_model, p, np, funa=funa

common modelinfo, delta_rv_index, h2o_depth_index, co2ch4_depth_index, delta_wl_index, gh0_coeff_index, gh1_coeff_index, other_index, lin_switch, wl_telluric, h2o, co2ch4, npixels, int_lab, wl_lab, template_over, temp_wl_over, oversamp, npixselect, firstpix, int_obs, err, fmode, visit, last_guess, visual, wl_soln, wl_soln_select, wl_soln_over, wl_soln_over_select, x, xx, x_select, xx_select, bcv, delta_bcv, parscale_all, wl_start, wl_index, sig_start, rv_index, param_m, param_b, ghga_coeff_index, otherga_index

;fmode='ga'
ga_gen=tic()

param_m=funa.param_m
param_b=funa.param_b

;p_scaled=param_m*p+param_b



chi2_vec=dblarr(np)
for model_num=0, np-1 do begin
    p_scaled=param_m*p[*,model_num]+param_b
    chi2=rvmodel(p_scaled)
    if chi2 lt 1 then chi2=1d10
    chi2_vec[model_num]=chi2
endfor

help, chi2_vec
print, minmax(chi2_vec)
;stop

ga_gen_time=toc(ga_gen)
;print, 'one generation of ' + strtrim(np,2) + ' models takes ' + strtrim(ga_gen_time,2) + ' seconds.'
;stop

return, chi2_vec

end



pro ircsrv_ga, epoch=epoch, object=object, trace=trace, visualize=visualize, first_pix=first_pix, npix_select=npix_select, mode=mode, run=run, n_bases_lsf=n_bases_lsf, n_other=n_other, s_iter=s_iter, current_tag=current_tag, initial_file=initial_file, s_template_num=s_template_num, n_lsf_ga=n_lsf_ga, _extra=ex

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
if n_elements(npix_select) eq 0 then npix_select=128L;npix_select=npix-first_pix
if n_elements(mode) eq 0 then mode='ga'


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




;;;Tuned variables

;Iteration number: 0 - initial run, 1 - run that fits for RV    
if n_elements(run) eq 0 then run=0  

;Iteration number of stellar template    
if n_elements(s_iter) eq 0 then s_iter=0
	n_run=run+(s_iter*2) ;The total number of runs (used to label files)


initial_tag='sign' ;Tag for reading in guesses for very first run
;current_tag='Mar16bigpen' ;Tag for marking files in this set of iterations
if n_elements(current_tag) eq 0 then current_tag='Mar17'

	;Setting up filenames and paths (basically not tuned)
	model_id_base=strjoin([object, epoch,'AB', strtrim(first_pix,2), strtrim(first_pix+npix_select-1,2)], '_')
        model_id_full=strjoin([model_id_base, current_tag ,'r'+strtrim(n_run,2)], '_')
        file_base=rvshift1outpath + model_id_base

;Controls RV guess: 0 - Start with just bcv-based guess, 1 - Read in guess from result of previous run
if run eq 0 and s_iter eq 0 then rv_read=0 else rv_read=1

;Which input file
if run eq 0 and s_iter eq 0 then begin
    ;input_file=strjoin([file_base, initial_tag, 'r1'], '_') + '.fits' 
    input_file=strjoin([file_base, initial_tag], '_') + '.fits'
    input_file=rvshift1outpath + 'GJ273_18Jan2011_AB_256_383_Mar20_mpfit2_g1pen_r0.fits'

;'GJ273_18Jan2011_AB_0_127_Mar16bigpen_r0.fits'
    if n_elements(initial_file) ne 0 then input_file=initial_file
endif else begin
    input_file=strjoin([file_base, current_tag, 'r'+strtrim(n_run-1,2)], '_') + '.fits'
endelse

;Output file
;output_file=strjoin([file_base, current_tag, 'r'+strtrim(n_run,2)], '_') + '.fits'
output_file=strjoin([file_base, current_tag, 'mpfit'], '_') + '.fits'

;Stellar Template
;s_template=s_iter ;Controls which stellar template to use (0 for
;original iterative one, 1 for 7x oversamp fitted one) 
if n_elements(s_template_num) eq 0 then s_template_num=0L


if mode eq 'one_call' then begin ;This allows us to call one_call on a run without switching the I/O tags
    oc_output_file=output_file+'_spec.fits'
endif

telluric_option=0 ;0-original model, 1-telluric_test model, 2-h2o_test model

if n_elements(n_bases_lsf) eq 0 then n_bases_lsf=2L

if n_elements(n_other) eq 0 then n_other=13L


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

common modelinfo, delta_rv_index, h2o_depth_index, co2ch4_depth_index, delta_wl_index, gh0_coeff_index, gh1_coeff_index, other_index, lin_switch, wl_telluric, h2o, co2ch4, npixels, int_lab, wl_lab, template_over, temp_wl_over, oversamp, npixselect, firstpix, int_obs, err, fmode, visit, last_guess, visual, wl_soln, wl_soln_select, wl_soln_over, wl_soln_over_select, x, xx, x_select, xx_select, bcv, delta_bcv, parscale_all, wl_start, wl_index, sig_start, rv_index, param_m, param_b, ghga_coeff_index, otherga_index

npixels=1024L
visual=visualize
npixselect=npix_select
firstpix=first_pix
fmode=mode

;Assume lin_switch is off
lin_switch=0

;;;Parameters from earlier calibration

wl_coeff=model_par.wl_result
wl_scale=model_par.wl_scale

other_guess=model_par.other_result
other_scale=model_par.other_scale

;;;

;Other parameters and constants
oversamp=model_par.oversamp
npix_lsf=(oversamp*10L)+1L
bigc=299792458.D
npix_model=npix_select*oversamp


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

if s_template_num eq 0 then begin
    templatefile=templatepath+object+'_template.dat'
    readcol, templatefile, temp_wl, temp_spec, format='D,D'
    if n_elements(temp_wl) ne 7L*npix then message, "You were wrong about the old template sampling"
    if n_elements(temp_spec) ne 7L*npix then message, "You were wrong about the old template sampling"
    temp_oversamp=7L
    norm_start=(10*temp_oversamp) > (first_pix*temp_oversamp)
    norm_end= (first_pix*temp_oversamp + npix_select*temp_oversamp - 1) < ((npix-10)*temp_oversamp-1)
    tempnorm=(max(temp_spec[norm_start:norm_end]))
    temp_spec=temp_spec/tempnorm
;FIX the normalization

    ;print, max(temp_spec[norm_start:norm_end])
    openw, lun, templatepath+object+'_template'+current_tag+'_'+strtrim(n_run,2)+'.dat', /get_lun
    for line=0,n_elements(temp_spec)-1 do begin
        printf, lun, temp_wl[line], temp_spec[line], format='(D,D)'
    endfor
    close, lun
    free_lun, lun
    
    ;stop
endif else begin
    ;templatefile='/home/stgilhool/RV_projects/IRCS_rv/data/smooth_penalty_test/test16/penaltytest.fits'
    ;temp_ext=30 ;THIS IS THE OLD 3x OVERSAMP TEMPLATE
    ;temp_oversamp=3L
    ;if s_template eq 1 then templatefile='/home/stgilhool/RV_projects/IRCS_rv/data/smooth_penalty_test/test22/penaltytest.fits'  ;USE THIS ONE, IT'S BETTER
    ;if s_template eq 2 then
    ;templatefile='/home/stgilhool/RV_projects/IRCS_rv/data/smooth_penalty_test/test21/penaltytest.fits'  ;USE THIS ONE, IT'S BETTER
    templatefile='/home/stgilhool/RV_projects/IRCS_rv/data/smooth_penalty_test/test'+strtrim(s_template_num, 2)+'/penaltytest.fits'
    temp_ext=1
    templatestr=mrdfits(templatefile, temp_ext)
    temp_spec=templatestr.template_spec
    temp_wl=templatestr.template_wl
    temp_oversamp=7L
    tempnorm=1D0
endelse


if oversamp ne temp_oversamp then begin 
    xx_temp=(dindgen(npix*temp_oversamp)-(temp_oversamp/2))/temp_oversamp
    temp_wl_over=interpol(temp_wl, xx_temp, xx)
    template_over=interpol(temp_spec, temp_wl, temp_wl_over)

endif else begin
    template_over=temp_spec
    temp_wl_over=temp_wl
endelse

template_over_select=template_over[first_pix*oversamp:(first_pix+npix_select)*oversamp-1]


;;;-----Finished reading in template----


;READ IN TELLURIC



case telluric_option of
    
    0: begin ;This is the orginal, highly-sampled theoretical model
        h2ostr=mrdfits(modelpath+'sky_h2o.fits',1)
        co2ch4str=mrdfits(modelpath+'sky_co2andch4.fits', 1)
        
        wl_telluric_long=h2ostr.wave
        h2o_long=h2ostr.trans
        co2ch4_long=co2ch4str.trans
        
        m24index=where(wl_telluric_long gt 2.275 and wl_telluric_long lt 3.365)
        wl_telluric=wl_telluric_long[m24index]
        h2o=h2o_long[m24index]
        co2ch4=co2ch4_long[m24index]
    end

    1: begin ;This is the model that resulted from the telluricfit code (3x oversamp, h2o free, co2ch4 free)
        tell_file='/home/stgilhool/RV_projects/IRCS_rv/data/telluric_test/test5/tellurictest.fits'
        tell_ext=1
        tellstr=mrdfits(tell_file, tell_ext)
        h2o=tellstr.h2o
        co2ch4=tellstr.co2ch4
        wl_telluric=tellstr.wl_telluric
    end

    2: begin ;This model is from h2ofit code (7x oversamp, h2o free, co2ch4 fixed)
        tell_file='/home/stgilhool/RV_projects/IRCS_rv/data/h2o_test/test2/h2otest.fits'
        tell_ext=3
        tellstr=mrdfits(tell_file, tell_ext)
        h2o=tellstr.h2o
        co2ch4=tellstr.co2ch4
        wl_telluric=tellstr.wl_telluric        
    end

endcase
;;;-----Finished reading in smoothed telluric

;FIX TO MAKE MORE ELEGANT BCV of template
stru=mrdfits(ABobj_file[0], 1)
head1=stru.header
;get bcv correction
bcvcorr_ircs, head1, params1
;bcv0=1000d0 * params1[0]
bcv0=params1[0]

;Now fit each observation
if visualize eq 1 then window, 0, xsize=1500, ysize=1000
!p.multi=[0,1,5]


;for visit=0, n_ABobj-1 do begin
for visit=0,0 do begin

    ;get spectrum for given visit
    int_obs=ABspec_arr[visit,*]
    mjd=ABmjd_arr[visit]
    
    ;get error for this observation
    err=ABerr_arr[visit,*]

    ;Mask entries
    ;Cut off pixels on either end of full 1024 pixel spectrum
    npix_trim_start=10L
    npix_trim_end=10L
    bigerr=1d10                 ;error value for masked pixels
    if n_elements(err) eq 1024L then begin
        err[0:npix_trim_start-1]=bigerr
        err[-1L*npix_trim_end:-1]=bigerr
    endif
    
    ;Set error high for masked pixels
    if n_elements(mask) ne 0 then begin
        if n_elements(mask) ne n_elements(err) then $
          message, "mask must have same length as error vector" $
        else begin
            mask_index=where(mask eq 1, flag_count)
            if flag_count gt 0 then err[mask_index]=bigerr
        endelse
    endif
    
    ;Correct for BCV
    ;get header
    struct=mrdfits(ABobj_file[visit], 1)
    head=struct.header
    ;Get bcv correction
    bcvcorr_ircs, head, params
    bcv=params[0]
    delta_bcv=bcv-bcv0
    
    ;directly compare them using AMOEBA
    ;;;;;;;;;;;;;;;;;;;;;;
    ;;;  Run Modeling Function   ;;;;
    ;;;;;;;;;;;;;;;;;;;;;;
    ;FIX ADD KEYWORD

    
    clock=tic()
    
    ;Define inputs to modeling function
    ftol=1d-10

    rv_guess=[-1*delta_bcv]
    h2o_depth_guess=[1d0]
    co2ch4_depth_guess=[1d0]
    delta_wl_guess=[1d-6]

    gh0_guess=replicate(1d-1, n_bases_lsf)


    other_guess=replicate(1d0, n_other)


    gh0_scale=dblarr(n_bases_lsf)
    rv_scale=[5d4]
    h2o_depth_scale=[0.5d0]
    co2ch4_depth_scale=[0.5d0]
    delta_wl_scale=[5d-6]
    other_scale=dblarr(n_other)


    ;indices
    delta_rv_index=0
    h2o_depth_index=1
    co2ch4_depth_index=2
    delta_wl_index=lindgen(n_elements(delta_wl_guess))+co2ch4_depth_index+1
    gh0_coeff_index=lindgen(n_elements(gh0_guess))+delta_wl_index[-1]+1
    gh1_coeff_index=!values.d_nan
    other_index=lindgen(n_elements(other_guess))+gh0_coeff_index[-1]+1
    nh3_depth_index=other_index[0]

    if n_elements(other_index) gt 1 then $
      norm_index=other_index[1:*]
    nparam_total=other_index[-1]+1
    

    if file_test(input_file) then begin    
        str=mrdfits(input_file, visit+1)
        if rv_read eq 0 then rv_guess=[-1d0*delta_bcv] else rv_guess=str.delta_rv
        n_bases_lsf_in=n_elements(str.gh0_coeff_index)
        
        h2o_depth_guess=str.result[str.h2o_depth_index]
        co2ch4_depth_guess=str.result[str.co2ch4_depth_index]
        
        other_guess_all=str.result[str.other_index]
        other_guess=other_guess_all[0:n_other-1]
    endif else begin
        print, "File doesn't exist"
        openw, lun, '/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/errors.dat', /append, /get_lun
        printf, lun, input_file + ' does not exist. Possible problem.'
        close, lun
        free_lun, lun

        rv_guess=[-1d0*delta_bcv]
        other_guess_all=replicate(1d0, 13)
        other_guess=other_guess_all[0:n_other-1]
    endelse
        

                                ;; SCALE THE PARAMETERS TO BE OF SAMEISH ORDER
       
    
    ;Readin best wl and sig from before
    readcol, 'startingparams.txt', wl_best_list, sig_best_start, chibest, rv_best_list, format='(D,D,D,D)'
    sig_start=sig_best_start[visit]
    wl_start=wl_best_list[visit]
        

    
            
                       
            ;Define mpfit constraints
            parinfo = replicate({parname:'null' , value:0d0 , step:0d0 , relstep:0d0 ,fixed:0, limited:[0,0], $
                                 limits:[0.D,0.D], mpside:0}, nparam_total)
            
  ;;; --- DELTA_RV ---
            ;Name
            parinfo[delta_rv_index].parname='delta_rv' ;in km/s
            
            ;Scaling
            delta_rv_parscale=1d0
            
            ;Value
            parinfo[delta_rv_index].value = rv_guess[0]
            
            ;Step
            parinfo[delta_rv_index].step = 1d-5 ;Step size is 0.01m/s
            
            ;Contraints
            if run eq 1 then begin
                parinfo[delta_rv_index].limited = [1,1]
                if rv_guess[0] eq 0 then rv_guess[0]=1d-3
                parinfo[delta_rv_index].limits = [rv_guess-1d0, rv_guess+1d0]
                parinfo[delta_rv_index].mpside = 2    
            endif else parinfo[delta_rv_index].fixed = 1 ;Fix delta_rv on initial runs
            
            
  ;;; --- TELLURIC ---
            ;Name
            parinfo[h2o_depth_index].parname = 'h2o_depth'
            parinfo[co2ch4_depth_index].parname = 'co2ch4_depth'
            
            ;Scaling
                                ;;Unnecessary because they are
                                ;;necessarily between 0 and 1
            h2o_depth_parscale=1d0
            co2ch4_depth_parscale=1d0

            ;Values
            parinfo[h2o_depth_index].value = h2o_depth_guess
            parinfo[co2ch4_depth_index].value = co2ch4_depth_guess
            
            ;Step
            parinfo[h2o_depth_index].relstep = 1d-2
            parinfo[co2ch4_depth_index].relstep = 1d-2
            
            ;Constraints
            parinfo[h2o_depth_index].limited=[1,1]
            parinfo[co2ch4_depth_index].limited=[1,1]
            
            parinfo[h2o_depth_index].limits=[0d0,1d0]
            parinfo[co2ch4_depth_index].limits=[0d0,1d0]
            
            
  ;;; --- DELTA WL ---
            ;Name
            for pnum = 0, n_elements(delta_wl_index)-1 do begin
                parinfo[delta_wl_index[pnum]].parname = 'delta_wl_coeff_'+strtrim(pnum,2)
            endfor
            
            ;Scaling
            delta_wl_parscale=1d0/wl_start
    
    
            ;Values
            delta_wl_guess=wl_start*delta_wl_parscale
            parinfo[delta_wl_index].value = delta_wl_guess
            
            ;Step
            parinfo[delta_wl_index].step = 1d-8*delta_wl_parscale
            
            
            ;Constraints
            
            parinfo[delta_wl_index].mpside = 2
                                ;;Reasonable limit, if you want to use
                                ;;parinfo[delta_wl_index].limited = [1,1]
                                ;;parinfo[delta_wl_index].limits = [-1d-5, 1d-5]*abs(delta_wl_parscale)
            
            
  ;;; --- GH0 ---
            ;Name
            for pnum = 0, n_elements(gh0_coeff_index)-1 do begin
                if pnum eq 0 then parinfo[gh0_coeff_index[pnum]].parname = 'Sigma' else $
                parinfo[gh0_coeff_index[pnum]].parname = 'gh0_'+strtrim(pnum,2)
            endfor
            ;Scaling
            gh0_coeff_parscale=replicate(1d0, n_elements(gh0_coeff_index))

            ;Values
            for pnum = 0, n_elements(gh0_coeff_index)-1 do begin                
                if pnum eq 0 then parinfo[gh0_coeff_index[pnum]].value = sig_start else $
                  parinfo[gh0_coeff_index[pnum]].value = gh0_guess[pnum]
            endfor

            ;Step
            parinfo[gh0_coeff_index].step = replicate(1d-5, n_elements(gh0_coeff_index))
            
            ;Constraints
                                ;;Generous limits on sigma
            parinfo[gh0_coeff_index[0]].limited=[1,1]
            parinfo[gh0_coeff_index[0]].limits=[0.2d0,1.5d0]
                                ;;parinfo[gh0_coeff_index].fixed=1
                                ;;parinfo[gh0_coeff_index].mpside=2

            
  ;;; --- OTHER ---
            ;Name
            for pnum = 0, n_elements(other_index)-1 do begin
                parinfo[other_index[pnum]].parname = 'other_'+strtrim(pnum,2)
            endfor
            
            ;Scaling
            other_parscale=replicate(1d0, n_elements(other_index))
            
            ;Values
            parinfo[other_index].value = other_guess * other_parscale
            
            ;Step
            parinfo[other_index].relstep = 1d-3
            
            ;Constraints
            parinfo[other_index].limited=[1,1]
            parinfo[other_index].limits=[0,1.5]
    

  ;;; --- ALL TOGETHER NOW ---
            ;Scaling
            parscale_all=replicate(1d0, nparam_total)
            
            parscale_all[delta_rv_index] = delta_rv_parscale
            parscale_all[h2o_depth_index] = h2o_depth_parscale
            parscale_all[co2ch4_depth_index] = co2ch4_depth_parscale
            parscale_all[delta_wl_index] = delta_wl_parscale
            parscale_all[gh0_coeff_index] = gh0_coeff_parscale
            parscale_all[other_index] = other_parscale
            
    
            guess=parinfo.value
            scale=[rv_scale, h2o_depth_scale, co2ch4_depth_scale, delta_wl_scale, gh0_scale, other_scale]
            
            ;scale=scale*parscale_all
        
            
  ;;;  ----  Check Starting Values
            start_check=1
            start_check_print=0
            
            if start_check_print eq 1 then start_check=1
            
            if start_check eq 1 then begin
                
                for i=0,nparam_total-1 do begin
                    if start_check_print eq 1 then print, parinfo[i].parname
                    if start_check_print eq 1 then print, parinfo[i].value
                    if parinfo[i].limited[0] eq 1 then begin
                        if parinfo[i].value lt ((parinfo[i].limits)[0]) or $
                          parinfo[i].value gt ((parinfo[i].limits)[1]) then begin
                            if start_check_print ne 1 then print, parinfo[i].parname
                            print, parinfo[i].limits
                            print, "Illegal starting value"
                            print, "Changing limits for illegal case"
                            if parinfo[i].value lt ((parinfo[i].limits)[0]) then begin
                                diff=((parinfo[i].limits)[0])-parinfo[i].value
                                tmp=parinfo[i].limits
                                tmp[0]=tmp[0]-2*diff
                                parinfo[i].limits=tmp
                            endif else if parinfo[i].value gt ((parinfo[i].limits)[1]) then begin
                                diff=parinfo[i].value-((parinfo[i].limits)[1])
                                tmp=parinfo[i].limits
                                tmp[1]=tmp[1]+2*diff
                                parinfo[i].limits=tmp
                            endif
                            print, "limits changed, checking once more"
                            if parinfo[i].value lt ((parinfo[i].limits)[0]) or $
                              parinfo[i].value gt ((parinfo[i].limits)[1]) then begin
                                if start_check_print ne 1 then print, parinfo[i].parname
                                print, parinfo[i].limits
                                message, "Illegal starting value"
                            endif else print, "Limit successfully changed"
                            
                        endif
                    endif
                endfor
                
                if start_check_print eq 1 then begin
                    print, ' '
                    print, '-------------------------------'
                    print, 'Parscale: '
                    print, parscale_all
                    print, ' '
                    print, 'Scaled Guess: '
                    print, guess
                    print, ' '
                    print, 'Original Guess: '
                    print, guess/parscale_all
                    stop
                endif
            endif
            
            
            
            
            
;;;Run Minimization scheme
            
            
            if fmode eq 'amoeba' then begin
                r=amoeba3(ftol, scale=scale, p0=guess,function_name='rvmodel', $
                          function_value=fval, nmax=150000L)
                chi2=fval[0]
            endif else if fmode eq 'mpfit' then begin
                r=mpfit('rvmodel', bestnorm=chi2, ftol=ftol, parinfo=parinfo, status=status, /quiet)
                fmode='one_call'
                ;print, "R before onecall: ", r
                ocstr=rvmodel(r)
                lsf_temp=ocstr.lsf
                lsf=lsf_temp[*,0]
                wl_grid=ocstr.wl_grid
                model_select=ocstr.model
                obs_select=ocstr.obs
                err_select=ocstr.error
                chi2_nopen=total(((model_select-obs_select)/err_select)^2, /double)
                ;print, "R after oncecall: ", r
                fmode='mpfit'
            endif else if fmode eq 'tnmin' then begin
                ;r=tnmin('rvmodel', guess, autoderivative=1, bestmin=chi2,
                ;parinfo=parinfo, status=status, /quiet)
                r=tnmin('rvmodel',  autoderivative=1, bestmin=chi2, parinfo=parinfo, status=status, /quiet)
            endif else if fmode eq 'one_call' then begin
                fmode='mpfit'
                chi1_vec=rvmodel(guess)
                chi2=total(chi1_vec^2, /double)
                r=last_guess
                fmode='one_call'
                ocstr=rvmodel(guess)
                status=-111
                ;if visit eq 5 or visit eq 6 or visit eq 7 then stop
            endif else if fmode eq 'ga' then break
            
            
            
            
            
            ;stop clock
            process_time=toc(clock)
            
            
            ;OUTPUT RESULTS
            
            
            if n_elements(r) eq 1 then begin
                print, "Minimization scheme failed to converge"
                ;rescale
                r=last_guess/parscale_all
                chi2=-1
                
            endif else begin
                
                ;rescale
                ;print, "R as is:  ", r
                r=r/parscale_all
                ;print, "R scaled: ", r
                ;stop

                print, "---------------------------------"
                print, "Model ID: ", model_id_full
                print, "---------------------------------"
                
                print, chi2
                
                
            endelse
            
            print, "---------------------------------"
            print, "Processing took ", process_time, " seconds."
            
            
            if file_test(output_file) then begin
                fits_info, output_file, n_ext=prev_ext, /silent
                ;if prev_ext gt n_ABobj then begin
                ;    print, "Possibly too many extensions"
                ;    stop
                ;endif 
                extension=prev_ext+1
            endif else extension=1
            
            
            
            ;Flag bad runs
            if status eq 2 then status_flag=0 else status_flag=1
            if chi2 le chi2_nopen then penalty_flag=0 else penalty_flag=1
            if chi2 gt 100d0 and chi2 lt 10d4 then chi2_flag=0 else chi2_flag=1
            
            
            
            output_str={VISIT:visit+1L, $
                        GUESS:guess/parscale_all, $
                        RESULT:r, $
                        STATUS:status, $
                        DELTA_RV_INDEX:delta_rv_index, $
                        H2O_DEPTH_INDEX:h2o_depth_index, $
                        CO2CH4_DEPTH_INDEX:co2ch4_depth_index, $
                        DELTA_WL_INDEX:delta_wl_index, $
                        GH0_COEFF_INDEX:gh0_coeff_index, $
                        OTHER_INDEX:other_index, $
                        OVERSAMP:oversamp, $
                        FIRST_PIX:first_pix, $
                        NPIX_SELECT:npix_select, $
                        TIME:process_time, $
                        MODE:fmode, $
                        CHI2:chi2, $
                        CHI2_NOPEN:chi2_nopen, $
                        WL_GRID:wl_grid, $
                        MODEL_SELECT:model_select, $
                        OBS_SELECT:obs_select, $
                        ERR_SELECT:err_select, $
                        LSF:lsf, $
                        ;WL_START:wl_start, $
                        ;WL_RESULT:r[delta_wl_index], $
                        ;SIG_START:sig_start, $
                        ;SIG_RESULT:r[gh0_coeff_index[0]], $
                        MJD:mjd, $
                        BCV:bcv, $
                        DELTA_BCV:delta_bcv, $
                        DELTA_RV:r[0], $
                        FINAL_RV:r[0]+delta_bcv, $
                        STATUS_FLAG:status_flag, $
                        PENALTY_FLAG:penalty_flag, $
                        CHI2_FLAG:chi2_flag $
                       }
            
            

    if fmode ne 'one_call' then mwrfits, output_str, output_file

endfor

;for visit=1, n_ABobj-1 do begin
for visit=0, 0 do begin
    
    fmode='ga'
                                ;;READ IN MPFIT RUN RESULTS
    input_str=mrdfits(output_file, visit+1)
    guess=input_str.result
    
    if n_elements(n_lsf_ga) eq 0 then n_lsf_ga=4 ;MAKE KEYWORD
    gh0_coeff_index=input_str.gh0_coeff_index
    n_lsf_diff=n_lsf_ga-n_elements(gh0_coeff_index)
    n_other=n_elements(input_str.other_index)
    
    ga_guess=dblarr(n_elements(guess)+n_lsf_diff)
    ga_guess[0:gh0_coeff_index[0]]=guess[0:gh0_coeff_index[0]]
    if n_lsf_diff ge 1 then begin
        ga_guess[gh0_coeff_index[0]+1:gh0_coeff_index[0]+n_lsf_diff]=replicate(0.5d0, n_lsf_diff)
        ga_guess[gh0_coeff_index[0]+n_lsf_diff+1:*]=guess[gh0_coeff_index[0]+1:*]
    
;        help, guess
;        print, guess
;        help, ga_guess
;        print, ga_guess
;        stop
        
        
        ghga_coeff_index=lindgen(n_lsf_ga)+gh0_coeff_index[0]
        help, gh0_coeff_index
        print, gh0_coeff_index
        
        otherga_index=lindgen(n_other)+ ghga_coeff_index[-1]+1
        help, otherga_index
        print, otherga_index
        
    endif else begin
        ga_guess=guess
        ghga_coeff_index=gh0_coeff_index
        otherga_index=other_index
    endelse
    
    ndim=n_elements(ga_guess)
;    help, ndim
;    print, ndim



;save and rename old common block variables
parscale_all_save=parscale_all
parscale_all=replicate(1d0, ndim)
gh0_coeff_index_save=gh0_coeff_index
gh0_coeff_index=ghga_coeff_index
other_index_save=other_index
other_index=otherga_index

    ;;SCALE SO THAT THE PARAMETERS ARE have limits at 0 and 1
    

  ;;; --- DELTA_RV ---
            
            ;Scaling
            delta_rv_m=0.8d0
            delta_rv_b=ga_guess[delta_rv_index]-(delta_rv_m/2d0)
            
            
  ;;; --- TELLURIC ---
            
            ;Scaling

            h2o_m=0.01d0
            co2ch4_m=0.01d0
            h2o_b=ga_guess[h2o_depth_index]-(h2o_m/2d0)
            co2ch4_b=ga_guess[co2ch4_depth_index]-(co2ch4_m/2d0)
            
            
  ;;; --- DELTA WL ---
            
            ;Scaling
            delta_wl_m=2d-5
            delta_wl_b=ga_guess[delta_wl_index]-(delta_wl_m/2d0)
            
  ;;; --- GH0 ---

            ;Scaling
            sig_m=0.1d0
            sig_b=ga_guess[ghga_coeff_index[0]]-(sig_m/2d0)
            gh0_m=replicate(1d0, n_lsf_diff)
            gh0_b=replicate(-0.5d0, n_lsf_diff)
            
  ;;; --- OTHER ---
            
            ;Scaling
            nh3_m=0.01d0
            nh3_b=ga_guess[otherga_index[0]]-(nh3_m/2d0)
            other_m=0.02d0
            other_b=ga_guess[otherga_index[1:*]]-(other_m/2d0)

  ;;; --- ALL TOGETHER NOW ---

            param_m=[delta_rv_m, $
                     h2o_m, $
                     co2ch4_m, $
                     delta_wl_m, $
                     sig_m, $
                     gh0_m, $
                     nh3_m, $
                     other_m $
                     ]
            
            param_b=[delta_rv_b, $
                     h2o_b, $
                     co2ch4_b, $
                     delta_wl_b, $
                     sig_b, $
                     gh0_b, $
                     nh3_b, $
                     other_b $
                     ]
            
            lim=rebin([0,1], 2, ndim)
            fitness=0d0

            funa={param_m:param_m, param_b:param_b}
            save_gfit=dblarr(ndim)
            save_gfit=0d0
            save_igen=0
	    save_gbest=0d0
	    save_time=0d0

            ;sol=solber('ga_model',ndim, npop=npop, funa=funa, lim=lim, ngen_max=ngen_max, gfit_best=fitness, term_flag=0, term_fit=2000L, plot_flag=1, new_save_gen=save_gen, new_save_gfit=save_gfit, new_save_igen=save_igen, new_save_gbest=save_gbest, new_save_time=save_time)
            sol=solber('ga_model',ndim, funa=funa, lim=lim,gfit_best=fitness, term_flag=0, term_fit=2000L, plot_flag=1, new_save_gen=save_gen, new_save_gfit=save_gfit, new_save_igen=save_igen, new_save_gbest=save_gbest, new_save_time=save_time, _extra=ex)
            help, sol
            print, sol

            real_param=param_m*sol+param_b
            help, real_param
            print, real_param

            chi2=fitness

            help, chi2
            print, chi2
            print, input_str.chi2
            print, input_str.result

            npop=ex.npop
            ngen_max=ex.ngen_max
            
            output_str={raw_result:sol, $
                        scaled_result:real_param, $
                        chi2:chi2, $
			npop:npop, $
			ngen_max:ngen_max, $
			save_gen:save_gen, $
			save_gfit:save_gfit, $
			save_igen:save_igen, $
			save_gbest:save_gbest, $
			save_time:save_time $
                        }

            n_top=15L
            help, save_gen
            print, reform(save_gen[0,0:n_top-1], n_top)
            help, save_gfit
            print, save_gfit[0:n_top-1]
            help, save_igen
            print, save_igen[0:n_top-1]
            help, save_gbest
            help, save_time
            print, save_time[-1]

            fmode='one_call'
            octest=rvmodel(real_param)
            device, decomposed=1
            window, 1, xs=700, ys=500
            plot, octest.lsf[*,0], /xs
            window, 2, xs=700, ys=500
            plot, octest.residuals, ps=3, /xs
            window, 3, xs=700, ys=500
            plot, octest.obs, ps=6, /xs
            oplot, octest.model, ps=6, color=200
; fmode='ga'           
;stop
            
            ;if file_test(output_file) then begin
            ;    fits_info, output_file, n_ext=prev_ext, /silent
            ;    extension=prev_ext+1
            ;endif else extension=1
            output_file_new=strjoin([file_base, current_tag, 'ga'], '_') + '.fits'
            print, output_file
            mwrfits, output_str, output_file_new
            
        endfor
   stop     
!p.multi=0

end
