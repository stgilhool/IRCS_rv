
function templatefit, p, fakekey=fakekey

common funcargs, n_exp, $
  npix_model, $
  npix_select, $
  first_pix, $
  oversamp, $
  temp_display, $
  obs_select, $
  err_select, $
  temp_wl_over, $
  wl_template_shifted, $
  wl_grid_over, $
  tell_nh3_spectrum, $
  lsf_arr, $
  norm_factor, $
  mode, $
  template_over, $
  guess, $
  pscale, $
  last_guess, $
  chi2_nopen, $
  model_final

    clock=tic()

    last_guess=p

    ;Graph Template
    ;plot, temp_display, /xs, yr=[0.2,1.3]
;    wait,0.1

    npix_over=n_elements(template_over)

    model_over=dblarr(npix_model, n_exp)
    
    template_long=template_over
    
    template_long[first_pix*oversamp:(first_pix+npix_select)*oversamp-1]=p[0:-2]
    
    norm_scale=p[-1]

    ;Spline
    smoothspline=0
    if smoothspline eq 1 then begin

        template_spline=cspline(temp_wl_over, template_long, temp_wl_over)
        template_arr=rebin(template_spline, npix_over, n_exp)
        stellar_trial=interpol(template_arr, wl_template_shifted, wl_grid_over, /spline)

    endif else template_arr=rebin(template_long, npix_over, n_exp)

                                ;Stellar Template: RV shift the stellar template
    stellar_trial=interpol(template_arr, wl_template_shifted, wl_grid_over)
    
                                ;Now convolve to make full model
    product_spectrum=tell_nh3_spectrum * stellar_trial
    
        
            
;            product=product_spectrum[*,index]
    for exp=0, n_exp-1 do begin

        ;convclockt=ic()
        model_over[*,exp]= ircs_convolve(product_spectrum[*,exp], lsf_arr[*,*,exp], oversamp=oversamp, first_pix=first_pix, npix_select=npix_select)
        ;toc, convclock
    endfor

        
                                ;Downsample to IRCS resolution
    model = ircsrv_downsample(model_over, oversamp)

        
        
                                ;final IRCS sampled 2d array of models
                                ;for all exposures (n_exp x npix_select)
    model_final=model*norm_factor*norm_scale
    model_longvec=reform(transpose(model_final), n_elements(model_final))
    obs_longvec=reform(transpose(obs_select), n_elements(obs_select))
;    stop    


    ;Add roughness penalty
    roughpenalty=1
    if roughpenalty ge 1 then begin
                                ;Take first derivative of p
        xv=temp_wl_over[first_pix*oversamp:(first_pix+npix_select)*oversamp-1]
        yv=template_long[first_pix*oversamp:(first_pix+npix_select)*oversamp-1]
        dx=xv[1:*]-xv[0:-2]
        dy=yv[1:*]-yv[0:-2]
        der=dy/dx
        
        ;pscale=1d-5
        
        rpenalty=pscale*total(der[10*oversamp:*]^2, /double) ;FIX careful with this
        penstopon=0
        if penstopon eq 0 then penstop=[0] else penstop=randomu(seed,1)
        if penstop[0] gt .999 then stop
    endif

                ;Calculate chi2
    if mode eq 'amoeba' then begin
        chi2_allexp=((model_final-obs_select)/err_select)^2
        chi2_arr=total(chi2_allexp, 2, /double)
        chi2=total(chi2_arr, /double)
        chi2_nopen=chi2
        if roughpenalty ge 1 then chi2=chi2+rpenalty
        print, chi2_nopen, chi2
    endif
    if mode eq 'mpfit' or mode eq 'one_call' then begin
        chi_allexp=(obs_select-model_final)/err_select
        ;chi_arr=reform(transpose(chi_allexp), n_elements(chi_allexp))
        chi_arr=reform(chi_allexp, n_elements(chi_allexp))
        ;chi_arr=total(chi_allexp, 2, /double)
        chi2_nopen=total((chi_arr)^2, /double)
        if roughpenalty eq 1 then begin ; this is kinda wrong, but okay
            pen=rpenalty/n_elements(chi_arr)
            pen=sqrt(pen)
            poschi=where(chi_arr ge 0, poscount)
            negchi=where(chi_arr lt 0, negcount)
            if poscount gt 0 then chi_arr[poschi]=chi_arr[poschi]+pen
            if negcount gt 0 then chi_arr[negchi]=chi_arr[negchi]-pen
            chi2_tot=total(chi_arr^2, /double)
        endif else if roughpenalty eq 2 then begin

            negchi=where(chi_arr lt 0, negcount)
            chi2_tot=chi2_nopen+rpenalty
            chi2perpix=chi2_tot/n_elements(chi_arr)
            chiperpix=sqrt(chi2perpix)
            chi_arr=replicate(chiperpix, n_elements(chi_arr))
            if negcount gt 0 then chi_arr[negchi]=-1d0 * chi_arr[negchi]
;            if finite(chi2_nopen) eq 0 then stop
                
        endif
        
        print, chi2_nopen, total(chi_arr^2, /double)
        ;help, chi_arr
        ;help, p
        ;stop
    endif
    
    res_longvec=reform(transpose((obs_select-model_final)), n_elements(obs_select))
    ;i=1
    ;plot, obs_select[*,i], /xs, yr=[0,1.3]
    ;oplot, model_final[*,i], ps=3, color=200
    ;plot, obs_select[*,i]-model_final[*,i], ps=3, /xs
    ;plot, obs_longvec, /xs, yr=[0,1.3]
    ;oplot, model_longvec, ps=6, color=200
    ;plot, res_longvec, /xs, ps=6
    ;oplot, p, ps=3, color=200
    ;wait, 0.01
;    stop
    randomon=1
    if randomon eq 1 then begin 
        
        randomplot=randomu(seed, 1) 
    endif else randomplot=[0]
    if randomplot[0] gt 0.995 then begin
        
        !p.multi=[0,3,n_exp]
        for i=0, n_exp-1 do begin
            plot, obs_select[*,i], /xs, yr=[0,1.3]
            oplot, model_final[*,i], ps=3, color=200
            plot, obs_select[*,i]-model_final[*,i], ps=3, /xs, yr=[-0.1,0.1]
            plot, guess, /xs, yr=[0,1.3]
            oplot, template_arr[first_pix*oversamp:(first_pix+npix_select)*oversamp-1, i], color=200, ps=3
        endfor
    endif


    ;randomstop=randomu(seed, 1)
    randomstop=[0]
    if randomstop[0] gt 0.9 then begin
        for i=0, n_exp-1 do begin
            plot, obs_select[*,i], /xs, yr=[0,1.3]
            oplot, model_final[*,i], ps=3, color=200
            plot, obs_select[*,i]-model_final[*,i], ps=3, /xs
            print, "Obs: ", i, " Chi2: ", total((chi_allexp[*,i]^2), /double)
            stop
        endfor
    endif

    if mode eq 'amoeba' then return, chi2 
    if mode eq 'mpfit' then return, chi_arr
    if mode eq 'one_call' then begin
        ocstr={model_final:model_final, chi2_nopen:chi2_nopen, chi2_tot:chi2_tot}
;        stop
        return, ocstr
    endif


end



pro ircsrv_templatefit, epoch=epoch, object=object, trace=trace, visualize=visualize, first_pix=first_pix, npix_select=npix_select, mode=mode, input_file=input_file, s_template_numout=s_template_numout

common funcargs, n_exp, $
  npix_model, $
  npixselect, $
  firstpix, $
  oversamp, $
  temp_display, $
  obs_select, $
  err_select, $
  temp_wl_over, $
  wl_template_shifted, $
  wl_grid_over, $
  tell_nh3_spectrum, $
  lsf_arr, $
  norm_factor, $
  modemode, $
  template_over, $
  guess, $
  pscale, $
  last_guess, $
  chi2_nopen, $
  model_final


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
if n_elements(npix_select) eq 0 then npix_select=128L ;npix_select=npix-first_pix
if n_elements(mode) eq 0 then mode='mpfit'
modemode=mode


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
oversamp_old=model_par.oversamp
oversamp=7L
npix_lsf=(oversamp*10L)+1L
bigc=299792458.D
npix_model=npix_select*oversamp
npix_over=npix*oversamp
npix_model_old=npix_select*oversamp_old
npix_over_old=npix*oversamp_old


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
s_template=0 ;Tuned var
if s_template eq 0 then begin
    templatefile=templatepath+object+'_template.dat'
    readcol, templatefile, temp_wl, temp_spec, format='D,D'
    if n_elements(temp_wl) ne 7L*npix then message, "You were wrong about the old template sampling"
    if n_elements(temp_spec) ne 7L*npix then message, "You were wrong about the old template sampling"
    temp_oversamp=7L
    norm_start=(10*temp_oversamp) > (first_pix*temp_oversamp)
    norm_end= (first_pix*temp_oversamp + npix_select*temp_oversamp - 1) < ((npix-10)*temp_oversamp-1)
    tempnorm=(max(temp_spec[norm_start:norm_end]))
    temp_spec=temp_spec/tempnorm
    
endif else begin
    templatefile='/home/stgilhool/RV_projects/IRCS_rv/data/smooth_penalty_test/test19/penaltytest.fits'  
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
print, tempnorm
wait, 5


;;;-----Finished reading in template----






;OLD READIN TEMPLATE
; templatefile=templatepath+object+'_template.dat'
; readcol, templatefile, temp_wl, temp_spec, format='D,D'


; ;oversample the template and template_wl
; if n_elements(temp_wl) eq npix then begin
;     temp_wl_over=interpol(temp_wl, x, xx)
;     template_over=interpol(temp_spec, temp_wl, temp_wl_over)
; endif else if oversamp eq 7 then begin
;     temp_wl_over=temp_wl
;     template_over=temp_spec
;     template_over_select=template_over[first_pix*oversamp:(first_pix+npix_select)*oversamp-1]
; endif else if oversamp ne 7 then begin
;     if n_elements(temp_wl) ne 7*npix then message, 'you are wrong about template npix'
;     xx_old=(dindgen(npix*oversamp_old)-(oversamp_old/2))/oversamp_old
;     temp_wl_over=interpol(temp_wl, xx_old, xx)
;     template_over=interpol(temp_spec, temp_wl, temp_wl_over)
;     template_over_select=template_over[first_pix*oversamp:(first_pix+npix_select)*oversamp-1]
; endif



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


;if visit eq 5 or visit eq 7 then continue

;varray=[0,1,2,3,4,6,8,9,10,11,12,13]
if first_pix ne 349L then varray=[0,1,2,3,4,5,6,8,9,10,12,13] $
  else varray=[0,1,3,4,5,6,7,8,9,12]
;varray=[0,1,2,3]


;Define some arrays
n_exp=n_elements(varray)
if n_elements(input_file) eq 0 then rfile='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_349_379_Mar16bigpen_r1.fits' else rfile=input_file
f=mrdfits(rfile, 1)
re=f.result

params=dblarr(n_elements(re), n_exp)
lsf=dblarr(npix_lsf, npix_model, n_exp)
tell_nh3_spectrum=dblarr(npix_over, n_exp)

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



foreach visit, varray, index do begin

;Make observation array
    obs[*, index]=ABspec_arr[visit, *]
;get error for this observation
    err[*, index]=ABerr_arr[visit, *]
;Mask entries
;Cut off pixels on either end of full 1024 pixel spectrum
    npix_trim_start=10L
    npix_trim_end=10L
    bigerr=1d10                 ;error value for masked pixels
    if n_elements(err[*,index]) eq 1024L then begin
        err[0:npix_trim_start-1, index]=bigerr
        err[-1L*npix_trim_end:-1, index]=bigerr
    endif



;READ IN PARAMETERS FROM 2-basis 2-run FITS
    ;rfile='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_127_run2sign.fits'

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
    params[index]=re


;;;BEGIN MODELING

;Perturb wl soln
    delta_wl=poly(xx, delta_wl_coeff)
    wl_grid_over[*, index]=wl_soln_over+delta_wl

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
    wl_template_shifted[*, index]=temp_wl_over*(1d0 + 1000d0*(delta_rv[index])/bigc)


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
        if npix_select ge 128L then norm_factor[*,index]=interpol(norm_pts_y, norm_pts_x, norm_pts_xx, /spline) $
          else norm_factor[*,index]=interpol(norm_pts_y, norm_pts_x, norm_pts_xx)
        
    endif

endforeach

                                ;Truncate observations and errors to model length
                                ;Indices that we are using
select_x = lindgen(npix)
select_index = select_x[first_pix:first_pix+npix_select-1]
select_xx= lindgen(npix_over)
select_index_over = select_xx[first_pix*oversamp:(first_pix*oversamp)+npix_model-1]

temp_display=template_over[select_index_over]
;first_template=rebin(template_over_select, npix_over, n_exp)





;NOW DO THE BILLIONS OF TRIALS

n_trials=1e2
;maxiter=1e9 - 1


                                ;Make 3D arrays out of some of the
                                ;arrays
; wl_template_shifted_arr=rebin(wl_template_shifted, npix_over, n_exp, n_trials)
; wl_grid_over_arr=rebin(wl_grid_over, npix_over, n_exp, n_trials)
; tell_nh3_spectrum_arr=rebin(tell_nh3_spectrum, npix_over, n_exp, n_trials)
; norm_factor_arr=rebin(norm_factor, npix_select, n_exp, n_trials)    
obs_select=rebin(obs[select_index, *], npix_select, n_exp)
err_select=rebin(err[select_index, *], npix_select, n_exp)
;lsf_arr=rebin(lsf, npix_lsf, npix_model, n_exp, n_trials)
lsf_arr=temporary(lsf)

;stop

;FIX this part about normalization


guess=template_over_select
;tempnorm=max(template_over_select[oversamp*10:*]) ;careful here!
;print, tempnorm
;tempnorm_new= 1.0089513d0
;print, tempnorm
;stop
;guess=guess/tempnorm_new

;tempnorm=max(guess[oversamp*10:*]) ;careful here!

;guess=guess/tempnorm

;Shift template values by small random amts across template
guess_orig=guess ;backup
guess=(1d0+randomn(seed, n_elements(guess))*0.005d0)*guess
;guess[150:180]=replicate(0.6d0, 31)+randomn(seed, 31)*0.01


guess=[guess,1d0]

toohigh=where(guess gt 1, highcount)
toolow=where(guess lt 0, lowcount)
if highcount gt 0 then guess[toohigh]=1d0
if lowcount gt 0 then guess[toolow]=0d0
ftol=1d-10



norm_factor=norm_factor;*tempnorm
plot, guess
;stop
!p.multi=[0,1,2]


maxiter=0L
for scaleiter=0, maxiter do begin
    
;scalelist=[5d-5, 2.5d-5, 1d-5, 7.5d-6, 5d-6]  
;pscale=scalelist[scaleiter]
    pscale=(scaleiter+1)*3d-6
;pscale=(scaleiter)*1d-4+(1d-3 - 1d-4/(maxiter/2.)) 
;pscale=1d-1^(scaleiter+1)
;pscale=1d-1^(scaleiter+5)
;pscale=0d0
    if modemode eq 'mpfit' then begin
        
        parinfo = replicate({fixed:0, limited:[1,1], $
                             limits:[0.D,1.01D0]}, n_elements(guess))
        
        if first_pix lt 10 then parinfo[0:(10-first_pix)*oversamp].fixed=1
        
        r=mpfit('templatefit',guess, bestnorm=chi2, ftol=ftol, $
                parinfo=parinfo, status=status, nfev=ncalls, niter=niter, /quiet)
        print, status
        
        final_template=template_over
        final_template[first_pix*oversamp:(first_pix+npix_select)*oversamp-1]=r[0:-2]
        
        modemode='one_call'
        oc=templatefit(r)
        
        model_final=oc.model_final
        chi2_nopen=oc.chi2_nopen
        chi2_tot=oc.chi2_tot
        modemode='mpfit'
    endif
    
    if modemode eq 'amoeba' then begin
        window, 1, xsize=1200, ysize=650
        scale=replicate(1d-2, n_elements(template_over_select))
        
        r=amoeba3(ftol, scale=scale, p0=guess,function_name='templatefit', $
                  function_value=fval, ncalls=ncalls, nmax=150000L)
        niter=ncalls
        status=-888
        if n_elements(r) eq 1 then begin
            r=last_guess
            status=-999
        endif
        chi2=fval[0]
        final_template=template_over
        final_template[first_pix*oversamp:(first_pix+npix_select)*oversamp-1]=r[0:-2]
        modemode='one_call'
        oc=templatefit(r)
        chi2_tot=oc.chi2_tot
        modemode='amoeba'
    endif
    
    ;Add rms to calculation
    res=model_final-obs_select
    rmean=rebin(reform(mean(res[10:*,*], dimension=1), 1, n_exp), npix_select, n_exp)
    var=(res-rmean)^2
    rms=sqrt(mean(var[10:*,*], dimension=1))
    
    
    penstr={template:r[0:-2], $
            guess:guess, $
            params:params, $
            chi2:chi2, $
            chi2_tot:chi2_tot, $
            template_spec:final_template, $
            template_wl:temp_wl_over, $
            pscale:pscale, $
            status:status, $
            first_pix:first_pix, $
            npix_select:npix_select, $
            model_arr:model_final, $
            obs_arr:obs_select, $
            err_arr:err_select, $
            chi2_nopen:chi2_nopen,$
            rms:rms, $
            ncalls:ncalls, $
            niter:niter $
           }
    
;if scaleiter eq 0 then stop
    
    if n_elements(s_template_numout) eq 0 then s_template_numout=23
    outfil='/home/stgilhool/RV_projects/IRCS_rv/data/smooth_penalty_test/test'+strtrim(s_template_numout,2)+'/penaltytest.fits'
    
    mwrfits, penstr, outfil
    
    
    
endfor

compare_template=0

if compare_template eq 1 then begin
;Do one final call with the original template
    modemode='one_call'
    oc=templatefit(guess_orig)
    model_final=oc.model_final
    chi2_tot=oc.chi2_tot
    chi2_nopen=oc.chi2_nopen
    chi2=chi2_tot
    status=-111
    ncalls=1
    niter=1
    
  
    ;Add rms to calculation
    res=model_final-obs_select
    rmean=rebin(reform(mean(res[10:*,*], dimension=1), 1, n_exp), npix_select, n_exp)
    var=(res-rmean)^2
    rms=sqrt(mean(var[10:*,*], dimension=1))
    
    
    penstr={template:r, $
            guess:guess, $
            params:params, $
            chi2:chi2, $
            chi2_tot:chi2_tot, $
            template_spec:final_template, $
            template_wl:temp_wl_over, $
            pscale:pscale, $
            status:status, $
            first_pix:first_pix, $
            npix_select:npix_select, $
            model_arr:model_final, $
            obs_arr:obs_select, $
            err_arr:err_select, $
            chi2_nopen:chi2_nopen,$
            rms:rms, $
            ncalls:ncalls, $
            niter:niter $
           }

;if scaleiter eq 0 then stop
    
    mwrfits, penstr, outfil
endif

    
;stop
        
        

!p.multi=0



;bestfile= '/home/stgilhool/RV_projects/IRCS_rv/data/GJ273/templatefit.fits'
;outstr={template:final_template, $
;        chi2:chi2 $
;        }

;mwrfits, outstr, bestfile

;stop


end
