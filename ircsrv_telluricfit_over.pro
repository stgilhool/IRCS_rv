function telluricfit, p, fakekey=fakekey

common funcargs, n_exp, $
  npix_model, $
  npix_select, $
  npix_over, $
  first_pix, $
  oversamp, $
  obs_select, $
  err_select, $
  obs_select_over, $
  err_select_over, $
  wl_grid_over, $
  temp_wl_over, $
  template_over, $
  stellar, $
  nh3_spectrum, $
  h2o_over, $
  co2ch4_over, $
  wl_soln_tell, $
  wl_soln_tell_select, $
  select_tell_index, $
  lsf_arr, $
  norm_factor, $
  norm_factor_over, $
  mode, $
  guess, $
  pscale, $
  last_guess, $
  chi2_nopen, $
  model_final

    clock=tic()

    last_guess=p
    bigc=299792458.D

    model_over=dblarr(npix_model, n_exp)
    ;tell_penalty=dblarr(n_elements(select_tell_index), n_exp)
    tell_penalty_h2o=dblarr(n_elements(select_tell_index), n_exp)
    tell_penalty_co2ch4=dblarr(n_elements(select_tell_index), n_exp)
    ;stellar=dblarr(npix_over, n_exp)
    ;telluric_model=dblarr(n_elements(
    
    ;modify telluric with free parameters
    h2o_over[select_tell_index]=p[0:n_elements(wl_soln_tell_select)-1]
    co2ch4_over[select_tell_index]=p[n_elements(wl_soln_tell_select):2*n_elements(wl_soln_tell_select)-1]
    h2o_depth=p[2*n_elements(wl_soln_tell_select):2*n_elements(wl_soln_tell_select)+n_exp-1]
    co2ch4_depth=p[2*n_elements(wl_soln_tell_select)+n_exp:2*n_elements(wl_soln_tell_select)+2*n_exp-1]
    delta_rv=p[2*n_elements(wl_soln_tell_select)+2*n_exp:2*n_elements(wl_soln_tell_select)+3*n_exp-1]



    

; help, p
; help, h2o_over
; help, co2ch4_over
; help, h2o_depth
; help, co2ch4_depth
; stop
    
;Loop through exposures for multiply and convolve
    for exp=0, n_exp-1 do begin


            ;Shift stellar wl
        wl_template_shifted=temp_wl_over*(1d0 + delta_rv[exp]/bigc)
        stellar=interpol(template_over, wl_template_shifted, wl_grid_over[*, exp])
        
        ;scale optical depth of telluric
        h2o_scaled=h2o_over^h2o_depth[exp]
        co2ch4_scaled=co2ch4_over^co2ch4_depth[exp]
        
                                ;save telluric "template" for roughness penalty
;        tell_penalty[*,exp]=h2o_scaled[select_tell_index]*co2ch4_scaled[select_tell_index]
        tell_penalty_h2o[*,exp]=h2o_scaled[select_tell_index]
        tell_penalty_co2ch4[*,exp]=co2ch4_scaled[select_tell_index]

        ;put telluric onto wavelength grid
        h2o=interpol(h2o_scaled, wl_soln_tell, wl_grid_over[*,exp])
        co2ch4=interpol(co2ch4_scaled, wl_soln_tell, wl_grid_over[*,exp])

        ;telluric_model[*,exp]=interpol(h2o*co2ch4, wl_grid_over[*,exp], wl_grid_select[*,exp])
        

        ;Product of stell, tell and nh3
        product_spectrum=stellar*nh3_spectrum[*,exp]*h2o*co2ch4

        model_over[*,exp]= ircs_convolve(product_spectrum, lsf_arr[*,*,exp], oversamp=oversamp, first_pix=first_pix, npix_select=npix_select)

    endfor

        
                                ;Downsample to IRCS resolution
;    model = ircsrv_downsample(model_over, oversamp)
        
        
                                ;final IRCS sampled 2d array of models
                                ;for all exposures (n_exp x npix_select)
;    model_final=model*norm_factor
    model_final=model_over*norm_factor_over
;    model_longvec=reform(transpose(model_final), n_elements(model_final))
;    obs_longvec=reform(transpose(obs_select), n_elements(obs_select))
;    stop    


    ;Add roughness penalty
    roughpenalty=1
    if roughpenalty ge 1 then begin
                                ;Take first derivative of p
        xv=wl_soln_tell_select
        ;yv=tell_penalty[*,0]
        yv_h2o=tell_penalty_h2o[*,0]
        yv_co2ch4=tell_penalty_co2ch4[*,0]
        dx=xv[1:*]-xv[0:-2]
        ;dy=yv[1:*]-yv[0:-2]
        dy_h2o=yv_h2o[1:*]-yv_h2o[0:-2]
        dy_co2ch4=yv_co2ch4[1:*]-yv_co2ch4[0:-2]
        ;der=dy/dx
        der_h2o=dy_h2o/dx
        der_co2ch4=dy_co2ch4/dx
        
        ;pscale=1d-5
        
;        rpenalty=pscale*total(der[10*oversamp:*]^2, /double) ;FIX
;        careful with this
        rpenalty_h2o=pscale*total(der_h2o[10*oversamp:*]^2, /double) ;FIX careful with this
        rpenalty_co2ch4=pscale*total(der_co2ch4[10*oversamp:*]^2, /double) ;FIX careful with this
        rpenalty=rpenalty_h2o+rpenalty_co2ch4
        penstopon=0
        if penstopon eq 0 then penstop=[0] else penstop=randomu(seed,1)
        if penstop[0] gt .999 then stop
    endif

                ;Calculate chi2
    if mode eq 'amoeba' then begin
        chi2_allexp=((model_final-obs_select_over)/err_select_over)^2
        chi2_arr=total(chi2_allexp, 2, /double)
        chi2=total(chi2_arr, /double)
        chi2_nopen=chi2
        if roughpenalty ge 1 then chi2=chi2+rpenalty
        print, chi2_nopen, chi2
    endif
    if mode eq 'mpfit' or mode eq 'one_call' then begin
        chi_allexp=(obs_select_over-model_final)/err_select_over
;         help, chi_allexp
;         help, obs_select_over
;         help, model_final
;         help, err_select_over
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

; help, chi2_tot
; help, chi_arr
; stop    

;    res_longvec=reform(transpose((obs_select_over-model_final)), n_elements(obs_select_over))
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
    if randomplot[0] gt 0.99 then begin
        
        !p.multi=[0,3,n_exp]
        for i=0, n_exp-1 do begin
            plot, obs_select_over[*,i], /xs, yr=[0,1.3]
            oplot, model_final[*,i], ps=3, color=200
            plot, obs_select_over[*,i]-model_final[*,i], ps=3, /xs, yr=[-0.1,0.1]
            plot, tell_penalty_h2o[*,i], /xs, yr=[0,1.3]
            oplot, tell_penalty_co2ch4[*,i], color=200
        endfor
    endif


    ;randomstop=randomu(seed, 1)
    randomstop=[0]
    if randomstop[0] gt 0.9 then begin
        for i=0, n_exp-1 do begin
            plot, obs_select_over[*,i], /xs, yr=[0,1.3]
            oplot, model_final[*,i], ps=3, color=200
            plot, obs_select_over[*,i]-model_final[*,i], ps=3, /xs
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



pro ircsrv_telluricfit_over, epoch=epoch, object=object, trace=trace, visualize=visualize, first_pix=first_pix, npix_select=npix_select, mode=mode


common funcargs, n_exp, $
  npix_model, $
  npixselect, $
  npix_over, $
  firstpix, $
  oversamp, $
  obs_select, $
  err_select, $
  obs_select_over, $
  err_select_over, $
  wl_grid_over, $
  temp_wl_over, $
  template_over, $
  stellar, $
  nh3_spectrum, $
  h2o_over, $
  co2ch4_over, $
  wl_soln_tell, $
  wl_soln_tell_select, $
  select_tell_index, $
  lsf_arr, $
  norm_factor, $
  norm_factor_over, $
  modemode, $
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
oversamp=3L
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
templatefile='/home/stgilhool/RV_projects/IRCS_rv/data/smooth_penalty_test/test16/penaltytest.fits'
temp_ext=30
templatestr=mrdfits(templatefile, temp_ext)
temp_spec=templatestr.template_spec
temp_wl=templatestr.template_wl

if oversamp ne 3 then begin
    temp_oversamp=3L
    xx_temp=(dindgen(npix*temp_oversamp)-(temp_oversamp/2))/temp_oversamp
    temp_wl_over=interpol(temp_wl, xx_temp, xx)
    template_over=interpol(temp_spec, temp_wl, temp_wl_over)

endif else begin
    template_over=temp_spec
    temp_wl_over=temp_wl
endelse

template_over_select=template_over[first_pix*oversamp:(first_pix+npix_select)*oversamp-1]



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
varray=[0,1,2,3,4,6,10,11,12,13]
;varray=[0,1,2,3]


;Define some arrays
n_exp=n_elements(varray)
lsf=dblarr(npix_lsf, npix_model, n_exp)
nh3_spectrum=dblarr(npix_over, n_exp)
;h2o_spectrum=dblarr(npix_over, n_exp)
;co2ch4_spectrum=dblarr(npix_over, n_exp)

wl_grid_over=dblarr(npix_over, n_exp)
err=dblarr(npix, n_exp)
final_spectrum=dblarr(npix_model)
delta_rv=dblarr(n_exp)
wl_grid_over=dblarr(npix_over, n_exp)
wl_grid_down=dblarr(npix, n_exp)
wl_grid_over_select=dblarr(npix_model, n_exp)
wl_grid_down_select=dblarr(npix_select, n_exp)
wl_template_shifted=dblarr(npix_over, n_exp)
norm_factor=dblarr(npix_select, n_exp)
norm_factor_over=dblarr(npix_model, n_exp)
obs=dblarr(npix, n_exp)
h2o_depth=dblarr(n_exp)
co2ch4_depth=dblarr(n_exp)


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
;    rfile='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_127_run2sign.fits'
    rfile='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_127_smooth30_r2.fits'

    f=mrdfits(rfile, visit+1)
    re=f.result

    delta_rv[index]=re[f.delta_rv_index]
    delta_wl_coeff=re[f.delta_wl_index]
    h2o_depth[index]=re[f.h2o_depth_index]
    co2ch4_depth[index]=re[f.co2ch4_depth_index]
    delta_wl_coeff=re[f.delta_wl_index]
    gh0_coeff=re[f.gh0_coeff_index]
    gh1_coeff=re[f.gh1_coeff_index]
    other=re[f.other_index]

    tau_scale=other[0]	


;;;BEGIN MODELING

;Perturb wl soln
    delta_wl=poly(xx, delta_wl_coeff)
    wl_grid_over[*, index]=wl_soln_over+delta_wl

    delta_wl_down=poly(x, delta_wl_coeff)
    wl_grid_down[*, index]=wl_soln+delta_wl_down

    wl_grid_down_select[*, index]=wl_grid_down[first_pix:first_pix+npix_select-1, index]
    wl_grid_over_select[*, index]=wl_grid_over[first_pix*oversamp:(first_pix+npix_select)*oversamp-1 ,index]

;NH3 Lab Spectrum: Interpolate onto the oversampled trial grid
;adjust optical depth first
    int_lab_depth=int_lab^tau_scale
    nh3_spectrum[*,index]=interpol(int_lab_depth, wl_lab, wl_grid_over[*, index])
    



;Make LSF

    lsf[*,*,index] = ircsrv_lsf(gh0_coeff, gh1_vec=gh1_coeff, oversamp=oversamp, first_pix=first_pix, npix_select=npix_select, neg_penalty=penalty)


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
        norm_factor_over[*,index]=interpol(norm_factor[*,index], x[0:npix_select-1], xx[0:npix_model-1])
        
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






obs_select=rebin(obs[select_index, *], npix_select, n_exp)
err_select=rebin(err[select_index, *], npix_select, n_exp)


obs_select_over=dblarr(npix_model, n_exp)
err_select_over=dblarr(npix_model, n_exp)
;obs_select_over=interpol(obs_select, wl_grid_down_select, wl_grid_over_select)
for exp=0, n_exp-1 do begin
    obs_select_over[*,exp]=interpol(obs_select[*,exp], wl_grid_down_select[*,exp], wl_grid_over_select[*,exp], /spline)
    err_select_over[*,exp]=reform(rebin(reform(err_select[*,exp], 1, npix_select), oversamp, npix_select), npix_model)
endfor

; help, obs_select
; help, obs_select_over
; help, err_select
; help, err_select_over
; help, wl_grid_down
; help, wl_grid_down_select
; help, wl_grid_over
; help, wl_grid_over_select
; stop


;lsf_arr=rebin(lsf, npix_lsf, npix_model, n_exp, n_trials)
lsf_arr=temporary(lsf)



!p.multi=[0,1,2]
extra_pix=4L
xx_tell=(dindgen((npix+extra_pix)*oversamp)-(oversamp/2))/oversamp-(extra_pix/2)
xx_tell_select=(dindgen((npix_select+extra_pix)*oversamp)-(oversamp/2))/oversamp-(extra_pix/2 + first_pix)

select_tell_xx=lindgen(n_elements(xx_tell))
select_tell_index=select_tell_xx[first_pix*oversamp:(first_pix*oversamp)+(npix_select+extra_pix)*oversamp - 1]

wl_soln_tell=poly(xx_tell, wl_coeff)
wl_soln_tell_select=poly(xx_tell_select, wl_coeff)

h2o_over=interpol(h2o, wl_telluric, wl_soln_tell)
co2ch4_over=interpol(co2ch4, wl_telluric, wl_soln_tell)


h2o_guess=interpol(h2o, wl_telluric, wl_soln_tell_select)
co2ch4_guess=interpol(co2ch4, wl_telluric, wl_soln_tell_select)

guess_orig=[h2o_guess, co2ch4_guess, h2o_depth, co2ch4_depth, delta_rv] ;backup
;Shift guess values a bit so that mpfit works
guess=guess_orig
guess[0:n_elements(wl_soln_tell_select)*2-1]=(1d0+randomn(seed, n_elements(guess_orig)-2*n_exp)*0.001d0)*guess_orig[0:n_elements(wl_soln_tell_select)*2 -1]
;guess[150:180]=replicate(0.6d0, 31)+randomn(seed, 31)*0.01


; help, wl_soln_tell
; help, wl_soln_tell_select
; help, guess_orig
; help, guess
; help, h2o_guess
; help, co2ch4_guess
; help, h2o_depth
; help, co2ch4_depth
; stop

maxiter=3L
for scaleiter=0, maxiter do begin

        ;scalelist=[5d-5, 2.5d-5, 1d-5, 7.5d-6, 5d-6]  
        ;pscale=scalelist[scaleiter]
    ;pscale=(scaleiter+1)*1d-5
	;pscale=(scaleiter)*1d-4+(1d-3 - 1d-4/(maxiter/2.)) 
	;pscale=1d-1^(scaleiter+1)
	pscale=1d-1^(scaleiter+3)
	;pscale=0d0
    if modemode eq 'mpfit' then begin

        parinfo = replicate({fixed:0, limited:[1,1], $
                             limits:[0.D,1.001D0]}, n_elements(guess))
        
        parinfo[0:10*oversamp].fixed=1
        parinfo[n_elements(wl_soln_tell_select):n_elements(wl_soln_tell_select)+(10*oversamp)].fixed=1
        parinfo[2*n_elements(wl_soln_tell_select)+2*n_exp:2*n_elements(wl_soln_tell_select)+3*n_exp-1].limits=[-2000D, 2000D]
        parinfo[2*n_elements(wl_soln_tell_select)+2*n_exp:2*n_elements(wl_soln_tell_select)+3*n_exp-1].fixed=1

        r=mpfit('telluricfit',guess, bestnorm=chi2, ftol=ftol, $
                parinfo=parinfo, status=status, nfev=ncalls, niter=niter, /quiet)
        print, "Status: ", status
        
        
        h2o_final=h2o_over
        h2o_final[select_tell_index]=r[0:n_elements(wl_soln_tell_select)-1]

        co2ch4_final=co2ch4_over
        co2ch4_final[select_tell_index]=r[n_elements(wl_soln_tell_select):2*n_elements(wl_soln_tell_select)-1]
        
        h2o_depth_final=r[2*n_elements(wl_soln_tell_select):2*n_elements(wl_soln_tell_select)+n_exp-1]
        co2ch4_depth_final=r[2*n_elements(wl_soln_tell_select)+n_exp:2*n_elements(wl_soln_tell_select)+2*n_exp-1]
        ;delta_rv_final=r[2*n_elements(wl_soln_tell_select)+2*n_exp:2*n_elements(wl_soln_tell_select)+3*n_exp-1]

        modemode='one_call'
        oc=telluricfit(r)
        
        model_final=oc.model_final
        chi2_nopen=oc.chi2_nopen
        chi2_tot=oc.chi2_tot
        modemode='mpfit'
    endif
    
    if modemode eq 'amoeba' then begin
        window, 1, xsize=1200, ysize=650
        scale=replicate(1d-2, n_elements(template_over_select))
        
        r=amoeba3(ftol, scale=scale, p0=guess,function_name='telluricfit', $
                  function_value=fval, ncalls=ncalls, nmax=150000L)
        niter=ncalls
        status=-888
        if n_elements(r) eq 1 then begin
            r=last_guess
            status=-999
        endif
        chi2=fval[0]

        h2o_final=h2o_over
        h2o_final[select_tell_index]=r[0:n_elements(wl_soln_tell_select)-1]

        co2ch4_final=co2ch4_over
        co2ch4_final[select_tell_index]=r[n_elements(wl_soln_tell_select):2*n_elements(wl_soln_tell_select)-1]
        
        h2o_depth_final=r[2*n_elements(wl_soln_tell_select):2*n_elements(wl_soln_tell_select)+n_exp-1]
        co2ch4_depth_final=r[2*n_elements(wl_soln_tell_select)+n_exp:2*n_elements(wl_soln_tell_select)+2*n_exp-1]
        ;delta_rv_final=r[2*n_elements(wl_soln_tell_select)+2*n_exp:2*n_elements(wl_soln_tell_select)+3*n_exp-1]

        modemode='one_call'
        oc=telluricfit(r)
        chi2_tot=oc.chi2_tot
        modemode='amoeba'
    endif
    
    
    penstr={h2o:h2o_final, $
            co2ch4:co2ch4_final, $
            delta_rv_guess:delta_rv, $
            delta_rv_final:delta_rv, $
            h2o_depth_guess:h2o_depth, $
            h2o_depth_final:h2o_depth_final, $
            co2ch4_depth_guess:co2ch4_depth, $
            co2ch4_depth_final:co2ch4_depth_final, $
            wl_telluric:wl_soln_tell, $
            wl_telluric_select:wl_soln_tell_select, $
            wl_grid:wl_grid_over, $
            guess:guess, $
            chi2:chi2, $
            chi2_tot:chi2_tot, $
            pscale:pscale, $
            status:status, $
            first_pix:first_pix, $
            npix_select:npix_select, $
            model_arr:model_final, $
            obs_arr:obs_select_over, $
            err_arr:err_select_over, $
            chi2_nopen:chi2_nopen,$
            ncalls:ncalls, $
            niter:niter $
           }
    
if scaleiter eq 0 then stop
    
    outfil='/home/stgilhool/RV_projects/IRCS_rv/data/telluric_test/test4/tellurictest.fits'
    
    mwrfits, penstr, outfil
    
    
    
endfor

compare_template=1

if compare_template eq 1 then begin
;Do one final call with the original template
    modemode='one_call'
    oc=telluricfit(guess_orig)
    model_final=oc.model_final
    chi2_tot=oc.chi2_tot
    chi2_nopen=oc.chi2_nopen
    chi2=chi2_tot
    status=-111
    ncalls=1
    niter=1
    
 
    
    penstr={h2o:h2o_final, $
            co2ch4:co2ch4_final, $
            delta_rv_guess:delta_rv, $
            delta_rv_final:delta_rv, $
            h2o_depth_guess:h2o_depth, $
            h2o_depth_final:h2o_depth_final, $
            co2ch4_depth_guess:co2ch4_depth, $
            co2ch4_depth_final:co2ch4_depth_final, $
            wl_telluric:wl_soln_tell, $
            wl_telluric_select:wl_soln_tell_select, $
            wl_grid:wl_grid_over, $
            guess:guess, $
            chi2:chi2, $
            chi2_tot:chi2_tot, $
            pscale:pscale, $
            status:status, $
            first_pix:first_pix, $
            npix_select:npix_select, $
            model_arr:model_final, $
            obs_arr:obs_select_over, $
            err_arr:err_select_over, $
            chi2_nopen:chi2_nopen,$
            ncalls:ncalls, $
            niter:niter $
           }

    
    mwrfits, penstr, outfil
endif

    
stop
        
        

!p.multi=0



end
