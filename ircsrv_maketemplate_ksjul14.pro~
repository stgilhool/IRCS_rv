function rvfit, p

common rvinfo, wl_soln, wl_soln_full, samp_index, wl_soln_over_full, int_lab_over_full, int_obs, lsf, npix_lsf, npix_model, first_pix, npix_select, oversamp, nnorm, err, visualize, npix_trim_start, npix_trim_end, min_type, nparam_other, last_guess, nlines, COwl_shift, h2o, co2ch4, tellwl, visit, iter, template_residuals,res, stellar_template, stellar_template_over, finalplot, int_lab, wl_lab, lab_depth_switch, norm_switch

lab_depth=lab_depth_switch
norm=norm_switch

last_guess=p

delta_rv = 0d0

co2ch4_tau=p[0]
h2o_tau=p[1]
if lab_depth eq 1 then nh3_tau=p[2]
if norm eq 1 then norm_pts_y=p[-1*nnorm:-1]


bigc=299792458d0 ;m/s 



;;;;;;;;;;;;;;;
wl_shifted=wl_soln_over_full*(1d0+(delta_rv/bigc))




;stop

;Deal with telluric
;scale telluric depth
h2o_depth=h2o^h2o_tau
co2ch4_depth=co2ch4^co2ch4_tau
if lab_depth eq 1 then nh3_depth=int_lab^nh3_tau

h2o_shift=interpol(h2o_depth, tellwl, wl_soln_over_full)
co2ch4_shift=interpol(co2ch4_depth, tellwl, wl_soln_over_full)
if lab_depth eq 1 then int_lab_over_full=interpol(nh3_depth, wl_lab, wl_soln_over_full)

;;;;;;;
;multiply and convolve
product_model=stellar_template_over*int_lab_over_full*h2o_shift* co2ch4_shift
obs_noconv=product_model
;convolve

;;;;
;Convolve model with lsf
;;;;

int_conv=dblarr(npix_model)
int_conv_mtx=dblarr(npix_lsf, npix_model)


for index=0,npix_model-1 do begin
    model_index=index+first_pix*oversamp
    if model_index lt (npix_lsf/2) or model_index gt (n_elements(obs_noconv)-1-(npix_lsf/2)) then $
      int_conv_mtx[*,index]=replicate(1d0, npix_lsf) $
    else int_conv_mtx[*,index]=obs_noconv[model_index-(npix_lsf/2):model_index+(npix_lsf/2)]
endfor
int_conv_temp=int_conv_mtx*lsf
int_conv=total(int_conv_temp, 1, /double)


;;;;
;Down-sample the model to IRCS resolution
;;;;

tophat=replicate(1d0/oversamp,oversamp)
int_avg=convol(int_conv, tophat)

stellar_template_avg=convol(stellar_template_over, tophat)

int_lab_over_avg=convol(int_lab_over_full, tophat)
samp_index_full=lindgen(n_elements(wl_soln_full))*oversamp+(oversamp/2)
stellar_template=stellar_template_avg[samp_index_full]
int_lab_down=int_lab_over_avg[samp_index_full]
h2o_shift_avg=convol(h2o_shift, tophat)
h2o_shift_down=h2o_shift_avg[samp_index_full]
h2o_shift_dim=h2o_shift_down[first_pix:first_pix+npix_select-1]
co2ch4_shift_avg=convol(co2ch4_shift, tophat)
co2ch4_shift_down=co2ch4_shift_avg[samp_index_full]
co2ch4_shift_dim=co2ch4_shift_down[first_pix:first_pix+npix_select-1]
stellar_template_dim=stellar_template[first_pix:first_pix+npix_select-1]
int_lab_dim=int_lab_down[first_pix:first_pix+npix_select-1]


int_model=int_avg[samp_index] 


;;;;
;Correct subtle variations in continuum
;;;;
 
if norm gt 0 then begin    
    
    ;X vector of same dimensions as ircs-resolution spectrum
    norm_pts_xx=lindgen(npix_select)+first_pix
    
    ;Create and populate vector with x-coord of node pts
    norm_pts_x=dblarr(nnorm)
    
    if nnorm ge 2 then begin
        norm_pts_x[0]=first_pix
        norm_pts_x[nnorm-1]=first_pix+npix_select-1
        if nnorm ge 3 then begin
            x_increment=npix_select/(nnorm-1)
            for i=1, nnorm-2 do norm_pts_x[i]=norm_pts_x[i-1]+x_increment
        endif
        
        ;Interpolate between nodes with a spline
        if npix_select ge 128L and nnorm ge 5 then norm_factor=interpol(norm_pts_y, norm_pts_x, norm_pts_xx, /spline) $
        else norm_factor=interpol(norm_pts_y, norm_pts_x, norm_pts_xx)
    endif else if nnorm eq 1 then norm_factor=replicate(norm_pts_y, npix_select) 

int_model=int_model*norm_factor   

endif
    




;;;;
;;;Get chi^2
;;;;

int_obs_dim=int_obs[first_pix:first_pix+npix_select-1]
err_dim=err[first_pix:first_pix+npix_select-1]
nan=where(finite(int_model) eq 0, nancount)
if nancount gt 0 then begin
    err_dim[nan]=1d-5
    int_model[nan]=1d0
endif

res=(int_obs_dim-int_model)
res2=res^2
dev=res/err_dim
dev2=dev^2

chi1_vec=dev;[npix_trim_start:-1*(npix_trim_end+1)]
chi2_vec=dev2;[npix_trim_start:-1*(npix_trim_end+1)]
chi2=total(chi2_vec, /double, /nan)



;;;;
;If visualize is set, do animation of fitting spectrum for full spectrum
;;;;

if (visualize eq 1) then begin
    
    ;title_str="RV shift: "+strtrim(delta_rv,2)+" | tell_rv: " + strtrim(tell_rv,2)+" | visit: " + strtrim(visit,2) + " | chi2: "+strtrim(chi2,2)
;    title_str="RV shift: "+strtrim(delta_rv,2)+" | visit: " + strtrim(visit,2)+" | depth: " + strtrim(temp_depth_coeff,2) + " | chi2: "+strtrim(chi2,2)

;Set up big dots
phi=findgen(32)*(!PI*2/32.)
phi = [ phi, phi(0) ]
usersym, cos(phi), sin(phi), /fill

;Set up plot environment
;;Make Postscript
;set_plot, 'ps
;device, /encapsulated, filename = 'rv_fit.eps'
;device, /color, bits=8
;loadct, 12

set_plot, 'x'


;Plot one Butler-style graph

if finalplot eq 1 then begin
    !p.multi=0
    ;Make Postscript
    set_plot, 'ps'
    device, /encapsulated, filename = 'nessf.eps'
    device, /color, bits=8
    loadct, 12

    ind=lindgen(500)+10
     plot, wl_soln[ind], int_lab_dim[ind] + 2.0d0, yr=[0,3.2], /ys, /xs, ytickname=replicate(' ', 5), yticks=1,yminor=1, xtitle="Wavelength (microns)", charsize=1.2 , title="Ammonia Cell RV Modeling Process"
     oplot, wl_soln[ind], stellar_template_dim[ind]+1.2d0
     oplot, wl_soln[ind], (h2o_shift_dim[ind]*co2ch4_shift_dim[ind])+0.6d0
;     oplot, wl_soln, co2ch4_shift_dim+0.6d0, color=200
     oplot, wl_soln[ind], int_obs_dim[ind]
     oplot, wl_soln[ind], int_model[ind], color=200, ps=8, symsize=0.3 
     oplot, wl_soln[ind], res[ind] + 0.25, ps=3

     device, /close
     set_plot, 'x'
     !p.multi=[0,1,5]
endif

;Plot multiplot of relevant spectra

plot, wl_soln, int_lab_dim, /xs, yr=[0.3, 1.1]

plot, wl_soln, stellar_template_dim, /xs ,yr=[0.3, 1.1]

plot, wl_soln, h2o_shift_dim*co2ch4_shift_dim, /xs, yr=[0.3, 1.1]

plot, wl_soln, int_obs_dim, /xs, yr=[0.3, 1.1]
oplot, wl_soln, int_model, ps=8, color=200, symsize=0.5

;plot, res, /xs, ps=3, yr=[-0.1, 0.1]
plot, dev, /xs, ps=3

;Close Postscript device
;device, /close
    
endif

;stop

;;;;
;Report results
;;;;



print, chi2
print, last_guess
print, "----------"


if min_type eq 'amoeba' then return, chi2 $
else if min_type eq 'mpfit' then return, chi1_vec $
else message, 'incorrect min_type... cannot return a chi2 value'


end


pro ircsrv_maketemplate_ksjul13, epoch=epoch, object=object, trace=trace, visualize=visualize, first_pix=first_pix, npix_select=npix_select, min_type=min_type, lab_depth=lab_depth, norm=norm, model_tag=model_tag, telluric_option=telluric_option, temp_file=temp_file, chi2_tol=chi2_tol, smooth_opt=smooth_opt

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;MAIN BODY ;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Define keywords that I just added
if n_elements(epoch) eq 0 then epoch='18Jan2011'
if n_elements(object) eq 0 then object='GJ273'
if n_elements(trace) eq 0 then trace='AB'
if n_elements(model_tag) eq 0 then model_tag=strjoin((strsplit(systime(),' ',/extract))[1:2])+'temporary'
if n_elements(telluric_option) eq 0 then telluric_option=0 ;original telluric

if n_elements(temp_file) ne 0 then begin
    template_start=1
    if file_test(temp_file) eq 0 then begin
        print, temp_file + " does not exist.  Using phoenix_template.fits"
        temp_file='../data/epoch/18Jan2011/temp_results/phoenix_template.fits'
    endif
endif else template_start=0
if n_elements(chi2_tol) eq 0 then chi2_tol=0.07
if n_elements(smooth_opt) eq 0 then smooth_opt=0

;CAUTION: I'm fixing the following variable temporarily
;npix_select=150L

npix=1024L
if n_elements(visualize) eq 0 then visualize=0
;Define keywords for fitting just a range
if n_elements(first_pix) eq 0 then first_pix=0
if n_elements(npix_select) eq 0 then npix_select=npix-first_pix
if n_elements(min_type) eq 0 then min_type='mpfit'
if n_elements(lab_depth) eq 0 then lab_depth=0
if n_elements(norm) eq 0 then norm=0


;set file paths
rootpath='/home/stgilhool/RV_projects/IRCS_rv/data/'
epochpath=rootpath+'epoch/'+epoch+'/'
objectpath=epochpath+'final_spectra/'
flatpath=epochpath+'final_spectra/'
calibpath=epochpath+'calib_results/'
outpath=epochpath+'temp_results/'
modelpath=rootpath+'supplemental/'
outputpath=epochpath+'rv_results/'

;p=paths(epoch=epoch, object=object, trace=trace)
;rootpath=p.dataroot
;epochpath=p.epochfold
;objectpath=p.specfold
;flatpath=p.specfold
;outpath=p.tempfold
;modelpath=p.dataroot+'supplemental/'


;Define flat file names and full path file names
ABflat_filename=strjoin(['NH3', object, epoch, 'AB1'], '_')+'.fits'
;BAflat_filename=strjoin(['NH3', object, epoch, 'BA1'], '_')+'.fits'
ABflat_file=flatpath+strjoin(['NH3', object, epoch, 'AB1'], '_')+'.fits'
;BAflat_file=flatpath+strjoin(['NH3', object, epoch, 'BA1'], '_')+'.fits'


;Define object file names and full path file names
ABobj_listname=strjoin([object, epoch, 'AB'], '_')+'.list'
;BAobj_listname=strjoin([object, epoch, 'BA'], '_')+'.list'
ABobj_list=objectpath+strjoin([object, epoch, 'AB'], '_')+'.list'
;BAobj_list=objectpath+strjoin([object, epoch, 'BA'], '_')+'.list'


;Read in object filenames and then spectra
readcol, ABobj_list, ABobj_filename, format='A'
;readcol, BAobj_list, BAobj_filename, format='A'

ABobj_file=objectpath+ABobj_filename
;BAobj_file=objectpath+BAobj_filename

n_ABobj=n_elements(ABobj_file)
;n_BAobj=n_elements(BAobj_file)


;Initialize arrays of spectra
temp_str=mrdfits(ABobj_file[0], 1)
temp_spectrum=temp_str.spectrum
npix=n_elements(temp_spectrum)
ABspec_arr=dblarr(n_ABobj, npix)
;BAspec_arr=dblarr(n_BAobj, npix)
ABflat_arr=dblarr(n_ABobj, npix)
;BAflat_arr=dblarr(n_BAobj, npix)


ABmjd_arr=dblarr(n_ABobj)
;BAmjd_arr=dblarr(n_BAobj)

;define error arrays
ABerr_arr=dblarr(n_ABobj, npix)
;BAerr_arr=dblarr(n_BAobj, npix)

;Set continuum normalization params
low_reject=0.53
high_reject=3.0


;Loop through each object file and flip and normalize it
for f=0, n_ABobj-1 do begin
    ABobj_spec=mrdfits(ABobj_file[f], 1)
    ;BAobj_spec=mrdfits(BAobj_file[f], 1)

    ABspec_backwards=ABobj_spec.spectrum
    ;BAspec_backwards=BAobj_spec.spectrum

    ;Flip the spectra
    ABspec=reverse(ABspec_backwards)
    ;BAspec=reverse(BAspec_backwards)

        ;Read in and flip errors
    ABerr=reverse(ABobj_spec.sigma)
    ;BAerr=reverse(BAobj_spec.sigma)

    ;Get MJD
    ABhead=ABobj_spec.header
    ;BAhead=BAobj_spec.header

    ;Normalize the model spectrum and the observed spectrum
    ABnorm=continuum_fit(dindgen(npix), ABspec, low_rej=low_reject, high_rej=high_reject)
    ;BAnorm=continuum_fit(dindgen(npix), BAspec, low_rej=low_reject, high_rej=high_reject)
    
    ;Populate arrays
    ;spectra
    ABspec_arr[f,*]=ABspec/ABnorm
    ;BAspec_arr[f,*]=BAspec/BAnorm
    ;errors
    ABerr_arr[f,*]=ABerr/ABnorm
    ;BAerr_arr[f,*]=BAerr/BAnorm
    ;MJD
    ABmjd_arr[f]=sxpar(ABhead, 'MJD')
    ;BAmjd_arr[f]=sxpar(BAhead, 'MJD')

endfor

;Read flat files and flip and normalize them too    

; ABflat_spec=mrdfits(ABflat_file, 1)
; BAflat_spec=mrdfits(BAflat_file, 1)

; ABflat_backwards=ABflat_spec.NH3_spectrum
; BAflat_backwards=BAflat_spec.NH3_spectrum

; Flip
; ABflat=reverse(ABflat_backwards)
; BAflat=reverse(BAflat_backwards)

; Normalize
; ABflatnorm=continuum_fit(dindgen(npix), ABflat, low_rej=low_reject, high_rej=high_reject)
; BAflatnorm=continuum_fit(dindgen(npix), BAflat, low_rej=low_reject, high_rej=high_reject)

; Make array
; ABflat_arr=rebin(reform(ABflat/ABflatnorm, 1, npix), n_ABobj, npix)
; BAflat_arr=rebin(reform(BAflat/BAflatnorm, 1, npix), n_BAobj, npix)


; ;;;;;;;
; Divide spec by flat
; ABtemp1=ABspec_arr/ABflat_arr
; BAtemp1=BAspec_arr/BAflat_arr

; Take median
; ABtemp_rough=median(ABtemp1, dimension=1, /double)
; BAtemp_rough=median(BAtemp1, dimension=1, /double)

;;;;;;;

;Read in calibration results
calib_file='GJ273_18Jan2011_AB1_5_7_7_13.fits'
calib_ext=13
model_file=calibpath+calib_file
model_par=mrdfits(model_file, calib_ext)


common rvinfo, wl_soln, wl_soln_full, samp_index, wl_soln_over_full, int_lab_over_full, int_obs, lsf, npix_lsf, npix_model, firstpix, npixselect, oversamp, nnorm, err, visual, npix_trim_start, npix_trim_end, mintype, nparam_other, last_guess, nlines, COwl_shift, h2o, co2ch4, tellwl, visit, iter, template_residuals, res, stellar_template, stellar_template_over, finalplot, int_lab, wl_lab, lab_depth_switch, norm_switch



lab_depth_switch=lab_depth
norm_switch=norm
visual=visualize
npixselect=npix_select
firstpix=first_pix
mintype=min_type

;Turn off finalplot
finalplot=0

;Assume lin_switch is on
lin_switch=1

;;;Parameters from amoeba
wl_coeff=model_par.wl_result

gh_coeff=model_par.gh0_result
gh_coeff_out=gh_coeff
sigma=gh_coeff[0]
;gh_coeff[0]=1d0


if lin_switch eq 1 then begin
    gh_lin_coeff=model_par.gh1_result
    gh_lin_coeff_out=gh_lin_coeff
    ;sigma_lin=gh_lin_coeff[0]
    ;gh_lin_coeff[0]=0d0
    other=model_par.other_result
endif else begin
    ;sigma_lin=0d0
endelse

tau_scale=other[0]
if n_elements(other) gt 1 then norm_pts=other[1:n_elements(other)-1]

nparam=model_par.nparam
nparam_wl=n_elements(wl_coeff)
nparam_gh=n_elements(gh_coeff)
nparam_other=n_elements(other)
nparam_gh1=nparam-nparam_wl-nparam_gh-nparam_other


;if n_elements(other) gt 1 then begin
;    nparam_norm=nparam_other-1
;    x_x=lindgen(npix_select)+first_pix
;    x_norm=dblarr(nparam_norm)
;    x_norm[0]=first_pix
;    x_norm[nparam_norm-1]=first_pix+npix_select-1
;    x_increment=npix_select/(nparam_norm-2)
;    for i=1, nparam_norm-2 do x_norm[i]=x_norm[i-1]+x_increment

;    norm_factor=interpol(norm_pts, x_norm, x_x, /spline)    
    
;endif
;;;

;Other parameters and constants
oversamp=model_par.oversamp
npix_lsf=(oversamp*10L)+1L
c=299792458.D

;Grab other structure variables just in case
;npix_trim_start=model_par.npix_trim_start
;npix_trim_end=model_par.npix_trim_end
npix_trim_start=5L
npix_trim_end=5L

npix_model=npix_select*oversamp







;Read in LAB SPECTRUM
modelfile= modelpath+'NH3_model.dat'
readcol, modelfile, wl_lab, int_lab, format='D,D'

norm_lab=continuum_fit(wl_lab, int_lab, low_rej=low_reject, high_rej=high_reject)
int_lab_copy=int_lab
int_lab=int_lab_copy/norm_lab


;Adjust the optical depth of the lab spectrum (I think this is okay to
;do for the high-resolution spectrum even though the wl scale isn't linear)
if lab_depth eq 0 then int_lab_depth=int_lab^tau_scale



;read in telluric spectrum
;tell=mrdfits(modelpath+'skytable.fits', 1)
;tellwl=tell.lam
;telluric=tell.trans


if telluric_option eq 0 then begin
    h2ostr=mrdfits(modelpath+'sky_h2o.fits',1)
    co2ch4str=mrdfits(modelpath+'sky_co2andch4.fits', 1)
    tell_wl_long=h2ostr.wave
    m24index=where(tell_wl_long gt 2.275 and tell_wl_long lt 2.365)
    
    h2o_long=h2ostr.trans
    co2ch4_long=co2ch4str.trans
endif else begin

    h2ostr=mrdfits(modelpath+'tapas_h2ovac.fits',1)
    co2ch4str=mrdfits(modelpath+'tapas_co2ch4vac.fits', 1)

;    h2ostr=mrdfits(modelpath+'tapas_h2onovac.fits',1)
;    co2ch4str=mrdfits(modelpath+'tapas_co2ch4novac.fits', 1)


    tell_wl_long=h2ostr.wavelength/1000d0
    tell_wl_long=reverse(tell_wl_long)
    
    m24index=where(tell_wl_long gt 2.275 and tell_wl_long lt 2.365)
    
    h2o_long=h2ostr.transmittance
    h2o_long=reverse(h2o_long)

    co2ch4_long=co2ch4str.transmittance
    co2ch4_long=reverse(co2ch4_long)
endelse

tellwl=tell_wl_long[m24index]
h2o=h2o_long[m24index]
co2ch4=co2ch4_long[m24index]


;;;Contstruct wl array and oversampled wl array from coeffs
x=dindgen(npix)
xx=(dindgen(npix*oversamp)-(oversamp/2))/oversamp

wl_soln=poly(x, wl_coeff)
;samp_index=(lindgen(npix)*oversamp)+(oversamp/2)

wl_soln_over=poly(xx, wl_coeff)

;Interpolate the lab spectrum onto the oversampled trial grid
if lab_depth eq 0 then int_lab_over_full=interpol(int_lab_depth, wl_lab, wl_soln_over)

;Save full length vectors
wl_soln_full=wl_soln
wl_soln_over_full=wl_soln_over


;Redefine those vectors to reflect only the specified range


x_select=dindgen(npix_select)+first_pix
xx_select=(dindgen(npix_model)-(oversamp/2))/oversamp+first_pix

wl_soln=poly(x_select, wl_coeff)
samp_index=(lindgen(npix_select)*oversamp)+(oversamp/2)
wl_soln_over=poly(xx_select, wl_coeff)








;;;;
;Make LSF
;;;;

lsf=ircsrv_lsfold(gh_coeff, gh1_vec=gh_lin_coeff, oversamp=oversamp, npix_select=npix_select, first_pix=first_pix, /verbose)

;Now fit each observation
;n_obs=n_elements(spectra[0,*])

;for visit=0, n_ABobj-1 do begin
for visit=0, n_ABobj-1 do begin
    
;Define Template
;temp=ABtemp_rough
    bigc=299792458d0 ;m/s
    
;get header
    struct=mrdfits(ABobj_file[visit], 1)
    head=struct.header
;get bcv correction
    bcvcorr_ircs, head, params
    bcv=params[0]
    if visit eq 0 then bcv0=bcv
    delta_bcv=bcv0-bcv

;read in observation
    int_obs=ABspec_arr[visit,*]
    ;date and
    mjd=ABmjd_arr[visit]
    ;errors
    err=ABerr_arr[visit, *]
    ;don't count first and last few
    err[0:9]=1d8
    err[-10:-1]=1d8
    ;or the middle few
    err[502:521]=1d8
    
    mask_ind=where(err ge 1d8, bcount)
    
    
;directly compare them using AMOEBA
;;;;;;;;;;;;;;;;;;;;;;
;;;  Run Modeling Function   ;;;;
;;;;;;;;;;;;;;;;;;;;;;
    
    clock=tic()
    
;Define inputs to modeling function
    ftol=1d-11
    
    if telluric_option eq 0 then begin
        co2ch4_depth_guess=[0.3d0]
        h2o_depth_guess=[0.3d0]
    endif else begin
        co2ch4_depth_guess=[0.05d0]
        h2o_depth_guess=[0.5d0]
    endelse
    
    co2ch4_depth_scale=co2ch4_depth_guess;[0.45d0]
    h2o_depth_scale=h2o_depth_guess;[0.45d0]
    
    guess=[co2ch4_depth_guess, h2o_depth_guess]
    scale=[co2ch4_depth_scale, h2o_depth_scale]
    
    if lab_depth eq 1 then begin
        nh3_depth_guess=tau_scale
        nh3_depth_scale=0.1*tau_scale
        guess=[guess, nh3_depth_guess]
        scale=[scale, nh3_depth_scale]
    endif
    
    
    if norm ne 0 then begin
        nnorm=2 > npix_select/50
        norm_guess=replicate(1d0,nnorm)
        norm_scale=replicate(0.1d0,nnorm)
        
        guess=[guess, norm_guess]
        scale=[scale, norm_scale]
    endif else nnorm=0
    
    
    
    if visualize eq 1 then begin
        set_plot, 'x'
        window, 0, xsize=1500, ysize=650
    endif
    !p.multi=[0,1,5]
;!p.multi=0
;;;Run Minimization scheme
    
;read in last guess
;guess=[-2799.7870,-1081.8196,      0.99741826,      0.20462968,      0.12492732,      0.40048046,   8.1290099d-05,   0.00011054021,   2.2794023d-05]
;r=rvfit(guess)
;stop
    
    nparam_total=n_elements(guess)
    
    parinfo = replicate({limited:[1,1], $
                         limits:[0.D,1.D], mpside:0}, nparam_total)
;parinfo[2:2+nnorm-1].limits=[0.85d0,1.15d0]            

;READ IN PHOENIX TEMPLATE AS A STARTING GUESS
if template_start ne 0 then begin
    ttt=mrdfits(temp_file, 1)
    
    stellar_template_over=ttt.template
    stellar_wl=ttt.wl_soln
    
    ;Shift phoenix template to bcv rest frame
    wl_shifted=stellar_wl*(1d0+(delta_bcv*1000d0/bigc))
    stellar_template_over=interpol(stellar_template_over, wl_shifted, wl_soln_over_full)
    
    ;Downsample to IRCS res
    tophat=replicate(1d0/oversamp,oversamp)
    stell_avg=convol(stellar_template_over, tophat)
    samp_index_full=lindgen(n_elements(wl_soln_full))*oversamp+(oversamp/2)
    stellar_template=stell_avg[samp_index_full]
    
endif else stellar_template_over=replicate(1d0, n_elements(wl_soln_over))
best_template=stellar_template_over
    iter=0
    lowprobcount=0

;START RECORDING
openw, lun1, '../data/ksresults_'+model_tag+'_'+strtrim(first_pix,2)+'_'+strtrim(visit,2)+'.txt', /get_lun
    printf, lun1, '______________________________________'
    printf, lun1, 'Visit: ' + strtrim(visit,2)
    printf, lun1, 'Iter | Chi2 | Chi2_frac | prob | D'
    printf, lun1, ''



repeat begin ;    for iter=0, 8 do begin
        
        
;if iter eq 9 then begin 
;    finalplot=1
;    epscall=rvfit(r)
;    stop
;endif
        
;Make stellar template on oversampled wavelength grid
    if iter ne 0 then begin
        stellar_template_over=interpol(template_residuals, wl_soln_full, wl_soln_over_full)
        
        ;snorm=continuum_fit(wl_soln_over_full, stellar_template_over, low_rej=low_reject, high_rej=high_reject)
        ;stellar_template_over=stellar_template_over/snorm
    endif
    
    if min_type eq 'amoeba' then begin
        r=amoeba3(ftol, scale=scale, p0=guess,function_name='rvfit', $
                  function_value=fval, nmax=150000L)
        chi2=fval[0]
    endif else if min_type eq 'mpfit' then begin
        r=mpfit('rvfit', guess, parinfo=parinfo, bestnorm=chi2, ftol=ftol, /quiet)
    endif



    
    
    ;   !p.multi=0
    
    ;stop clock
    process_time=toc(clock)
        
        
        
    if iter eq 0 then begin
        ;initialize template_residuals vector
        template_residuals=stellar_template
        chi2_last=2*chi2
    endif
    
    chi2_frac=(chi2_last-chi2)/chi2_last
    chi2_last=chi2
    
    template_residuals_dim=template_residuals[first_pix:first_pix+npix_select-1]
    if smooth_opt eq 1 then smoothres=poly_smooth(res, gh_coeff[0]*4) else smoothres=res

    template_residuals_dim=template_residuals_dim+smoothres
    
    ;Check KS statistic
    pixel_x=first_pix+lindgen(npix_select)
    err_y=err[pixel_x]
    rmask_ind=where(err_y le 1d8, gcount)
    if gcount ne 0 then begin
        resdist=smoothres[rmask_ind]
        errdist=err_y[rmask_ind]    
        pixel_x=pixel_x[rmask_ind]
    endif else begin
        resdist=smoothres
        errdist=err_y
    endelse

    normresdist=resdist/errdist
    
    normdist=randomn(seed, n_elements(normresdist))

    kstwo, normresdist, normdist, D, prob

    if iter eq 0 then begin
        best_prob=prob
        best_iter=iter
        best_chi2=chi2
        best_chi2_frac=chi2_frac
    endif
    if prob ge best_prob then begin
        lowprobcount=0
        best_prob=prob
        best_template=stellar_template_over
        best_iter=iter
        best_chi2=chi2
        best_chi2_frac=chi2_frac
    endif else lowprobcount++
        
    printf, lun1, iter, chi2, chi2_frac, prob, D

;    plot, pixel_x, normresdist, title='Iter: '+strtrim(iter,2)+' | chi2: ' + strtrim(chi2,2) + ' | chi2_frac: '+ strtrim(chi2_frac,2) + ' | prob: ' + strtrim(prob, 2) + ' | D: ' + strtrim(D,2), ps=6, symsize=0.8, /xs, charsize=3.0

;    stop

    ;!p.multi=[0,1,2]
    ;plot, res, ps=6, symsize=0.5, title='Iteration: ' + strtrim(iter,2), charsize=2.5
    ;oplot, smoothres, ps=6, symsize=0.5, color=200
    
    ;plot, template_residuals_dim, title='chi2: ' + strtrim(chi2,2), charsize=2.5
    terminate=0
    
    if (chi2_frac lt chi2_tol or lowprobcount ge 10) and iter gt 4 then begin 
;    if iter eq 35 then begin 
        print, iter
        n_iter=iter
        terminate=1
        ;stop
    endif
    
    ;!p.multi=[0,1,5]
    
    ;tnorm=continuum_fit(wl_soln[first_pix:first_pix+npix_select-1],
    ;template_residuals_dim, low_rej=low_reject, high_rej=high_reject)
    ;tnorm=max(template_residuals_dim[npix_trim_start:-1*(npix_trim_end+1)])
    ;template_residuals_dim=template_residuals_dim/tnorm
    template_residuals[first_pix:first_pix+npix_select-1]=template_residuals_dim
    iter++
    ;if iter eq 11 then terminate=1 
endrep until terminate 


close, lun1
free_lun, lun1

;    endfor 
    
; openw, tlun, 'template.dat', /get_lun
; for line=0, n_elements(stellar_template_over)-1 do begin
;     printf, tlun, wl_soln_over_full[line], stellar_template_over[line]
; endfor
; close, tlun
; free_lun, tlun
    
;stop
    
;OUTPUT RESULTS
    
    model_id=strjoin([object, epoch,'AB'+strtrim(visit,2), strtrim(first_pix,2), strtrim(first_pix+npix_select-1,2), model_tag],'_')
    
    if n_elements(r) eq 1 then begin
        print, "Minimization scheme failed to converge"
        result_vec=last_guess
        rvfinal=last_guess[0]
        chi2=-1
        
    endif else begin
        
        print, "---------------------------------"
        print, "Model ID: ", model_id
        print, "---------------------------------"
        
        print, chi2
        
        result_vec=r
        rvfinal=r[0]
    endelse
    
    print, "---------------------------------"
    print, "Processing took ", process_time, " seconds."
    
    
    
    outfile=model_id+".fits"
    
    outfinal=outpath+outfile
    
    if file_test(outfinal) then begin
        fits_info, outfinal, n_ext=prev_ext, /silent
        extension=prev_ext+1
    endif else extension=1
    
    
    output_str={OBJ:object, $
                EPOCH:epoch, $
                TRACE:trace, $
                CALIB_FILE:calib_file, $
                CALIB_EXT:calib_ext, $
                WL_COEFF:wl_coeff, $
                GH0_COEFF:gh_coeff_out, $
                GH1_COEFF:gh_lin_coeff_out, $
                OTHER_COEFF:other, $
                OVERSAMP:oversamp, $
                FTOL:ftol, $
                NPIX_TRIM_START:npix_trim_start, $
                NPIX_TRIM_END:npix_trim_end, $
                FIRST_PIX:first_pix, $
                NPIX_SELECT:npix_select, $
                TIME:process_time, $
                MIN_TYPE:min_type, $
                EXT:extension, $
                CHI2:best_chi2, $
                CHI2_FRAC:best_chi2_frac, $
                MJD:mjd, $
                DELTA_RV:rvfinal, $
                N_ITER:n_iter, $
;                TEMPLATE:stellar_template_over, $
                TEMPLATE:best_template, $
                  BEST_ITER:best_iter, $
                  BEST_PROB:best_prob, $
                RESULT_VEC:result_vec $
               }
    
;stop
    mwrfits, output_str, outfinal
    
    
endfor

end
