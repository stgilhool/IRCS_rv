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

int_conv = ircs_convolve(obs_noconv, lsf, oversamp=oversamp, first_pix=first_pix, npix_select=npix_select)


;;;;
;Down-sample the model to IRCS resolution
;;;;

stellar_template	= downsample_tophat(stellar_template_over, oversamp)

int_model		= downsample_tophat(int_conv, oversamp)
int_lab_down		= downsample_tophat(int_lab_over_full, oversamp)
h2o_shift_down		= downsample_tophat(h2o_shift, oversamp)
co2ch4_shift_down	= downsample_tophat(co2ch4_shift, oversamp)

;And take just the selected pixel range

stellar_template_dim	= stellar_template[first_pix:first_pix+npix_select-1]
int_lab_dim		= int_lab_down[first_pix:first_pix+npix_select-1]
h2o_shift_dim		= h2o_shift_down[first_pix:first_pix+npix_select-1]
co2ch4_shift_dim	= co2ch4_shift_down[first_pix:first_pix+npix_select-1]



;;;;
;Correct subtle variations in continuum
;;;;

if norm gt 0 then begin    
    
    norm_factor 	= ircsrv_normcorr(norm_pts_y, first_pix=first_pix, $
                                          npix_select=npix_select)
 
    int_model		= int_model*norm_factor   

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


pro ircsrv_modeltemplate, epoch=epoch, object=object, trace=trace, visualize=visualize, first_pix=first_pix, npix_select=npix_select, min_type=min_type, lab_depth=lab_depth, norm=norm, model_tag=model_tag, telluric_option=telluric_option, temp_file=temp_file, chi2_tol=chi2_tol, smooth_opt=smooth_opt, calib_model=calib_model

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


if n_elements(calib_model) eq 0 then calib_model=calibpath+'5_2_0_6_best.fits'

ircsrv_readobs, object=object, epoch=epoch, trace=trace, spectra=spectra, errors=errors, mjds=mjds, headers=headers


common rvinfo, wl_soln, wl_soln_full, samp_index, wl_soln_over_full, int_lab_over_full, int_obs, lsf, npix_lsf, npix_model, firstpix, npixselect, oversamp, nnorm, err, visual, npix_trim_start, npix_trim_end, mintype, nparam_other, last_guess, nlines, COwl_shift, h2o, co2ch4, tellwl, visit, iter, template_residuals, res, stellar_template, stellar_template_over, finalplot, int_lab, wl_lab, lab_depth_switch, norm_switch

lab_depth_switch=lab_depth
norm_switch=norm
visual=visualize
npixselect=npix_select
firstpix=first_pix
mintype=min_type



;Read in calibration results
calibinfo_all = mrdfits(calib_model, 2)

used_ind = where(calibinfo_all.used eq 1, parcount)

if parcount eq 0 then message, "error reading in calibration results" else begin
    
    calibinfo	= calibinfo_all[used_ind]

    ind_str	= parinfo_getindex(calibinfo, /used)

    calib_r	= calibinfo.tbestfit

    n_wl	= ind_str.n_wl
    wli		= ind_str.wli
    n_gh	= ind_str.n_gh
    ghi		= ind_str.ghi
    n_d_gh	= ind_str.n_d_gh
    if n_d_gh gt 0 then begin
        lin_switch 	= 1
        d_ghi		= ind_str.d_ghi
        gh_lin_coeff	= calib_r[d_ghi]
    endif else begin
        lin_switch 	= 0
        gh_lin_coeff	= 0d0
    endelse
    n_t		= ind_str.n_t
    ti		= ind_str.ti
    n_k		= ind_str.n_k
    ki		= ind_str.ki

    
    wl_coeff	= calib_r[wli]
    gh_coeff	= calib_r[ghi]
    
    sigma	= gh_coeff[0]  
    tau_scale 	= calib_r[ti[0]]
    norm_pts 	= calib_r[ki]

endelse

;Turn off finalplot
finalplot=0

;Other parameters and constants
oversamp=model_par.oversamp
npix_lsf=(oversamp*10L)+1L
c=299792458.D

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


;Eliminate first observation
n_ABobj=n_ABobj-1
ABspec_arr=ABspec_arr[1:*,*]
ABerr_arr=ABerr_arr[1:*,*]
ABmjd_arr=ABmjd_arr[1:*,*]


;;;;
;Make LSF
;;;;

lsf=ircsrv_lsf(gh_coeff, gh1_vec=gh_lin_coeff, oversamp=oversamp, npix_select=npix_select, first_pix=first_pix, /verbose)

;Now fit each observation

for visit=0, n_ABobj-1 do begin
    
;Define Template
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
    ftol=1d-10
    
    if model_input eq 1 then begin
        model_in	= parinfo_readin(model_input)
    
    
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
    stell_avg=downsample_tophat(stellar_template_over, oversamp)
    samp_index_full=lindgen(n_elements(wl_soln_full))*oversamp+(oversamp/2)
    stellar_template=stell_avg[samp_index_full]
    
endif else stellar_template_over=replicate(1d0, n_elements(wl_soln_over))
best_template=stellar_template_over
    iter=0
    lowprobcount=0

;START RECORDING
openw, lun1, '../data/ksresults_'+model_tag+'_'+strtrim(first_pix,2)+'_'+strtrim(visit+1,2)+'.txt', /get_lun
    printf, lun1, '______________________________________'
    printf, lun1, 'Visit: ' + strtrim(visit+1,2)
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
    
    model_id=strjoin([object, epoch,'AB'+strtrim(visit+1,2), strtrim(first_pix,2), strtrim(first_pix+npix_select-1,2), model_tag],'_')
    
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
