function calib_fit, p, fake_key=fake_key

common fit_info, parinfo_all, parinfo,wl_fts, spec_fts, spec_obs, npix_ircs, npix, fpix, mask, err, oversamp, fit_measure, visualize, n_wl, n_gh, n_d_gh, n_t, n_k


;start clock
calibfunc_tic=tic()

;Constants
npix_lsf=(oversamp*10L)+1L
c=299792458.D


;Save and scale params
raw_guess = p
pscale = parinfo.parscale[0]
pshift = parinfo.parscale[1]
scaled_guess = param_scale_up(raw_guess, pscale, pshift)

;;;Parameters from amoeba
parnames=parinfo.parname

                                ;;Get wl vectors
wli	= where(strmatch(parnames, 'wl[0-4]*') eq 1, n_wl)
wl_coeff= scaled_guess[wli]

                                ;;Get constant lsf params
ghi	= where(strmatch(parnames, 'gh[0-9]*') eq 1, n_gh)
gh 	= scaled_guess[ghi]
                                ;;Get lsf shift params
d_ghi	= where(strmatch(parnames, 'd_gh[0-9]*') eq 1, n_d_gh)
if n_d_gh gt 0 then begin
    lin_switch 	= 1
    d_gh	= scaled_guess[d_ghi]
endif else lin_switch = 0

                                ;;Get tau params
ti	= where(strmatch(parnames, 't\_*') eq 1, n_t)
tau_nh3 = scaled_guess[ti]

                                ;;Get norm params
ki	= where(strmatch(parnames, 'k[0-9]*') eq 1, n_k)
k	= scaled_guess[ki]
;;;


;Adjust the optical depth of the lab spectrum (I think this is okay to
;do for the high-resolution spectrum even though the wl scale isn't linear)
spec_depth_fts	=	spec_fts^tau_nh3[0]
spec_alias_obs	=	spec_obs

;;;Contstruct wl array and oversampled wl array from coeffs
npix_over		= npix * oversamp

x_full			= dindgen(npix_ircs)
xx_full			= (dindgen(npix_ircs*oversamp)-(oversamp/2))/oversamp

x			= dindgen(npix)+fpix
x_over			= (dindgen(npix_over)-(oversamp/2))/oversamp+fpix

wl_soln_full		= poly(x_full, wl_coeff)
wl_soln_over_full	= poly(xx_full, wl_coeff)

wl_soln			= poly(x, wl_coeff)
wl_soln_over		= poly(x_over, wl_coeff)

;Interpolate the lab spectrum onto the oversampled trial grid
spec_over_full_fts = interpol(spec_depth_fts, wl_fts, wl_soln_over_full)
;Oversample the obs
;spec_over_full_obs = interpol(spec_alias_obs, wl_soln_full, wl_soln_over_full)


;;;;
;Make LSF
;;;;

lsf_tic=tic()

if lin_switch eq 1 then begin
    lsf = ircsrv_lsf(gh, gh1_vec=d_gh, oversamp=oversamp, first_pix=fpix, npix_select=npix, neg_penalty=lsf_penalty)
endif else $
  lsf = ircsrv_lsf(gh, oversamp=oversamp, first_pix=fpix, npix_select=npix, neg_penalty=lsf_penalty)

lsf_t=toc(lsf_tic)
print, "LSF construction took : ", lsf_t, " seconds."


;;;;
;Convolve model with lsf
;;;;

convol_tic=tic()
spec_over_model = ircs_convolve(spec_over_full_fts, lsf, oversamp=oversamp, first_pix=fpix, npix_select=npix)
convol_t=toc(convol_tic)
;print, "Convolution took : ", convol_t, " seconds."



;;;;
;Down-sample the model to IRCS resolution
;;;;

spec_dim_model = downsample_tophat(spec_over_model, oversamp)

;Redefine vectors to reflect only the specified range

spec_dim_obs	=	spec_alias_obs[x]
err_dim		=	err[x]
mask_dim	=	mask[x]


;;;;
;Correct subtle variations in continuum
;;;;

norm_factor = ircsrv_normcorr(k, first_pix=fpix, npix_select=npix)
    
spec_dim_model  = spec_dim_model*norm_factor

;;;


;;;;
;;;Get chi^2
;;;;
;Penalize undefined pixels
nan_ind=where(finite(spec_dim_model) eq 0, nancheck)
if nancheck gt 0 then begin
    unmaskednan=where(err_dim[nan_ind] lt 1d3, nancount)
    nan_penalty=nancount
    spec_dim_model[unmaskednan]=1d0
endif else nan_penalty=0d0

;Penalize negative lsf


if lsf_penalty gt 1d-7 then lsf_penalty= 1d0+lsf_penalty

penalty=1d0-exp(-1d0*(lsf_penalty+nan_penalty))

res = (spec_dim_obs - spec_dim_model)

mask_ind  = where(mask_dim eq 1, mcount, complement=good_ind, ncomplement=ndata)
if mcount gt 0 then err_dim[mask_ind] = 1d10
    
dev = res/err_dim

if ndata gt 0 then begin
    dev_masked = dev[good_ind]
    res_masked = res[good_ind]
    err_masked = err[good_ind]
endif else message, "must have at least 1 valid data point!"

;Reject maximally deviant pixels???

chi1_vec	= dev * exp(penalty)
chi1_vec_nopen	= dev

chi2		= total(chi1_vec^2, /double, /nan)
chi2_nopen	= total(chi1_vec_nopen^2, /double, /nan)

;for ks test, must use unmasked normalized residuals only
normdist = randomn(seed, n_elements(dev_masked))

kstwo, dev_masked, normdist, D, prob



;;;;
;If visualize is set, do animation of fitting spectrum for full spectrum
;;;;

if (visualize eq 1) then begin

;Set up some stuff to help display lsf across spectrum
    x1=npix/6
    x2=npix/2
    x3=(5*npix/6)
    
    xx1=x1*oversamp+(oversamp/2)
    xx2=x2*oversamp+(oversamp/2)
    xx3=x3*oversamp+(oversamp/2)
    if lin_switch eq 1 then begin
        lsf1=lsf[*,xx1]
        lsf2=lsf[*,xx2]
        lsf3=lsf[*,xx3]
        
    endif else begin
        lsf1=lsf
        lsf2=lsf
        lsf3=lsf
    endelse
    
    max_lsf=max(lsf1)>max(lsf2)>max(lsf3)

    x1_vec=(dindgen(npix_lsf)-(npix_lsf/2))+x1+fpix
    x2_vec=(dindgen(npix_lsf)-(npix_lsf/2))+x2+fpix
    x3_vec=(dindgen(npix_lsf)-(npix_lsf/2))+x3+fpix

run=0
niter=0
    title_str="NH3 Model Fitting Process | run "+strtrim(run+1,2)+"/"+strtrim(niter,2)+" | chi2 : "+strtrim(chi2,2)

    plot, wl_soln, spec_dim_obs,yr=[0.1,1.1], /xs, title=title_str,xtitle="Wavelength (microns)", ytitle="Relative Instensity", charsize=2.0 
    oplot, wl_soln, spec_dim_model, color=200 

    plot, wl_soln, res, title="Residuals | RMS : ", xtitle="Wavelength (microns)", yr=[-0.1, 0.1], ps=3, /xs

    plot, x1_vec, lsf1, yr=[-0.05*max_lsf,1.1*max_lsf], xr=[fpix,fpix+npix-1], /xs
    oplot, x2_vec, lsf2
    oplot, x3_vec, lsf3
;stop
;plot_t=toc(plot_tic)

endif



;;;;
;Report results
;;;;

calibfunc_t=toc(calibfunc_tic)

print, parnames
print, scaled_guess
print, ''
print, ['Chi2','Chi2_nopen', 'D', 'prob'] 
print, chi2, chi2_nopen, D, prob
print, "----------"



if fit_measure eq 'chi2' then return, chi2 $
else if fit_measure eq 'chi_vec' then return, chi1_vec $
else if fit_measure eq 'residuals' then return, res $
else if fit_measure eq 'prob' then return, prob $
else if fit_measure eq 'D' then return, D $
else message, 'incorrect fit_measure... cannot return a value'


end



function ga_fit, p, np, funa=funa

common fit_info, parinfo_all, parinfo,wl_fts, spec_fts, spec_obs, npix_ircs, npix, fpix, mask, err, oversamp, fit_measure, visualize, n_wl, n_gh, n_d_gh, n_t, n_k   


ga_gen=tic()

fit_vec=dblarr(np)
for model_num=0, np-1 do begin
    fitness=calib_fit(p[model_num])
    if fitness lt 1 then fitness=1d10
    fit_vec[model_num]=fitness
endfor

print, minmax(fit_vec)

ga_gen_time=toc(ga_gen)

return, fit_vec

end



;;;;;;;;;;;;;;;;;;;
;;;;;MAIN BODY;;;;;
;;;;;;;;;;;;;;;;;;;
pro ircsrv_calibrate_jul15, parinfo_all, first_pix=first_pix, npix_select=npix_select, object=object, epoch=epoch, trace=trace, fmode=fmode, fit_measure=fit_measure, visualize=visualize, oversamp=oversamp, outfile=outfile

npix_ircs=1024L
if n_elements(fmode) eq 0 then fmode="mpfit"
if n_elements(fit_measure) eq 0 then begin
    if fmode eq 'mpfit' then fit_measure = 'chi_vec' else fit_measure = 'chi2'
endif
if n_elements(visualize) eq 0 then visualize=0
if n_elements(oversamp) eq 0 then oversamp=7L
if n_elements(first_pix) eq 0 then first_pix=0L
if n_elements(npix_select) eq 0 then npix_select=npix_ircs-first_pix 
if n_elements(object) eq 0 then object='GJ273'
if n_elements(epoch) eq 0 then epoch='18Jan2011'
if n_elements(trace) eq 0 then trace='AB1'
if n_elements(outfile) eq 0 then outfile='defaultout.fits'



common fit_info, parinfo_all1, parinfo, wl_fts, spec_fts, spec_obs, npix_ircs1, npix, fpix, mask, err, oversamp1, fit_measure1, visualize1, n_wl, n_gh, n_d_gh, n_t, n_k
fit_measure1=fit_measure
parinfo_all1=parinfo_all
oversamp1=oversamp
visualize1=visualize


npix_ircs1=npix_ircs
npix=npix_select
fpix=first_pix




;start clock
starttime=systime()
clock=tic()

;initialize paths
paths=ircsrv_paths(epoch=epoch, object=object, trace=trace)

ftspath = paths.sup 	;;;data/supplemental
obspath = paths.spec 	;;;ep/final_spectra
outpath = paths.calib	;;;ep/calib_results

                                ;;Define file names
ftsfile = ftspath+'NH3_model.dat'
obsfile = obspath+'NH3_'+strjoin([object, epoch, trace], '_')+'.fits'


;read in both spectra

readcol, ftsfile, wl_fts, spec_nonorm_fts, format='D,D'

obs_fits=mrdfits(obsfile, 1)

spec_nonorm_obs=reverse(obs_fits.nh3_spectrum)
err=reverse(obs_fits.nh3_sigma)
mask=reverse(obs_fits.nh3_mask)
mask[0:9]		= 1
mask[502L:521L]		= 1
mask[1014L:1023L]	= 1



;Normalize the model spectrum and the observed spectrum
low_reject=0.53
high_reject=3.0

norm_fts = continuum_fit(wl_fts, spec_nonorm_fts, $
                         low_rej=low_reject, high_rej=high_reject)
spec_fts=spec_nonorm_fts/norm_fts

x=dindgen(npix_ircs)
x1=dindgen(npix_ircs/2)
x2=dindgen(npix_ircs/2)+npix_ircs/2

spec_nonorm_obs1=spec_nonorm_obs[x1]
spec_nonorm_obs2=spec_nonorm_obs[x2]

norm_obs1=continuum_fit(x1, spec_nonorm_obs1, low_rej=low_reject, high_rej=high_reject)
norm_obs2=continuum_fit(x1, spec_nonorm_obs2, low_rej=low_reject, high_rej=high_reject)
norm_obs=[norm_obs1, norm_obs2]
spec_obs=spec_nonorm_obs/norm_obs

;scale the error by same amount
err=err/norm_obs


;directly compare them using AMOEBA
;;;;;;;;;;;;;;;;;;;;;;
;;;  Run Modeling Function   ;;;;
;;;;;;;;;;;;;;;;;;;;;;

;;;Define amoeba inputs
ftol=1d-10

freepar=where(parinfo_all1.fixed eq 0, nfree)

if nfree eq 0 then message, "Must have at least 1 free param" else begin
    parinfo=parinfo_all1[freepar]
    guess=parinfo.value
    amoeba_scale=(parinfo.value-parinfo.limits[0]) < (parinfo.limits[1]-parinfo.value)
endelse

if visualize eq 1 then window, 0, xsize=1000, ysize=400
!p.multi=[0,1,3]

;;;Run Minimization scheme


    if fmode eq 'amoeba' then begin
        ;fit_measure1='chi2'
        r = amoeba3(ftol, scale=amoeba_scale, p0=guess,function_name='calib_fit', $
                  function_value=fval, nmax=5000L)
        fitness=fval[0]
        
    endif else if fmode eq 'mpfit' then begin
        ;fit_measure1='chi_vec'
        r = mpfit('calib_fit', guess, parinfo=parinfo, bestnorm=fitness, ftol=ftol, /quiet)
    
    endif else if fmode eq 'ga' then begin
        limit=parinfo.limits
        funa={nothing:0}
        r = solber('ga_fit', nfree, funa=funa, lim=limit, gfit_best=fitness, term_flag=0, term_fit=-1, plot_flag=1, npop=1000L, ngen_max=3);, new_save_gen=save_gen, new_save_gfit=save_gfit, new_save_igen=save_igen, new_save_gbest=save_gbest, new_save_time=save_time, _extra=ex)
        
    endif 

!p.multi=0

;stop clock
process_time=toc(clock)
endtime=systime()

;OUTPUT RESULTS
model_id 	= strjoin(strtrim([n_wl,n_gh,n_d_gh,n_k],2), '_')

;Rescale results

if n_elements(r) eq 1 then begin
    parinfo.convflag	= 0
    parinfo.bestfit 	= raw_guess
    parinfo.tbestfit 	= scaled_guess
    parinfo_all1[freepar]=parinfo
    print, 'WARNING: Minimization scheme failed to converge!'

endif else begin
    r_scaled 		= param_scale_up(r, parinfo.parscale[0], parinfo.parscale[1])
    parinfo.convflag	= 1
    parinfo.bestfit	= r
    parinfo.tbestfit	= r_scaled
    parinfo_all1[freepar]=parinfo

endelse

print, "---------------------------------"
print, "Model ID: ", model_id
print, "---------------------------------"
print, "PARNAMES: "
print, parinfo.parname
print, ''
print, "GUESSES: "
print, parinfo.tvalue
print, ''
print, "RESULTS: "
print, parinfo.tbestfit
print, ''
print, "FITNESS MEASURE AND VALUE: "
print, fit_measure1, fitness

print, "---------------------------------"
print, "Processing took ", process_time, " seconds."



;outfile=model_id+".fits"
;run=0

; repeat begin
    
;     run++
    
;     outfile=model_id+"_"+strtrim(run,2)+".dat"

;     outfinal=outpath+outfile

; endrep until file_test(outfinal) eq 0

outfinal= outpath+outfile
outlast = outpath+'lastfit.fits'
outbest = outpath+model_id+'_best.fits'

output_str={MODEL_ID:model_id, $
            OBJ:object, $
            EPOCH:epoch, $
            TRACE:trace, $
            OVERSAMP:oversamp, $
            STARTTIME:starttime, $
            ENDTIME:endtime, $
            TIME:process_time,$ 
            FIRST_PIX:fpix, $
            NPIX_SELECT:npix, $
            NPARAM:nfree, $
            fmode:fmode, $
            FIT_MEASURE:fit_measure1, $
            FITNESS:fitness $
           }

mwrfits, output_str, outfinal, /create
mwrfits, parinfo_all1, outfinal
mwrfits, parinfo, outfinal

mwrfits, output_str, outlast, /create
mwrfits, parinfo_all1, outlast
mwrfits, parinfo, outlast


if file_test(outbest) eq 0 then begin

    mwrfits, output_str, outbest, /create
    mwrfits, parinfo_all1, outbest
    mwrfits, parinfo, outbest
    
endif else begin
    prevbest = mrdfits(outbest, 1)
    if prevbest.fitness gt fitness and total(parinfo.convflag) ne 0 then begin
        mwrfits, output_str, outbest, /create
        mwrfits, parinfo_all1, outbest
        mwrfits, parinfo, outbest
    endif
    
endelse

;save, output_str, file=outfinal

end






