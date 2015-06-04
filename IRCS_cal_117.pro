function sg_hermite_coeff

coeff=dblarr(11, 11)
coeff[0,0]=1d0
coeff[1,1]=2d0
coeff[2,2]=4d0
coeff[0,2]=-2d0
coeff[1,3]=-12d0
coeff[3,3]=8d0
coeff[0,4]=12d0
coeff[2,4]=-48d0
coeff[4,4]=16d0
coeff[1,5]=120d0
coeff[3,5]=-160d0
coeff[5,5]=32d0
coeff[0,6]=-120d0
coeff[2,6]=720d0
coeff[4,6]=-480d0
coeff[6,6]=64d0
coeff[1,7]=-1680d0
coeff[3,7]=3360d0
coeff[5,7]=-1344d0
coeff[7,7]=128d0
coeff[0,8]=1680d0
coeff[2,8]=-13440d0
coeff[4,8]=13440d0
coeff[6,8]=-3584d0
coeff[8,8]=256d0
coeff[1,9]=30240d0
coeff[3,9]=-80649d0
coeff[5,9]=48384d0
coeff[7,9]=-9216d0
coeff[9,9]=512d0
coeff[0,10]=-30240d0
coeff[2,10]=302400d0
coeff[4,10]=-403200d0
coeff[6,10]=161280d0
coeff[8,10]=-23040d0
coeff[10,10]=1024d0


return, coeff

end


function continuum_fit, x_vec, int, low_rej=low_rej, high_rej=high_rej, PIX_MASK=mask, verbose=verbose

;Initialize arrays and loop variables
npix=n_elements(x_vec)
if n_elements(mask) eq 0 then mask=replicate(0,npix)
if n_elements(low_rej) eq 0 then low_rej=0.5
if n_elements(high_rej) eq 0 then high_rej=3.
if n_elements(verbose) eq 0 then verbose=0
n_iter=10
n_reject=0
n_reject_total=0
deg=5

;Apply pixel mask
unmask_index=where(mask eq 0)
x_vec=x_vec[unmask_index]
int=int[unmask_index]

;Copy vectors for iteration
x_vec_iter=x_vec
int_iter=int


;;;;Iteratively fit the spectrum to approximate continuum
for iteration=1,n_iter do begin
    
    ;Fit a robust 5th degree polynomial to data
    fit_coeff=robust_poly_fit(x_vec_iter, int_iter, deg, fit_iter, sig)
   
    ;Define upper and lower rejection bounds
    upper_bound=fit_iter+(high_rej*sig)
    lower_bound=fit_iter-(low_rej*sig)
    
    ;Find the indices where the data are within the bounds
    keep_index=where(int_iter le upper_bound and int_iter ge lower_bound)
    
    ;And keep track of the points that are rejected
    rej_index=where(int_iter gt upper_bound or int_iter lt lower_bound, rej_count)
    if rej_count gt 0 then begin
        n_reject=n_elements(rej_index)
        n_reject_total=n_reject_total+n_reject
    endif else n_reject=0

;     ;Stop iterating if code rejects half of the data or more
;                                 NOTE: This needs to be fixed so as to
;                                 output the last iteration's fit!
;     if n_reject_total ge 0.5*n_elements(x_vec) then begin
;         if verbose ne 0 then print, "Too many rejections. Exiting loop at iteration number: ", iteration
;         break
;         endif

    ;Calculate RMSE
    rmse=sqrt(mean((fit_iter-int_iter)^2))

    ;Display fit if verbose is set
    if verbose ne 0 then begin
                                ;Plot the fit and rejected points
        window, 1, xsize=1300, ysize=650
        plot, x_vec, int, /xs
        oplot, x_vec_iter, fit_iter, color=200
        if rej_count gt 0 then oplot, x_vec_iter[rej_index], int_iter[rej_index], ps=7
        
                                ;Print the results
        print, "Iteration number: ", iteration
        print, "Coefficients: ", fit_coeff
        print, "Sigma: ", sig
        print, "Number of rejected points (current iteration): ", n_reject
        print, "Number of rejected points (in total): ", n_reject_total
        print, "RMSE: ", rmse
        print, ""
        ;stop
    endif

    ;Change x_vec_iter and int_iter according to rejection
    x_vec_iter=x_vec_iter[keep_index]
    int_iter=int_iter[keep_index]

endfor

;create array to return
fit=dblarr(n_elements(int))
order=dblarr(n_elements(int))
for i=0,deg do begin
    order=fit_coeff[i]*(x_vec^i)
    fit=fit+order
endfor

;plot, x_vec, int, /xs
;oplot, x_vec, fit, color=200
;stop

return, fit

end


function amoebafunction, p



common amoeba_info, run, niter, visualize, oversamp, npix, npix_trim_start, npix_trim_end, model_func, dof, wl_lab, int_lab, int_obs, nparam_wl, nparam_gh, nparam_other, lin_switch, trial_lsf, err, H_coeff, min_type


amoebafunc_tic=tic()

;;;Parameters from amoeba
nparam=n_elements(p)
if lin_switch eq 1 then nparam_gh_lin=nparam_gh else nparam_gh_lin=0

wl_coeff=p[0:nparam_wl-1]

gh_coeff=p[nparam_wl:nparam_wl+nparam_gh-1]
sigma=gh_coeff[0]
gh_coeff[0]=1d0

if lin_switch eq 1 then begin
    gh_lin_coeff=p[nparam_wl+nparam_gh:nparam_wl+(2*nparam_gh)-1]
    sigma_lin=gh_lin_coeff[0]
    gh_lin_coeff[0]=0d0
    other=p[nparam_wl+(2*nparam_gh):nparam-1]
endif else begin
    other=p[nparam_wl+nparam_gh:nparam-1]
    sigma_lin=0d0
endelse

tau_scale=other[0]
if n_elements(other) gt 1 then norm_pts=other[1:n_elements(other)-1]


;Other parameters and constants
npix_lsf=71L
c=299792458.D



;Adjust the optical depth of the lab spectrum (I think this is okay to
;do for the high-resolution spectrum even though the wl scale isn't linear)
int_lab_depth=int_lab^tau_scale


;;;Contstruct wl array and oversampled wl array from coeffs
x=dindgen(npix)
xx=dindgen(npix*oversamp)

trial_soln=poly(x, wl_coeff)
samp_index=(indgen(npix)*oversamp)+(oversamp/2)
trial_soln_over=interpol(trial_soln, samp_index, xx)


;Interpolate the lab spectrum onto the oversampled trial grid
int_lab_over=interpol(int_lab_depth, wl_lab, trial_soln_over)





;;;;
;Make LSF
;;;;

lsf_tic=tic()
npix_model=n_elements(int_lab_over)
;populate lsf coordinates
x_lsf=rebin(dindgen(npix_lsf)-(npix_lsf/2), npix_lsf, npix_model)
pixel_2d=rebin(reform(dindgen(npix_model), 1, npix_model), npix_lsf, npix_model)
sig_arr=rebin(reform(sigma+(sigma_lin*dindgen(npix_model)), 1, npix_model), npix_lsf, npix_model)
x_sig_2d=rebin(x_lsf/sig_arr, npix_lsf, npix_model)
x_sig_3d=rebin(x_lsf/sig_arr, npix_lsf, npix_model, nparam_gh)

;Make non-normalized gh-polynomial cube (npix_lsf x npix_model x
;basis)

gh_cube_nonorm=dblarr(npix_lsf, npix_model, nparam_gh)
gh_cube_norm=dblarr(npix_lsf, npix_model, nparam_gh)
gh_cube=dblarr(npix_lsf, npix_model, nparam_gh)
for basis=0, nparam_gh-1 do begin
    ;make gh_polynomial for given basis
    gh_cube_nonorm[*,*,basis]=poly(x_sig_2d, H_coeff[*,basis])*exp(-0.5*(x_sig_2d^2))
    ;normalize
    norm=rebin(reform(total(abs(gh_cube_nonorm[*,*,basis]), 1, /double), 1, npix_model), npix_lsf, npix_model)
    gh_cube_norm[*,*,basis]=gh_cube_nonorm[*,*,basis]/norm
    ;Adjust amplitudes according to gh parameters
    if lin_switch eq 1 then basis_coeff=gh_coeff[basis]+(gh_lin_coeff[basis]*pixel_2d) $
      else if lin_switch eq 0 then basis_coeff=replicate(gh_coeff[basis], npix_lsf, npix_model)
    gh_cube[*,*,basis]=basis_coeff*gh_cube_norm[*,*,basis]
endfor

;make linear combination of bases
if nparam_gh gt 1 then trial_lsf_nonorm=total(gh_cube, 3, /double) $
  else trial_lsf_nonorm=gh_cube[*,*,0]

;account for unrealistically negative lsfs

min_vec=min(trial_lsf_nonorm, dimension=1)
neg_ind=where(min_vec lt 0, neg_count)
if neg_count gt 0 then begin
    neg_vec=dblarr(npix_model)
    neg_vec[neg_ind]=min_vec[neg_ind]

    neg_arr=rebin(reform(neg_vec, 1, npix_model), npix_lsf, npix_model)

    trial_lsf_nonorm_copy=trial_lsf_nonorm
    trial_lsf_nonorm=trial_lsf_nonorm-neg_arr
endif

;normalize lsf
lsf_norm=rebin(reform(total(trial_lsf_nonorm, 1, /double), 1, npix_model), npix_lsf, npix_model)
trial_lsf=trial_lsf_nonorm/lsf_norm
if lin_switch eq 0 then trial_lsf=trial_lsf[*,0]

lsf_t=toc(lsf_tic)
print, "LSF construction took : ", lsf_t, " seconds."



;decompose the lsf if model_func is set
if model_func ne 0 then begin
    lsf_norm_3d=rebin(lsf_norm, npix_lsf, npix_model, nparam_gh)

;decompose the lsf
    if neg_count gt 0 then begin
        neg_arr_3d=rebin(neg_arr, npix_lsf, npix_model, nparam_gh)
        lsf_decomp=(gh_cube-neg_arr_3d)/lsf_norm_3d
    endif else begin
        lsf_decomp=gh_cube/lsf_norm_3d
    endelse
endif




;;;;
;Convolve model with lsf
;;;;

convol_tic=tic()
if lin_switch eq 0 then begin
    int_lab_conv=convol(int_lab_over, trial_lsf) 
endif else begin
    int_lab_conv=dblarr(npix_model)
    int_lab_conv_mtx=dblarr(npix_lsf, npix_model)
    
    for index=0,npix_model-1 do begin
        if (index lt npix_lsf/2 or index gt (npix_model-1-(npix_lsf/2))) then $
            int_lab_conv_mtx[*,index]=replicate(1d0, npix_lsf) $
          else int_lab_conv_mtx[*,index]=int_lab_over[index-(npix_lsf/2):index+(npix_lsf/2)]
    endfor
    int_lab_conv_temp=int_lab_conv_mtx*trial_lsf
    int_lab_conv=total(int_lab_conv_temp, 1, /double)
endelse
        
convol_t=toc(convol_tic)
;print, "Convolution took : ", convol_t, " seconds."



;;;;
;Down-sample the model to IRCS resolution
;;;;

npix_tophat=oversamp
tophat=replicate(1d0/npix_tophat,npix_tophat)
int_lab_avg=convol(int_lab_conv, tophat)
int_model=int_lab_avg[samp_index] 



;;;;
;Correct subtle variations in continuum
;;;;

if n_elements(other) gt 1 then begin
    nparam_norm=nparam_other-1
    x_x=lindgen(npix)
    x_norm=dblarr(nparam_norm)
    x_norm[0]=0L
    x_norm[nparam_norm-1]=npix-1
    x_increment=npix/(nparam_norm-2)
    for i=1, nparam_norm-2 do x_norm[i]=x_norm[i-1]+x_increment

    norm_factor=interpol(norm_pts, x_norm, x_x, /spline)
    
    int_model=int_model*norm_factor
endif
;;;


;;;;
;;;Get chi^2
;;;;

chi_tic=tic()	;Start clock

res=(int_obs-int_model)
res2=res^2
dev=res/err
dev2=dev^2

chi1_vec=dev[npix_trim_start:-1*(npix_trim_end+1)]
chi2_vec=dev2[npix_trim_start:-1*(npix_trim_end+1)]
chi2=total(chi2_vec, /double, /nan)
chi2perdegree=chi2/dof

chi_t=toc(chi_tic);Stop clock



;;;;
;If visualize is set, do animation of fitting spectrum for full spectrum
;;;;

if (visualize eq 1) or (model_func ne 0) then begin

;Set up some stuff to help display lsf across spectrum
    x1=100L
    x2=500L
    x3=900L
    
    xx1=x1*oversamp+(oversamp/2)
    xx2=x2*oversamp+(oversamp/2)
    xx3=x3*oversamp+(oversamp/2)
    if lin_switch eq 1 then begin
        lsf1=trial_lsf[*,xx1]
        lsf2=trial_lsf[*,xx2]
        lsf3=trial_lsf[*,xx3]
        
    endif else begin
        lsf1=trial_lsf
        lsf2=trial_lsf
        lsf3=trial_lsf
    endelse
    
    max_lsf=max(lsf1)>max(lsf2)>max(lsf3)
    x1_vec=dindgen(npix_lsf)+x1
    x2_vec=dindgen(npix_lsf)+x2
    x3_vec=dindgen(npix_lsf)+x3
    
;plot_tic=tic()
    
    title_str="NH3 Model Fitting Process | run "+strtrim(run+1,2)+"/"+strtrim(niter,2)+" | chi2/DoF : "+strtrim(chi2perdegree,2)

    plot, trial_soln, int_obs,yr=[0.1,1.1], /xs, title=title_str,xtitle="Wavelength (microns)", ytitle="Relative Instensity", charsize=2.0 
    oplot, trial_soln, int_model, color=200 

    plot, trial_soln, res, title="Residuals | RMS : ", xtitle="Wavelength (microns)", yr=[-0.1, 0.1], ps=3, /xs

    plot, x1_vec, lsf1, yr=[-0.05*max_lsf,1.1*max_lsf], xr=[0,npix], /xs
    oplot, x2_vec, lsf2
    oplot, x3_vec, lsf3

;plot_t=toc(plot_tic)

endif



;;;;
;Report results
;;;;

amoebafunc_t=toc(amoebafunc_tic)
;print, "One full iteration of amoebafunc took : ",amoebafunc_t, " s"
;print, "LSF construction took ", (lsf_t/amoebafunc_t)*100d0, " percent of the time"
;print, "Convolution took ", (convol_t/amoebafunc_t)*100d0, " percent of the time"
;print, "Plotting took ", plot_t, " seconds and ", (plot_t/amoebafunc_t)*100d0, " percent of the time"

print, "Free params: ",nparam_wl, nparam_gh, nparam_gh_lin, nparam_other, min_type, " | Run : "+strtrim(run+1,2)+'/'+strtrim(niter,2)
print, p
print, chi2, chi2perdegree
print, "----------"

if model_func ne 0 then begin
    print, "Stopping for user analysis of model"
    print, "Define variable for outfilename"
    outfilename='tempout.fits'
    outstr = {lsf:trial_lsf, $
              lsf_decomp: lsf_decomp, $
              trial_soln: trial_soln, $
              trial_soln_over: trial_soln_over, $
              oversamp: oversamp, $
              int_model: int_model, $
              norm_factor: norm_factor, $
              int_lab_conv: int_lab_conv $
              }
    stop
    
    mwrfits, outstr, outfilename, /create
    print, "outfile written"
    stop

endif else begin
    if min_type eq 'amoeba' then return, chi2 $
    else if min_type eq 'mpfit' then return, chi1_vec $
    else message, 'incorrect min_type... cannot return a chi2 value'
endelse

end


;;;;;;;;;;;;;;;;;;;
;;;;;MAIN BODY;;;;;
;;;;;;;;;;;;;;;;;;;
pro IRCS_cal_117, wl_guess, wl_scale, gh_guess, gh_scale, gh_lin_guess, gh_lin_scale, other_guess, other_scale, run, niter, lin_switch=lin_switch, min_type=min_type, visualize=visualize, oversamp=oversamp, npix_trim_start=npix_trim_start, npix_trim_end=npix_trim_end, model_func=model_func


if n_elements(lin_switch) eq 0 then lin_switch=1
if n_elements(min_type) eq 0 then min_type="amoeba"
if n_elements(visualize) eq 0 then visualize=0
if n_elements(oversamp) eq 0 then oversamp=7L
if n_elements(npix_trim_start) eq 0 then npix_trim_start=50L
if n_elements(npix_trim_end) eq 0 then npix_trim_end=50L
if n_elements(model_func) eq 0 then model_func=0


common amoeba_info, runn, n_iter, visual, oversample, npix, npixtrim_start, npixtrim_end, modelfunc, dof, wl_lab, int_lab, int_obs, nparam_wl, nparam_gh, nparam_other, linear_switch, trial_lsf, err, H_coeff, mintype

;Define constants and stuff
version_number='1.16b'

npix=1024L
npixtrim_start=npix_trim_start
npixtrim_end=npix_trim_end
oversample=oversamp
linear_switch=lin_switch
runn=run
visual=visualize
modelfunc=model_func
n_iter=niter
mintype=min_type

starttime=systime()

H_coeff=sg_hermite_coeff()

;start clock
clock=tic()

;read in both spectra
readcol, 'NH3_model.dat', wl_lab, int_lab, format='D,D'
obs_fits=mrdfits('NH3_obs.fits', 1)
int_obs_nonorm_backwards=obs_fits.nh3_spec
err_backwards=obs_fits.nh3_err
mask_backwards=obs_fits.nh3_mask

;flip the spectrum to read low to high
int_obs_nonorm=reverse(int_obs_nonorm_backwards)
err=reverse(err_backwards)
mask=reverse(mask_backwards)

low_reject=0.53
high_reject=3.0
;Normalize the model spectrum and the observed spectrum
norm=continuum_fit(wl_lab, int_lab, low_rej=low_reject, high_rej=high_reject)
int_lab_copy=int_lab
int_lab=int_lab_copy/norm

wl_proxy=dindgen(n_elements(int_obs_nonorm))
norm_obs=continuum_fit(wl_proxy, int_obs_nonorm, low_rej=low_reject, high_rej=high_reject)
int_obs=int_obs_nonorm/norm_obs
;scale the error by same amount
err=err/norm_obs
mask_err=1d10*mask
err=err+mask_err


;directly compare them using AMOEBA
;;;;;;;;;;;;;;;;;;;;;;
;;;  Run Modeling Function   ;;;;
;;;;;;;;;;;;;;;;;;;;;;

;;;Define amoeba inputs
;ftol=1d-15
;ftol=1d-10
ftol=5d-11

;trial_lsf=dblarr(100)

;Numbers of parameters
nparam_wl=n_elements(wl_guess)
nparam_gh=n_elements(gh_guess)
nparam_other=n_elements(other_guess)
if lin_switch eq 0 then begin
    nparam_gh_lin=0 
    nparam_gh_tot=nparam_gh
endif else begin
    nparam_gh_lin=n_elements(gh_lin_guess)
    nparam_gh_tot=nparam_gh+nparam_gh_lin
endelse
nparam_tot=nparam_wl+nparam_gh_tot+nparam_other
dof=npix-nparam_tot-(npix_trim_start+npix_trim_end)

;;;Define guess and scale vectors
if lin_switch eq 1 then begin
    guess=[wl_guess, gh_guess, gh_lin_guess, other_guess]    
    scale=[wl_scale, gh_scale, gh_lin_scale, other_scale]
endif else begin
    guess=[wl_guess, gh_guess, other_guess]    
    scale=[wl_scale, gh_scale, other_scale]
endelse



;window,0,xsize=2050, ysize=900
if visualize eq 1 then window, 0, xsize=1000, ysize=400
!p.multi=[0,1,3]
;;;Run Minimization scheme

if model_func ne 0 then begin
    r=amoebafunction(guess)
endif else begin
    if min_type eq 'amoeba' then begin
        r=amoeba3(ftol, scale=scale, p0=guess,function_name='amoebafunction', $
                  function_value=fval, nmax=150000L)
        chi2=fval[0]
    endif else if min_type eq 'mpfit' then begin
        r=mpfit('amoebafunction', guess, bestnorm=chi2, ftol=ftol, maxiter=1000L, /quiet)
    endif
endelse

!p.multi=0

;stop clock
process_time=toc(clock)
endtime=systime()


;OUTPUT RESULTS

model_id=strtrim(nparam_wl,2)+"_"+strtrim(nparam_gh,2)+"_"+strtrim(nparam_gh_lin,2)+"_"+strtrim(nparam_other,2)

if n_elements(r) eq 1 then begin
    print, "Minimization scheme failed to converge"
    wl_result=-1
    gh_result=-1
    gh_lin_result=-1
    other_result=-1
    chi2=-1
    chi2perDOF=-1
endif else begin
    chi2perDOF=chi2/dof
    
    wl_result=dblarr(nparam_wl)
    gh_result=dblarr(nparam_gh)
    if lin_switch eq 1 then gh_lin_result=dblarr(nparam_gh_lin) 
    other_result=dblarr(nparam_other)

    print, "---------------------------------"
    print, "Model ID: ", model_id
    print, "---------------------------------"
    print, "Wavelength solution coeffiecients:"
    for i=0, nparam_wl-1 do begin
        wl_result[i]=r[i]
        print, r[i]
    endfor
    print, "---------------------------------"
    print, "Gauss-Hermite coefficients:"
    for i=nparam_wl, nparam_wl+nparam_gh-1 do begin
        gh_result[i-nparam_wl]=r[i]
        print, r[i]
    endfor
    print, "---------------------------------"

    if lin_switch eq 1 then begin
        print, "Gauss-Hermite linear trend coefficients: "
        for i=nparam_wl+nparam_gh, nparam_wl+(2*nparam_gh)-1 do begin
            gh_lin_result[i-nparam_wl-nparam_gh]=r[i] 
            print, r[i]
        endfor
        print, "---------------------------------"
    endif else gh_lin_result=!values.f_nan

    print, "Other coefficients:"
    for i=nparam_tot-nparam_other, nparam_tot-1 do begin
        other_result[i-nparam_tot+nparam_other]=r[i]
        print, r[i]
    endfor
    print, "---------------------------------"
    print, "Chi Squared, n_freeparams: "
    print, chi2
    print, nparam_tot
    
endelse

print, "---------------------------------"
print, "Processing took ", process_time, " seconds."


output_str={DATE:starttime, $
            ENDDATE:endtime, $
            WL_GUESS:wl_guess, $
            GH0_GUESS:gh_guess, $
            GH1_GUESS:gh_lin_guess, $
            OTHER_GUESS:other_guess, $
            WL_SCALE:wl_scale, $
            GH0_SCALE:gh_scale, $
            GH1_SCALE:gh_lin_scale, $
            OTHER_SCALE:other_scale, $
            WL_RESULT:wl_result, $
            GH0_RESULT:gh_result, $
            GH1_RESULT:gh_lin_result, $
            OTHER_RESULT:other_result, $
            OVERSAMP:oversamp, $
            FTOL:ftol, $
            NPIX_TRIM_START:npix_trim_start, $
            NPIX_TRIM_END:npix_trim_end, $
            NPARAM:nparam_tot, $
            DOF:DoF, $     
            VERSION:version_number, $
            TIME:process_time, $
            MIN_TYPE:min_type, $
            CHI2:chi2, $
            CHI2_PER_DOF:chi2perDOF $
           }


outpath_base='/home/stgilhool/RV_projects/IRCS_rv/cal_results/'
outpath_ext='Oct31/'
outpath=outpath_base+outpath_ext

outfile="model_"+model_id+".fits"

mwrfits, output_str, outpath+outfile
 
end






