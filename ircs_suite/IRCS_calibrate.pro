function amoebafunction, p

common amoeba_info, run, niter, visualize, oversamp, npix, npix_trim_start, npix_trim_end, model_func, dof, wl_lab, int_lab, int_obs, nparam_wl, nparam_gh, nparam_other, lin_switch, trial_lsf, err, H_coeff, min_type, first_pix, npix_select


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
npix_lsf=(oversamp*10L)+1L
c=299792458.D



;Adjust the optical depth of the lab spectrum (I think this is okay to
;do for the high-resolution spectrum even though the wl scale isn't linear)
int_lab_depth=int_lab^tau_scale
int_obs_alias=int_obs

;;;Contstruct wl array and oversampled wl array from coeffs
x=dindgen(npix)
xx=(dindgen(npix*oversamp)-(oversamp/2))/oversamp

trial_soln=poly(x, wl_coeff)
samp_index=(lindgen(npix)*oversamp)+(oversamp/2)
;trial_soln_over=interpol(trial_soln, samp_index, xx)
trial_soln_over=poly(xx, wl_coeff)

;Interpolate the lab spectrum onto the oversampled trial grid
int_lab_over=interpol(int_lab_depth, wl_lab, trial_soln_over)

;Save full length vectors
trial_soln_full=trial_soln
trial_soln_over_full=trial_soln_over
int_lab_over_full=int_lab_over

;Redefine those vectors to reflect only the specified range
npix_model=npix_select*oversamp

x_select=dindgen(npix_select)+first_pix
xx_select=(dindgen(npix_model)-(oversamp/2))/oversamp+first_pix

int_obs_dim=int_obs_alias[first_pix:(first_pix+npix_select-1)]
err_dim=err[first_pix:(first_pix+npix_select-1)]
trial_soln=poly(x_select, wl_coeff)
samp_index=(lindgen(npix_select)*oversamp)+(oversamp/2)
trial_soln_over=poly(xx_select, wl_coeff)

;;;;
;Make LSF
;;;;

lsf_tic=tic()

 



;populate lsf coordinates
x_lsf=rebin((dindgen(npix_lsf)-(npix_lsf/2))/oversamp, npix_lsf, npix_model)
pixel_2d=rebin(reform((((dindgen(npix_model)-(oversamp/2))/oversamp)+first_pix), 1, npix_model), npix_lsf, npix_model)
;sig_arr=rebin(reform(sigma+(sigma_lin*dindgen(npix_model)), 1, npix_model), npix_lsf, npix_model)
sig_arr=sigma+(sigma_lin*pixel_2d)
x_sig_2d=x_lsf/sig_arr
x_sig_3d=rebin(x_sig_2d, npix_lsf, npix_model, nparam_gh)

;Make non-normalized gh-polynomial cube (npix_lsf x npix_model x
;basis)

gh_cube_nonorm=dblarr(npix_lsf, npix_model, nparam_gh)
gh_cube_norm=dblarr(npix_lsf, npix_model, nparam_gh)
gh_cube=dblarr(npix_lsf, npix_model, nparam_gh)
basis_coeff=dblarr(npix_lsf, npix_model, nparam_gh)
for basis=0, nparam_gh-1 do begin
    ;make gh_polynomial for given basis
    gh_cube_nonorm[*,*,basis]=poly(x_sig_2d, H_coeff[*,basis])*exp(-0.5*(x_sig_2d^2))
    ;normalize each basis
    norm=rebin(reform(total(abs(gh_cube_nonorm[*,*,basis]), 1, /double), 1, npix_model), npix_lsf, npix_model)
    gh_cube_norm[*,*,basis]=gh_cube_nonorm[*,*,basis]/norm
    ;Adjust amplitudes according to gh parameters
    if lin_switch eq 1 then basis_coeff[*,*,basis]=gh_coeff[basis]+(gh_lin_coeff[basis]*pixel_2d) $
      else if lin_switch eq 0 then basis_coeff[*,*,basis]=replicate(gh_coeff[basis], npix_lsf, npix_model)
    gh_cube[*,*,basis]=basis_coeff[*,*,basis]*gh_cube_norm[*,*,basis]
endfor

;Sum the bases
if nparam_gh gt 1 then trial_lsf_nonorm=total(gh_cube, 3, /double) $
  else trial_lsf_nonorm=gh_cube[*,*,0]

;account for unrealistically negative lsfs

min_vec=min(trial_lsf_nonorm, dimension=1)
neg_ind=where(min_vec lt 0, neg_count)
if neg_count gt 0 then begin
    neg_vec=dblarr(npix_model)
    neg_vec[neg_ind]=min_vec[neg_ind]

    neg_arr=rebin(reform(neg_vec, 1, npix_model), npix_lsf, npix_model)
;TAKE OUT LATER
print, min(neg_vec)
    trial_lsf_nonorm_copy=trial_lsf_nonorm
    trial_lsf_nonorm=trial_lsf_nonorm-neg_arr
endif

;normalize lsf
lsf_norm=rebin(reform(total(trial_lsf_nonorm, 1, /double), 1, npix_model), npix_lsf, npix_model)
trial_lsf=trial_lsf_nonorm/lsf_norm
if lin_switch eq 0 then trial_lsf=trial_lsf[*,0]


;;;;
;Fix LSF center of mass
;;;;
lsf_com=rebin(reform(total(trial_lsf*x_sig_2d, 1, /double), 1, npix_model), npix_lsf, npix_model)
x_sig_2d_com=x_sig_2d+lsf_com
x_sig_3d_com=rebin(x_sig_2d_com, npix_lsf, npix_model, nparam_gh)

;;;Remake LSF
;Make non-normalized gh-polynomial cube (npix_lsf x npix_model x
;basis)

gh_cube_nonorm_com=dblarr(npix_lsf, npix_model, nparam_gh)
gh_cube_norm_com=dblarr(npix_lsf, npix_model, nparam_gh)
gh_cube_com=dblarr(npix_lsf, npix_model, nparam_gh)
for basis=0, nparam_gh-1 do begin
    ;make gh_polynomial for given basis
    gh_cube_nonorm_com[*,*,basis]=poly(x_sig_2d_com, H_coeff[*,basis])*exp(-0.5*(x_sig_2d_com^2))
    ;normalize each basis
    norm_com=rebin(reform(total(abs(gh_cube_nonorm_com[*,*,basis]), 1, /double), 1, npix_model), npix_lsf, npix_model)
    gh_cube_norm_com[*,*,basis]=gh_cube_nonorm_com[*,*,basis]/norm_com
    
    gh_cube_com[*,*,basis]=basis_coeff[*,*,basis]*gh_cube_norm_com[*,*,basis]
endfor

;Sum the bases
if nparam_gh gt 1 then trial_lsf_nonorm_com=total(gh_cube_com, 3, /double) $
  else trial_lsf_nonorm_com=gh_cube_com[*,*,0]

;account for unrealistically negative lsfs

min_vec_com=min(trial_lsf_nonorm_com, dimension=1)
neg_ind_com=where(min_vec_com lt 0, neg_count_com)
if neg_count_com gt 0 then begin
    neg_vec_com=dblarr(npix_model)
    neg_vec_com[neg_ind_com]=min_vec_com[neg_ind_com]

    neg_arr_com=rebin(reform(neg_vec_com, 1, npix_model), npix_lsf, npix_model)

    trial_lsf_nonorm_com_copy=trial_lsf_nonorm_com
    trial_lsf_nonorm_com=trial_lsf_nonorm_com-neg_arr_com
endif

;normalize lsf
lsf_norm_com=rebin(reform(total(trial_lsf_nonorm_com, 1, /double), 1, npix_model), npix_lsf, npix_model)
trial_lsf_com=trial_lsf_nonorm_com/lsf_norm_com
if lin_switch eq 0 then trial_lsf_com=trial_lsf_com[*,0]

;make copy of old lsf 
trial_lsf_copy=trial_lsf
trial_lsf=trial_lsf_com

;lsf_com_com=rebin(reform(total(trial_lsf*x_sig_2d, 1, /double), 1, npix_model), npix_lsf, npix_model)

lsf_t=toc(lsf_tic)
print, "LSF construction took : ", lsf_t, " seconds."


;;;I THINK THIS IS WRONG
;decompose the lsf if model_func is set
if model_func ne 0 then begin
    lsf_norm_3d=rebin(lsf_norm, npix_lsf, npix_model, nparam_gh)
    lsf_norm_3d_com=rebin(lsf_norm_com, npix_lsf, npix_model, nparam_gh)

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

int_lab_conv=dblarr(npix_model)
int_lab_conv_mtx=dblarr(npix_lsf, npix_model)



for index=0,npix_model-1 do begin
    model_index=index+first_pix*oversamp
    if model_index lt (npix_lsf/2) or model_index gt (n_elements(int_lab_over_full)-1-(npix_lsf/2)) then $
      int_lab_conv_mtx[*,index]=replicate(1d0, npix_lsf) $
    else int_lab_conv_mtx[*,index]=int_lab_over_full[model_index-(npix_lsf/2):model_index+(npix_lsf/2)]
endfor
int_lab_conv_temp=int_lab_conv_mtx*trial_lsf
int_lab_conv=total(int_lab_conv_temp, 1, /double)


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
    x_x=lindgen(npix_select)+first_pix
    x_norm=dblarr(nparam_norm)
    x_norm[0]=first_pix
    x_norm[nparam_norm-1]=first_pix+npix_select-1
    x_increment=npix_select/(nparam_norm-2)
    for i=1, nparam_norm-2 do x_norm[i]=x_norm[i-1]+x_increment

    norm_factor=interpol(norm_pts, x_norm, x_x, /spline)
    
    int_model=int_model*norm_factor
endif
;;;


;;;;
;;;Get chi^2
;;;;

chi_tic=tic()	;Start clock



res=(int_obs_dim-int_model)
res2=res^2
dev=res/err_dim
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
    x1=npix_select/6
    x2=npix_select/2
    x3=(5*npix_select/6)
    
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
;    if x1 gt npix_lsf/2 then begin
    x1_vec=(dindgen(npix_lsf)-(npix_lsf/2))+x1+first_pix
    x2_vec=(dindgen(npix_lsf)-(npix_lsf/2))+x2+first_pix
    x3_vec=(dindgen(npix_lsf)-(npix_lsf/2))+x3+first_pix
;    endif else begin
        


;plot_tic=tic()
    
    title_str="NH3 Model Fitting Process | run "+strtrim(run+1,2)+"/"+strtrim(niter,2)+" | chi2/DoF : "+strtrim(chi2perdegree,2)

    plot, trial_soln, int_obs_dim,yr=[0.1,1.1], /xs, title=title_str,xtitle="Wavelength (microns)", ytitle="Relative Instensity", charsize=2.0 
    oplot, trial_soln, int_model, color=200 

    plot, trial_soln, res, title="Residuals | RMS : ", xtitle="Wavelength (microns)", yr=[-0.1, 0.1], ps=3, /xs

    plot, x1_vec, lsf1, yr=[-0.05*max_lsf,1.1*max_lsf], xr=[first_pix,first_pix+npix_select-1], /xs
    oplot, x2_vec, lsf2
    oplot, x3_vec, lsf3
;stop
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
    ;IRCS_make_figs, chi2perdegree, trial_soln, trial_lsf, oversamp, npix_lsf, int_obs, int_model, npix_trim_start, npix_trim_end, lin_switch
    ;stop
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
pro IRCS_calibrate, wl_guess, wl_scale, gh_guess, gh_scale, gh_lin_guess, gh_lin_scale, other_guess, other_scale, run, niter, lin_switch=lin_switch, min_type=min_type, visualize=visualize, oversamp=oversamp, npix_trim_start=npix_trim_start, npix_trim_end=npix_trim_end, model_func=model_func, first_pix=first_pix, npix_select=npix_select, object=object, epoch=epoch, trace=trace


if n_elements(lin_switch) eq 0 then lin_switch=1
if n_elements(min_type) eq 0 then min_type="mpfit"
if n_elements(visualize) eq 0 then visualize=0
if n_elements(oversamp) eq 0 then oversamp=7L
if n_elements(npix_trim_start) eq 0 then npix_trim_start=15L
if n_elements(npix_trim_end) eq 0 then npix_trim_end=15L
if n_elements(model_func) eq 0 then model_func=0
if n_elements(first_pix) eq 0 then first_pix=0L
if n_elements(npix_select) eq 0 then npix_select=1024L-first_pix 
if n_elements(object) eq 0 then object='GJ273'
if n_elements(epoch) eq 0 then epoch='18Jan2011'
if n_elements(trace) eq 0 then trace='AB1'


common amoeba_info, runn, n_iter, visual, oversample, npix, npixtrim_start, npixtrim_end, modelfunc, dof, wl_lab, int_lab, int_obs, nparam_wl, nparam_gh, nparam_other, linear_switch, trial_lsf, err, H_coeff, mintype, firstpix, npixselect

;Define constants and stuff
version_number='1.2'

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
firstpix=first_pix
npixselect=npix_select

starttime=systime()

H_coeff=sg_hermite_coeff()

;start clock
clock=tic()

;read in both spectra
rootpath='/home/stgilhool/RV_projects/IRCS_rv/data/'
modelpath=rootpath+'supplemental/'
epochpath=rootpath+'epoch/'+epoch+'/'
outpath=epochpath+'calib_results/'
obspath=epochpath+'final_spectra/'

modelfile= modelpath+'NH3_model.dat'
obsfile=obspath+'NH3_'+strjoin([object, epoch, trace], '_')+'.fits'

readcol, modelfile, wl_lab, int_lab, format='D,D'
obs_fits=mrdfits(obsfile, 1)
int_obs_nonorm_backwards=obs_fits.nh3_spectrum
err_backwards=obs_fits.nh3_sigma
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
dof=npix_select-nparam_tot-(npix_trim_start+npix_trim_end)

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
        r=mpfit('amoebafunction', guess, bestnorm=chi2, ftol=ftol, /quiet)
    endif
endelse

!p.multi=0

;stop clock
process_time=toc(clock)
endtime=systime()


;OUTPUT RESULTS

model_id=strjoin([object, epoch, trace, strtrim(nparam_wl,2),strtrim(nparam_gh,2),strtrim(nparam_gh_lin,2),strtrim(nparam_other,2)],'_')

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



outfile=model_id+".fits"

outfinal=outpath+outfile

if file_test(outfinal) then begin
    fits_info, outfinal, n_ext=prev_ext, /silent
    extension=prev_ext+1
endif else extension=1


output_str={OBJ:object, $
            EPOCH:epoch, $
            TRACE:trace, $
            STARTTIME:starttime, $
            ENDTIME:endtime, $
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
            FIRST_PIX:first_pix, $
            NPIX_SELECT:npix_select, $
            NPARAM:nparam_tot, $
            DOF:DoF, $     
            VERSION:version_number, $
            TIME:process_time, $
            MIN_TYPE:min_type, $
            EXT:extension, $
            CHI2:chi2, $
            CHI2_PER_DOF:chi2perDOF $
           }




mwrfits, output_str, outfinal
 
end






