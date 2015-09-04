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


function cb_hermite, n, x_in, sigma

;n is order, x is the (x-x0) term in the Gaussian, and sigma is the
;sigma of the Gaussian, in the same units as x_in

n=n*1d0

x=x_in/sigma

h=0d

m=n/2l

for i=0l, m do h=h+(-1d0)^(i*1d0)*(2d0*x)^(n-2d0*i)/(factorial(i)*factorial(n-2*i))

output=(1./sqrt(sigma))*(2d0^n*factorial(n)*Sqrt(!pi))^(-0.5d0)*factorial(n)*h*exp(-x_in^2d0/(2d0*sigma^2d0))



return, output

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



pro IRCS_cal_display, nparam_wl, nparam_gh, nparam_gh_lin, nparam_other, extension=extension, short=short, pix0=pix0

;deal with keywords
nwl=strtrim(nparam_wl,2)
ngh0=strtrim(nparam_gh,2)
ngh1=strtrim(nparam_gh_lin,2)
nother=strtrim(nparam_other,2)
resultpath='/home/stgilhool/RV_projects/IRCS_rv/cal_results/Oct31/'
outputpath='/home/stgilhool/RV_projects/IRCS_rv/cal_figures/'


fileparams1=nwl+'_'+ngh0+'_'+ngh1+'_'+nother



if n_elements(short) eq 0 then begin
    short=0
    fileroot='model_' 
    endif else if short eq 1 then fileroot='model_short_'

resultfile=resultpath+fileroot+fileparams1+'.fits'
fits_info, resultfile, /silent, n_ext=last_ext

if n_elements(extension) eq 0 then extension=last_ext else begin
    if extension gt last_ext then message, "Not enough extensions in fits file"
endelse


fileparams=nwl+'_'+ngh0+'_'+ngh1+'_'+nother+'_'+strtrim(long(extension),2)





if nparam_gh eq nparam_gh_lin then lin_switch=1 else lin_switch=0

;read in both spectra
readcol, 'NH3_model.dat', wl_lab, int_lab, format='D,D'
;readcol, 'NH3_obs_nonorm.dat', trash, int_obs_nonorm_backwards, format='I,D'
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


;Read in results
a=mrdfits(resultfile, extension)
wl_res=a.wl_result
gh0_res=a.gh0_result
gh1_res=a.gh1_result
other_res=a.other_result

tau_scale=other_res[0]
sigma=gh0_res[0]
wl_coeff=wl_res
gh_coeff=[1,gh0_res[1d0:-1]]

H_coeff=sg_hermite_coeff()

if lin_switch eq 1 then begin
    sigma_lin=gh1_res[0]
    gh_lin_coeff=[0d0, gh1_res[1:-1]]
endif

if short eq 0 then npix=1024L else if short eq 1 then npix=512L
x=dindgen(npix)
    
;Use results
npix_lsf=71L
oversamp=7
c=299792458.D
npix_trim=50L
dof=a.dof
nparam=a.nparam



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




lsf_tic=tic()
;Make LSF
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
neg_vec=dblarr(npix_model)

if neg_count gt 0 then begin
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

lsf_norm_3d=rebin(lsf_norm, npix_lsf, npix_model, nparam_gh)

;decompose the lsf
if neg_count gt 0 then begin
    neg_arr_3d=rebin(neg_arr, npix_lsf, npix_model, nparam_gh)
    lsf_decomp=(gh_cube-neg_arr_3d)/lsf_norm_3d
endif else begin
    lsf_decomp=gh_cube/lsf_norm_3d
endelse



; !p.multi=0

; test1=(indgen(5)*(npix_model-1))/5
; for i=0, 4 do begin
;     test=test1[i]
;     num=string(test)
;     plot,trial_lsf[*,test], title=num
;     wait,0.1
; endfor
; stop

; !p.multi=[0,1,3]


;Convolve model with lsf
convol_tic=tic()
if lin_switch eq 0 then begin
    int_lab_conv=convol(int_lab_over, trial_lsf) ;, /norm) 
endif else begin
    int_lab_conv=dblarr(npix_model)
    ;int_lab_conv[0:npix_lsf/2-1]=1d0
    ;int_lab_conv[npix_model-1-(npix_lsf/2):npix_model-1]=1d0
                                ;Make matrix of pieces of the spectrum
                                ;that will be convolved with lsf at
                                ;each pixel
    ;int_lab_conv_mtx=dblarr(npix_model, npix_lsf)
    
    ;for index=0,npix_model-1-npix_lsf do begin
    ;    int_lab_conv_mtx[index, *]=int_lab_over[index:index+npix_lsf-1]
    ;endfor
    int_lab_conv_mtx=dblarr(npix_lsf, npix_model)
    
    for index=0,npix_model-1 do begin
        if (index lt npix_lsf/2 or index gt (npix_model-1-(npix_lsf/2))) then $
            int_lab_conv_mtx[*,index]=replicate(1d0, npix_lsf) $
          else int_lab_conv_mtx[*,index]=int_lab_over[index-(npix_lsf/2):index+(npix_lsf/2)]
    endfor
    int_lab_conv_temp=int_lab_conv_mtx*trial_lsf
    int_lab_conv=total(int_lab_conv_temp, 1, /double)
;help, int_lab_conv
;stop
endelse
        
convol_t=toc(convol_tic)
;print, "Convolution took : ", convol_t, " seconds."

;    for x_var=0,npix_model-1 do begin
;        x_coord_lsf=x_lsf+x_var ;shift the position of the lsf on the spectrum
;        lsf_trim_index=where(x_coord_lsf ge 0 and x_coord_lsf lt npix_model)
;        x_lsf_mod=x_coord_lsf[lsf_trim_index]
;        trial_lsf_mod=trial_lsf[lsf_trim_index,x_var]
;        int_overlap=int_lab_over[x_coord_lsf]
;        int_lab_conv[x_var]=total(int_overlap*trial_lsf_mod)
;stop
;stop

;Down-sample the model to IRCS resolution
npix_tophat=oversamp
tophat=replicate(1d0/npix_tophat,npix_tophat)
int_lab_avg=convol(int_lab_conv, tophat)
int_model=int_lab_avg[samp_index] 


;Correct subtle variations in continuum
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

;;;Get chi^2
chi_tic=tic()	;Start clock


res=(int_obs-int_model)
res2=res^2
dev=res/err
dev2=dev^2
chi1_vec=dev[npix_trim:-1*(npix_trim+1)]
chi2_vec=dev2[npix_trim:-1*(npix_trim+1)]
chi2=total(chi2_vec, /double, /nan)
chi2perdegree=chi2/dof

chi_t=toc(chi_tic);Stop clock

dof=npix-nparam-140L
chi2_dof=chi2/dof

;;;Do animation of fitting spectrum for full spectrum
; !p.multi=0
;Set up some stuff to help display lsf across spectrum
x1=100L
x2=500L
x3=900L
;x2=250L
;x3=450L

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

;get rms
residuals=int_obs-int_model
avg_residual=mean(residuals[49:-50])
mse=0d0
n=0d0
for pix=49, npix-50 do begin
    ;mse=mse+(residuals[pixel]-avg_residual)^2
    mse=mse+(residuals[pix]^2)
    n++
endfor
rmse=sqrt(mse/n)

chi2_dof=a.chi2_per_dof ;this is just thrown in

;title_str='n_wl: '+strtrim(nparam_wl, 2)+ $
;  ' | n_gh0: '+strtrim(nparam_gh, 2)+ $
;  ' | n_gh01 '+strtrim(nparam_gh_lin, 2)+ $
;  ' | other: '+strtrim(nparam_other, 2)+ $
;  ' | chi2: '+strtrim(chi2, 2)
title_str="NH3 Model Fitting Process | chi2/degree of freedom : "+strtrim(chi2_dof,2)
phi=findgen(32)*(!PI*2/32.)
phi = [ phi, phi(0) ]
usersym, cos(phi), sin(phi), /fill

;Start plotting
if n_elements(pix0) eq 0 then lsf_cent_lowres=x3 else lsf_cent_lowres=pix0
lsf_cent=lsf_cent_lowres*oversamp+(npix_lsf/2)

set_plot, 'ps'
device, /encapsulated, filename = outputpath+'lsfdecomp_'+fileparams+'.eps'
!p.font=0
device, /color, bits=8
loadct, 12


!p.multi=[0,2,6]

for base=0, nparam_gh-1 do begin
    plot, lsf_decomp[*,lsf_cent,base], /xs, color=18L*base, title=strtrim(base,2)
endfor

device, /close




set_plot, 'ps'
device, /encapsulated, filename =  outputpath+'lsfdecomp_oplot_'+fileparams+'.eps'
!p.font=0
device, /color, bits=8
loadct, 12


!p.multi=0


plot, lsf_decomp[*,lsf_cent,0], /xs, yr=[-0.01, 0.04], title="all bases"

for i=1, nparam_gh-1 do begin
    oplot, lsf_decomp[*,lsf_cent,i], color=18L*i
    oplot, [i,i], [!y.crange[0],!y.crange[1]], color=18L*i
endfor




device, /close




set_plot, 'ps'
device, /encapsulated, filename =  outputpath+'lsfdecomp_cumulative_'+fileparams+'.eps'
!p.font=0
device, /color, bits=8
loadct, 12


!p.multi=[0,2,6]

for base=0, nparam_gh-1 do begin
    lsf_cumulative=dblarr(npix_lsf)
    for maxbase=0, base do begin
        lsf_cumulative=lsf_cumulative+lsf_decomp[*,lsf_cent,maxbase]
    endfor
    plot, lsf_cumulative, /xs, color=18L*base, title="Sum of bases up to order "+ strtrim(base,2)
endfor

device, /close




; set_plot, 'x'

; !p.multi=[0,2,6]

; for base=0, nparam_gh-1 do begin
;     lsf_cumulative_even=dblarr(npix_lsf)
;     lsf_cumulative_odd=dblarr(npix_lsf)
;     for maxbase=0, base do begin
;         if (long(maxbase) and 1) eq 0 then lsf_cumulative_even=lsf_cumulative_even+lsf_decomp[*,lsf_cent,maxbase]
;         if (long(maxbase) and 1) eq 1 then lsf_cumulative_odd=lsf_cumulative_odd+lsf_decomp[*,lsf_cent,maxbase]
;     endfor
;     if (long(base) and 1) eq 0 then plot, lsf_cumulative_even, /xs, title="Even bases up to order "+ strtrim(base,2)
;     if (long(base) and 1) eq 1 then plot, lsf_cumulative_odd, /xs, title="Odd bases up to order "+ strtrim(base,2)
; endfor

; stop


; set_plot, 'ps'
; device, /encapsulated, filename =  outputpath+'lsfdecomp_evenodd_'+fileparams+'.eps'
; !p.font=0
; device, /color, bits=8
; loadct, 12
; !p.multi=[0,2,6]

; for base=0, nparam_gh-1 do begin
;     lsf_cumulative_even=dblarr(npix_lsf)
;     lsf_cumulative_odd=dblarr(npix_lsf)
;     for maxbase=0, base do begin
;         if (long(maxbase) and 1) eq 0 then lsf_cumulative_even=lsf_cumulative_even+lsf_decomp[*,lsf_cent,maxbase]
;         if (long(maxbase) and 1) eq 1 then lsf_cumulative_odd=lsf_cumulative_odd+lsf_decomp[*,lsf_cent,maxbase]
;     endfor
;     if (long(base) and 1) eq 0 then plot, lsf_cumulative_even, /xs, title="Even bases up to order "+ strtrim(base,2)
;     if (long(base) and 1) eq 1 then plot, lsf_cumulative_odd, /xs, title="Odd bases up to order "+ strtrim(base,2)
; endfor

; device, /close

; stop



set_plot, 'ps'
device, /encapsulated, filename =  outputpath+'fit_'+fileparams+'.eps'
!p.font=0
device, /color, bits=8
loadct, 12

;window, 0, xsize=900, ysize=600
!p.multi=[0,1,2]


plot, trial_soln, int_obs, /xs, title=title_str,xtitle="Wavelength (microns)", thick=5.0, xr=[2.31,2.33], yr=[0.4,1.1]
;oplot, trial_soln, int_model, color=200, linestyle=0
;oplot, trial_soln[dindgen(npix/3)*3], int_model[dindgen(npix/3)*3], ps=8, symsize=0.6, color=200 
oplot, trial_soln, int_model, ps=8, symsize=0.4, color=200 
oplot, trial_soln, residuals+0.6, ps=8, symsize=0.2
al_legend, ['NH3 observed', 'NH3 model', 'Residuals'], psym=[-0,88,88], color=['black','red','black']
xyouts,2.32, 0.56, 'RMS : '+strtrim(rmse,2)
;axgap, [0,0.7]
;, linsize=[0.5,1,1];, number=[1,1,5]
;yr=[-0.05, 0.05], ps=3, charsize=1.5, /xs
;plot, x1_vec, lsf1, yr=[-0.05*max_lsf,1.1*max_lsf], xr=[0,npix], /xs
;oplot, x2_vec, lsf2
;oplot, x3_vec, lsf3
;plot, x_x, norm_factor, yr=[0.5, 1.5], /xs

device, /close


!p.multi=0
device, /encapsulated, filename =  outputpath+'lsf_'+fileparams+'.eps'
!p.font=0
device, /color, bits=8
loadct, 12
;window,1,xsize=700, ysize=600

if lin_switch eq 1 then begin
    plot, dindgen(npix_lsf)/oversamp, lsf1, xtitle='Pixel', title="Line Spread Function", /xs
    oplot, dindgen(npix_lsf)/oversamp, lsf2, linestyle=2
    oplot, dindgen(npix_lsf)/oversamp, lsf3, linestyle=4
    al_legend, ['LSF @ pixel 100', 'LSF @ pixel 500', 'LSF @ pixel 900'], linestyle=[0,2,4]
endif else begin
    plot, dindgen(npix_lsf), trial_lsf
endelse
device, /close



;OUTPUT RESIDUALS - DISABLED


; ;set_plot, 'x'
; incr=5
; skip=dindgen(npix/incr)*incr

; ;window, 0, xsize=1500, ysize=300
; device, /encapsulated, filename = 'residuals_sparse_'+fileparams+'.eps'
; !p.font=0
; device, /color, bits=8
; loadct, 12
; plot, x[skip],residuals[skip], ps=8, /xs, yr=[-0.05, 0.05], title="Residuals (every 5th pixel"
; errplot, x[skip], residuals[skip]+err[skip], residuals[skip]-err[skip]
; device, /close

; device, /encapsulated, filename = 'residuals_'+fileparams+'.eps'
; !p.font=0
; device, /color, bits=8
; loadct, 12
; plot, x,residuals, ps=3, /xs, yr=[-0.05, 0.05], title="Residuals"
; errplot, x, residuals+err, residuals-err
; device, /close

; device, /encapsulated, filename = 'residuals_zoom_'+fileparams+'.eps'
; !p.font=0
; device, /color, bits=8
; loadct, 12
; plot, trial_soln,residuals, ps=3, /xs, yr=[-0.05, 0.05], title="Residuals", xr=[2.31,2.33]
; errplot, trial_soln, residuals+err, residuals-err
; device, /close


; set_plot, 'x'

; openw, lun, 'residuals.txt', /get_lun
; for i=0,npix-1 do begin
;     printf, lun, residuals[i], x[i], trial_soln[i]
; endfor
; close, lun
; free_lun, lun


; stop


end
