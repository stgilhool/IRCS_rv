function cb_hermite, n, x_in, sigma

;n is order, x is the (x-x0) term in the Gaussian, and sigma is the
;sigma of the Gaussian, in the same units as x_in

n=n*1d0

x=x_in/sigma

h=0d

m=n/2l

for i=0l, m do h=h+(-1d0)^(i*1d0)*(2d0*x)^(n-2d0*i)/(factorial(m)*factorial(n-2*i))

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



function IRCS_cal_bestrun, filename

;deal with keywords

fits_info, filename, /silent, n_ext=last_ext

N_ext=last_ext

chi2s=dblarr(N_ext)
chi2pers=dblarr(N_ext)
for ext=0,N_ext-1 do begin
    run_info=mrdfits(filename, ext+1)
    chi2s[ext]=run_info.chi2
    chi2pers[ext]=run_info.chi2_per_dof
endfor

fail_ind=where(chi2s lt 0, nfails)

chi2per_finals=chi2pers
if nfails gt 0 then chi2per_finals[fail_ind]=999999d0

bestrun=where(chi2per_finals eq min(chi2per_finals))

;bestrun_str=mrdfits(filename, bestrun+1)
help, bestrun

return, bestrun[0]+1L 


end



pro IRCS_cal_interpret, filelist

readcol, filelist, filename, format='(A)'

nfiles=n_elements(filename)

ntestparams=11L

model_arr=dblarr(ntestparams, nfiles)

for file=0, nfiles-1 do begin
    
    bestrun=IRCS_cal_bestrun(filename[file])
    bestrun_info=mrdfits(filename[file], bestrun)
    
    nparam_wl=n_elements(bestrun_info.wl_result)
    nparam_gh0=n_elements(bestrun_info.gh0_result)
    nparam_other=n_elements(bestrun_info.other_result)
    nparam_tot=bestrun_info.nparam
    nparam_gh1=nparam_tot-nparam_wl-nparam_gh0-nparam_other
    
    model_arr[0, file]=bestrun_info.chi2
    model_arr[1, file]=bestrun_info.dof
    model_arr[2, file]=bestrun_info.chi2/bestrun_info.dof
    model_arr[3, file]=nparam_tot
    model_arr[4, file]=nparam_wl
    model_arr[5, file]=nparam_gh0
    model_arr[6, file]=nparam_gh1
    model_arr[7, file]=nparam_other
    model_arr[8, file]=bestrun

endfor

;index list of rows sorted by chi2/dof
nparam_order=sort(model_arr[3,*])

;sorts each column by nparam_order index list
for col=0, ntestparams-1 do begin
    column=model_arr[col,*]
    model_arr[col, *]=column[nparam_order]
endfor

for file=0, nfiles-1 do begin
    if file eq 0 then begin
        model_arr[9, file]=-999999
        model_arr[10,file]=-999999
    endif else begin
        numer=(model_arr[0,file-1]-model_arr[0,file])/(model_arr[1,file-1]-model_arr[1,file])
        denom=model_arr[2,file]
        f_val=numer/denom
        if f_val ge 0 then $
          prob=mpftest(f_val, model_arr[1,file-1]-model_arr[1,file], model_arr[1,file])
        if f_val lt 0 then begin
            fval=!values.f_nan
            prob=!values.f_nan
        endif
        model_arr[9,file]=f_val
        model_arr[10,file]=prob
    endelse
endfor

print, 'chi2', 'dof', 'chi2perdof', 'nparam', 'nparam_wl', 'nparam_gh0', 'nparam_gh1', 'nparam_other', 'extension', 'f_val', 'prob'
print, model_arr



end
