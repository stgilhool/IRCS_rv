function continuum_fit, x_vec, int, low_rej=low_rej, high_rej=high_rej, PIX_MASK=mask, display=display

;Initialize arrays and loop variables
npix=n_elements(x_vec)
if n_elements(mask) eq 0 then mask=replicate(0,npix)
if n_elements(low_rej) eq 0 then low_rej=0.5
if n_elements(high_rej) eq 0 then high_rej=3.
if n_elements(display) eq 0 then display=0
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
;         if display ne 0 then print, "Too many rejections. Exiting loop at iteration number: ", iteration
;         break
;         endif

    ;Calculate RMSE
    rmse=sqrt(mean((fit_iter-int_iter)^2))

    ;Display fit if display is set
    if display ne 0 then begin
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
