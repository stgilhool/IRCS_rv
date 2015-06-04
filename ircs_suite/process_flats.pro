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

;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;

pro process_flats, gc_on_list, gc_off_list, nc_on_list, nc_off_list

;set file names
filepath='/home/stgilhool/RV_projects/IRCS_rv/IRCS_flats/'
maskname='mask1D.fits'
gc_on_list='gc_on_list.txt'
gc_off_list='gc_off_list.txt'
nc_on_list='nc_on_list.txt'
nc_off_list='nc_off_list.txt'

gc_on_list=filepath+gc_on_list
gc_off_list=filepath+gc_off_list
nc_on_list=filepath+nc_on_list
nc_off_list=filepath+nc_off_list
maskfile=filepath+maskname

;read in mask
mask=mrdfits(maskfile,0,/dscale)

;read in file names
readcol, gc_on_list, gc_on_file, format='A'
readcol, gc_off_list, gc_off_file, format='A'
readcol, nc_on_list, nc_on_file, format='A'
readcol, nc_off_list, nc_off_file, format='A'
;add filepath
gc_on_file=filepath+gc_on_file
gc_off_file=filepath+gc_off_file
nc_on_file=filepath+nc_on_file
nc_off_file=filepath+nc_off_file

;initialize pointer arrays to store spectra and error data
nflats=n_elements(gc_on_file)
temp=mrdfits(gc_on_file[0], 0)
ncols=n_elements(temp[0,*])
nrows=n_elements(temp[*,0])

gc_on_arr=dblarr(nflats, nrows, ncols)
gc_off_arr=dblarr(nflats, nrows, ncols)
nc_on_arr=dblarr(nflats, nrows, ncols)
nc_off_arr=dblarr(nflats, nrows, ncols)

gc_on_var=dblarr(nflats, nrows, ncols)
gc_off_var=dblarr(nflats, nrows, ncols)
nc_on_var=dblarr(nflats, nrows, ncols)
nc_off_var=dblarr(nflats, nrows, ncols)



;Error related info
gain=3.8d0 ;electrons/ADU
read_noise=68d0 ;electrons/NDR
dark_current=0.05d0 ;electrons/second
time=0.5 ;seconds

;Read in spectra in number of photo-electrons
for file=0,nflats-1 do begin
    gc_on_arr[file, *, *]=mrdfits(gc_on_file[file], 0, /dscale)*gain
    gc_off_arr[file, *, *]=mrdfits(gc_off_file[file], 0, /dscale)*gain
    nc_on_arr[file, *, *]=mrdfits(nc_on_file[file], 0, /dscale)*gain
    nc_off_arr[file, *, *]=mrdfits(nc_off_file[file], 0, /dscale)*gain
                               
endfor

;variance at all pixels in each pixel in each flat
gc_on_var=gc_on_arr + (read_noise^2) + (dark_current*time)
gc_off_var=gc_off_arr + (read_noise^2) + (dark_current*time)
nc_on_var=nc_on_arr + (read_noise^2) + (dark_current*time)
nc_off_var=nc_off_arr + (read_noise^2) + (dark_current*time)

;sum along columns of pixels and take median of the flats to get
;single 1D spectra for gc/nc , on/off
;gc_on_med=median(total(gc_on_arr, 3), dimension=1, /double)
;gc_off_med=median(total(gc_off_arr, 3), dimension=1, /double)
;nc_on_med=median(total(nc_on_arr, 3), dimension=1, /double)
;nc_off_med=median(total(nc_off_arr, 3), dimension=1, /double)

;get the mean in a similar way
;gc_on_mean=mean(total(gc_on_arr, 3), dimension=1, /double)
;gc_off_mean=mean(total(gc_off_arr, 3), dimension=1, /double)
;nc_on_mean=mean(total(nc_on_arr, 3), dimension=1, /double)
;nc_off_mean=mean(total(nc_off_arr, 3), dimension=1, /double)

;Sum the columns and spectra to get total counts per pixel
gc_on_tot=total(total(gc_on_arr, 3), 1)
gc_off_tot=total(total(gc_off_arr, 3), 1)
nc_on_tot=total(total(nc_on_arr, 3), 1)
nc_off_tot=total(total(nc_off_arr, 3), 1)

gc_on_tot_var=total(total(gc_on_var, 3),1)
gc_off_tot_var=total(total(gc_off_var, 3),1)
nc_on_tot_var=total(total(nc_on_var, 3),1)
nc_off_tot_var=total(total(nc_off_var, 3),1)
;SNR is reduced by the following factor when taking the median
;median_factor=sqrt(!PI/(2*nflats))

;calculate noise in the various spectra
;gc_on_err=median(sqrt(total(gc_on_var, 3)), dimension=1, /double)*median_factor
;gc_off_err=median(sqrt(total(gc_off_var, 3)), dimension=1, /double)*median_factor
;nc_on_err=median(sqrt(total(nc_on_var, 3)), dimension=1, /double)*median_factor
;nc_off_err=median(sqrt(total(nc_off_var, 3)), dimension=1, /double)*median_factor


;Now take on-off
gc_arr=gc_on_tot-gc_off_tot
nc_arr=nc_on_tot-nc_off_tot
;errors
;gc_err=sqrt((gc_on_err^2)+(gc_off_err^2))
;nc_err=sqrt((nc_on_err^2)+(nc_off_err^2))
gc_arr_var=gc_on_tot_var+gc_off_tot_var
nc_arr_var=nc_on_tot_var+nc_off_tot_var

;Finally, do gc/nc
observed_spectrum=gc_arr/nc_arr
;errors
var_final=(observed_spectrum^2)*((gc_arr_var/(gc_arr^2))+(nc_arr_var/(nc_arr^2)))
error=sqrt(var_final)

out_str={nh3_spec:observed_spectrum,nh3_err:error,nh3_mask:mask}
outfile='NH3_obs.fits'
mwrfits, out_str, outfile, /create

stop

end




;gc_off_arr=ptrarr(nflats, /allocate_heap)
;nc_on_arr=ptrarr(nflats, /allocate_heap)
;nc_off_arr=ptrarr(nflats, /allocate_heap)

;gc_on_err=ptrarr(nflats, /allocate_heap)
;gc_off_err=ptrarr(nflats, /allocate_heap)
;nc_on_err=ptrarr(nflats, /allocate_heap)
;nc_off_err=ptrarr(nflats, /allocate_heap)

;     *gc_on_arr[file]=mrdfits(gc_on_file[file], 0, /dscale)
;     *gc_off_arr[file]=mrdfits(gc_off_file[file], 0, /dscale)
;     *nc_on_arr[file]=mrdfits(nc_on_file[file], 0, /dscale)
;     *nc_off_arr[file]=mrdfits(nc_off_file[file], 0, /dscale)

;     *gc_on_err[file]=sqrt(total(*gc_on_arr[file],2)+$
;                           (ncols*(read_noise^2))+$
;                           (time*dark_current*ncols))
;     *gc_off_err[file]=sqrt(total(*gc_off_arr[file],2)+$
;                            (ncols*(read_noise^2))+$
;                            (time*dark_current*ncols))
;     *nc_on_err[file]=sqrt(total(*nc_on_arr[file],2)+$
;                           (ncols*(read_noise^2))+$
;                           (time*dark_current*ncols))
;     *nc_off_err[file]=sqrt(total(*nc_off_arr[file],2)+$
;                            (ncols*(read_noise^2))+$
;                            (time*dark_current*ncols))
