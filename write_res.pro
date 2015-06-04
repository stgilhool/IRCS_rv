pro write_res, model_arr, obs_arr, wl_arr, err_arr

file='residuals_all.dat'


res=model_arr[*,3]-obs_arr[*,3]

wl=wl_arr[*,3]

err=err_arr[*,3]

openw, lun, file, /get_lun, /append

for i=0, n_elements(wl)-1 do begin
    
    printf, lun, wl[i], res[i], err[i], format='(D,D,D)'

endfor

close, lun
free_lun, lun

end

