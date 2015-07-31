function ircsrv_normcorr, k, first_pix=first_pix, npix_select=npix_select

n_norm = n_elements(k)

if n_norm gt 1 then begin
    norm_pts_y = k

    ;X vector of same dimensions as ircs-resolution spectrum
    norm_pts_xx=lindgen(npix_select)+first_pix

    ;Create and populate vector with x-coord of node pts
    norm_pts_x=dblarr(n_norm)

    if n_norm ge 2 then begin
        norm_pts_x[0]=first_pix
        norm_pts_x[n_norm-1]=first_pix+npix_select-1
        if n_norm ge 3 then begin
            x_increment=npix_select/(n_norm-2)
            for i=1, n_norm-2 do norm_pts_x[i]=norm_pts_x[i-1]+x_increment
        endif
    
    ;Interpolate between nodes with a spline
        if npix_select ge 128L and n_norm ge 5 then norm_factor=interpol(norm_pts_y, norm_pts_x, norm_pts_xx, /spline) $
        else norm_factor=interpol(norm_pts_y, norm_pts_x, norm_pts_xx)
    endif else if n_norm eq 1 then norm_factor=replicate(norm_pts_y, npix_select)


endif else begin
    
    print, "WARNING: No normalization"
    
    norm_factor = replicate(1d0, npix_select)

endelse

return, norm_factor

end
