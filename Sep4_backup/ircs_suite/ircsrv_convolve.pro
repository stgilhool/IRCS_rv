function ircs_convolve, v1, lsf, oversamp=oversamp, first_pix=first_pix, npix_select=npix_select
;take the LSF array (lsf, pixel) and convolves it with an oversampled
;spectrum.  Returns the oversampled output, padded with 1's at each
;end, if necessary (because of edge effects)

;Keywords
npix=1024L ;CAUTION this needs to be tuned for non-IRCS data sets
if n_elements(oversamp) eq 0 then oversamp=7L
if n_elements(first_pix) eq 0 then first_pix=0L
if n_elements(npix_select) eq 0 then npix_select=npix-first_pix

;Constant defs
npix_model=npix_select*oversamp
npix_lsf=oversamp*10L+1L

;Initialize some arrays for re-sorting of intensity vector
int_conv=dblarr(npix_model)
int_conv_mtx=dblarr(npix_lsf, npix_model)


for index=0,npix_model-1 do begin
    model_index=index+first_pix*oversamp
    if model_index lt (npix_lsf/2) or model_index gt (n_elements(v1)-1-(npix_lsf/2)) then $
      int_conv_mtx[*,index]=replicate(1d0, npix_lsf) $
    else int_conv_mtx[*,index]=v1[model_index-(npix_lsf/2):model_index+(npix_lsf/2)]
endfor
int_conv_temp=int_conv_mtx*lsf
int_conv=total(int_conv_temp, 1, /double)


return, int_conv

end
