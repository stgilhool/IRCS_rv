function downsample_tophat, vec_over, npix_tophat
;Take an oversampled vector and downsample by convolving with a tophat
;function
;Oversampled vector should be offset so that 0th downsampled pixel is
;in the oversamp/2 element
;ie for oversamp=7
;Vec:xxx0xxxxxx1
;Ove:012345678910

npix_over = n_elements(vec_over)
npix_down = npix_over/npix_tophat

tophat=replicate(1d0/npix_tophat,npix_tophat)

avg_vec=convol(vec_over, tophat)

samp_index=(lindgen(npix_down)*npix_tophat)+(npix_tophat/2)

vec_down = avg_vec[samp_index] 

return, vec_down

end
