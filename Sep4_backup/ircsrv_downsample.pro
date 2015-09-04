function ircsrv_downsample, vec_over, npix_tophat
;Take an oversampled vector and downsample by convolving with a tophat
;function
;Oversampled vector should be offset so that 0th downsampled pixel is
;in the oversamp/2 element
;ie for oversamp=7
;Vec:xxx0xxxxxx1
;Ove:012345678910

vsize=size(vec_over)

n_dim=vsize[0]
dim_v=vsize[lindgen(n_dim)+1]

npix_over = dim_v[0]
npix_down = npix_over/npix_tophat

tophat=replicate(1d0/npix_tophat,npix_tophat)

avg_vec=convol(vec_over, tophat)

case n_dim of

    1: begin
        samp_index=(lindgen(npix_down)*npix_tophat)+(npix_tophat/2)
        vec_down = avg_vec[samp_index] 
    end
    
    2: begin
        samp_index=rebin((lindgen(npix_down)*npix_tophat)+(npix_tophat/2), npix_down, dim_v[1])
        index2=rebin(reform(lindgen(dim_v[1]), 1, dim_v[1]), npix_down, dim_v[1])
        vec_down = avg_vec[samp_index, index2]
    end
    
    3: begin
        samp_index=rebin((lindgen(npix_down)*npix_tophat)+(npix_tophat/2), [npix_down, dim_v[1:*]])
        index2=rebin(reform(lindgen(dim_v[1]), 1, dim_v[1],1), [npix_down, dim_v[1:*]])
        index3=rebin(reform(lindgen(dim_v[2]), 1, 1, dim_v[2]), [npix_down, dim_v[1:*]])
        vec_down = avg_vec[samp_index, index2, index3]
    end

endcase

return, vec_down

end
