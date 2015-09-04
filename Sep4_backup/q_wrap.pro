pro q_wrap, first_pix=first_pix, npix_select=npix_select, visualize=visualize, smooth_first=smooth_first

npix=1024L

if n_elements(first_pix) eq 0 then first_pix = 0L
if n_elements(npix_select) eq 0 then npix_select=128L < (npix-first_pix)
if n_elements(visualize) eq 0 then visualize=1
if n_elements(smooth_first) eq 0 then smooth_first = 0


current_tag1='Mar24_wlgh_test'
;This one had delta_rv.step = 1d-5, delta_wl.step = automatic (did not
;fit), sigma.step=1d-6

current_tag1='Mar24_wlgh_test2'
;delta_rv.step=1d-5, sigma.step=1d-5, delta_wl.relstep=1d-2, changed
;to sigma_start=[0.5, 1.0, 1.5], wl_start=[-1d-6, -1d-7, ..., 1d-6, 1d-7] 

current_tag1='Mar24_wlgh_test3'
;Same as test2, but with wl.relstep=1d-3, and sigma_start=0.75 only

current_tag1='Mar24_wlgh_test4'
;Same as test3, but with fixed rv


n_bases_lsf1=1L
  first_pix=310L
  npix_select=89L
  temp_num1=45L
  mode_1='mpfit'
  description='Many runs with different gh and wl guesses (1 parameter each)'


ircsrv_rvguesstest, first_pix=first_pix, npix_select=npix_select, visualize=1, run=0, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1


end