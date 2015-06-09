pro q_wrap3, first_pix=first_pix, npix_select=npix_select, visualize=visualize, smooth_first=smooth_first

npix=1024L

if n_elements(first_pix) eq 0 then first_pix = 0L
if n_elements(npix_select) eq 0 then npix_select=128L < (npix-first_pix)
if n_elements(visualize) eq 0 then visualize=1
if n_elements(smooth_first) eq 0 then smooth_first = 0


n_bases_lsf1=1L
  first_pix=310L
  npix_select=89L
  n_other1=2L ;1 norm parameter
    mode_1='ga'
    run1=1L



current_tag1='Jun09_gatest1'
;First test with genetic algorithm iterative optimization (mpfit w/
;gaussian, then solber genetic algorithm)


ircsrv_ga, first_pix=first_pix, npix_select=npix_select, visualize=0, run=0, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1, n_other=n_other1


end
