pro q_wrap2, first_pix=first_pix, npix_select=npix_select, visualize=visualize, smooth_first=smooth_first

npix=1024L

if n_elements(first_pix) eq 0 then first_pix = 0L
if n_elements(npix_select) eq 0 then npix_select=128L < (npix-first_pix)
if n_elements(visualize) eq 0 then visualize=1
if n_elements(smooth_first) eq 0 then smooth_first = 0


current_tag1='Mar26_wl_test'
;This one has many wl guesses, 1 normalization parameter, and fixed RV
;for exposures 10, 11, and 12 (11 was the big outlier), delta_wl.relstep=1d-3

current_tag1='Mar26_wl_test2'
;This one has many wl guesses, 1 normalization parameter, and FREE RV
;for exposures 10, 11, and 12 (11 was the big outlier), delta_wl.relstep=1d-3

current_tag1='Mar26_wl_allexp'
;Same as test2 (free rv) but for all exposures

current_tag1='Mar26_wl_allexp2'
;Same as previous but with larger (but less fine) wl steps.

current_tag1='Mar26_wl_allexp3'
;Same but smaller wl step 1d-4, only -1d-6 and 1d-6 as guesses, and
;mpside =2 for rv and wl

current_tag1='Mar26_wl_allexp4'
;Same as 3 but with RV fixed.  I'm specifically looking at visit 7
;which is the biggest outlier, but with the lowest chi2.  The wl shift
;is ~1d-7 and all surrounding exposures are ~3d-6. Also just turned
;off limits on wl

current_tag1='Mar26_wl_allexp5'
;RV free, wl fixed at 0.  WHOOPS, not fixed.  Results of this one were
;pretty bad, but it should have been the same (or close to the same)
;as allexp3.  I guess since the starting wl guess was 0, that screwed
;things up

current_tag1='Mar26_wl_allexp6'
;now with wl fixed at 0. Also turned off rv limits

current_tag1='Mar27_wlrv_allbrute'
;looping through fixed RV and wls to look at chi2 space

current_tag1='Mar31_wlrv_allbrute7'
;looping through fixed RV and wls to look at chi2 space, only on exp
;7.  limits are changed so that we can include both minima (it was cut
;off in the brute force run)
n_bases_lsf1=1L


current_tag1='Apr01_lsf_brute'
;loopling through fixed RV and wls for exp 1,2,3 and with a 2-basis
;lsf  

current_tag1='Apr03_lsf_brute'
;loopling through fixed RV and wls for all exp and with a 2-basis lsf

current_tag1='Apr03_lsf_brutemore'
;looping through fixed RV and wls for all exp and with a 2-basis lsf

current_tag1='Apr03_lsf3_brute'
;loopling through fixed RV and wls for all exp and with a 3-basis lsf

current_tag1='Apr04_lsf_sig'
;fixed RV and wls for all exp and with 3-basis lsf with better sig
;guess (0.65 instead of 0.75)

current_tag1='Apr04_lsf_sigplus'
;fixed RV and wls for all exp and with 2-basis lsf with better sig
;guess (0.65 instead of 0.75), and 1d-3 for other basis guess
 
current_tag1='Apr04_lsf_sigminus'
;fixed RV and wls for all exp and with 2-basis lsf with better sig
;guess (0.65 instead of 0.75), and -1d-3 for other basis guess

current_tag1='Apr07_lsf2_brute'
;looping through sigma and gh0 parameters fixed (NEW BRUTE FORCE)

current_tag1='Apr20_lsf2_fine'
;looping through a finer grid of sig and GH0, for just the first
;exposure
n_bases_lsf1=2L
  first_pix=310L
  npix_select=89L
  n_other1=2L ;1 norm parameter
    mode_1='mpfit'



current_tag1='Apr21_lsf3_brute'
;looping through gh0_1 and gh0_2

current_tag1='Apr22_lsf3_wide'
;looping through a wider range of gh0_1 and _2

current_tag1='Apr22_lsf3_more'
;looping through a wider range of gh0_1 and _2

current_tag1='Apr23_lsf3_many'
;looping through a wider range of gh0_1 and _2



n_bases_lsf1=3L
  first_pix=310L
  npix_select=89L
  n_other1=2L ;1 norm parameter
    mode_1='mpfit'
  


ircsrv_rvguesstest, first_pix=first_pix, npix_select=npix_select, visualize=0, run=0, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1, n_other=n_other1


end
