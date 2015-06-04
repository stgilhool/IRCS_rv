pro quick_wrapper, first_pix=first_pix, npix_select=npix_select, visualize=visualize, smooth_first=smooth_first

npix=1024L

if n_elements(first_pix) eq 0 then first_pix = 0L
if n_elements(npix_select) eq 0 then npix_select=128L < (npix-first_pix)
if n_elements(visualize) eq 0 then visualize=1
if n_elements(smooth_first) eq 0 then smooth_first = 0


; current_tag1='Mar19_tnmin2'
; n_bases_lsf1=2L
; first_pix=256L
; npix_select=128L
; temp_num1=34L
; mode_1='tnmin'
; description='This is a 2 basis run on 256-383 with tnmin instead of mpfit. '

; current_tag1='Mar19_mpfit2'
; n_bases_lsf1=2L
; first_pix=256L
; npix_select=128L
; temp_num1=35L
; mode_1='mpfit'
; description='This is a 2 basis run on 256-383 with mpfit and using the relative step and scaled parameters. '

;  current_tag1='Mar19_mpfit2_5d-3'
;  n_bases_lsf1=2L
;  first_pix=256L
;  npix_select=128L
;  temp_num1=36L
;  mode_1='mpfit'
;  description='This is a 2 basis run on 256-383 with mpfit and relstep of 5d-3 on gh0. '

;   current_tag1='Mar19_mpfit2_1d-2'
;   n_bases_lsf1=2L
;   first_pix=256L
;   npix_select=128L
;   temp_num1=37L
;   mode_1='mpfit'
;   description='This is a 2 basis run on 256-383 with mpfit and relstep of 1d-2 on gh0. '

;   current_tag1='Mar19_mpfit2_5d-2'
;   n_bases_lsf1=2L
;   first_pix=256L
;   npix_select=128L
;   temp_num1=38L
;   mode_1='mpfit'
;   description='This is a 2 basis run on 256-383 with mpfit and relstep of 5d-2 on gh0. '

;   current_tag1='Mar19_mpfit2_1d-1'
;   n_bases_lsf1=2L
;   first_pix=256L
;   npix_select=128L
;   temp_num1=39L
;   mode_1='mpfit'
;   description='This is a 2 basis run on 256-383 with mpfit and relstep of 1d-1 on gh0. '

;   current_tag1='Mar19_mpfit2_5d-1'
;   n_bases_lsf1=2L
;   first_pix=256L
;   npix_select=128L
;   temp_num1=40L
;   mode_1='mpfit'
;   description='This is a 2 basis run on 256-383 with mpfit and relstep of 5d-1 on gh0. '
;  current_tag1='Mar23_mpfit2_wl_gh'
current_tag1='display'
n_bases_lsf1=2L
  first_pix=290L
  npix_select=128L
  temp_num1=44L
  mode_1='mpfit'
  description='This is a 2 basis run on 290-419 with mpfit and relstep of 1d-1 on gh0, and relstep of 1d-1 on gh1.  Changed the negativity penalty because the numbers were blowing up. And this one has the wl coeff step sizes tweaked (1d-14, 1d-18). and initial sigma is set to 0.6 with new absolute step size'


ircsrv_rvfit, first_pix=first_pix, npix_select=npix_select, visualize=1, run=0, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1
stop
ircsrv_rvfit, first_pix=first_pix, npix_select=npix_select, visualize=0, run=1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1
;stop
s_input_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_'+strtrim(first_pix,2)+'_'+strtrim(first_pix+npix_select-1, 2)+'_'+current_tag1+'_r1.fits'


ircsrv_templatefit, first_pix=first_pix, npix_select=npix_select, input_file=s_input_file, s_template_numout=temp_num1



openw, lun, '/home/stgilhool/RV_projects/IRCS_rv/data/smooth_penalty_test/test'+strtrim(temp_num1,2)+'/description.txt', /get_lun
printf, lun, description + s_input_file
close, lun
free_lun, lun

ircsrv_rvfit, first_pix=first_pix, npix_select=npix_select, visualize=0, run=0, s_iter=1, s_template_num=temp_num1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1

ircsrv_rvfit, first_pix=first_pix, npix_select=npix_select, visualize=0, run=1, s_iter=1, s_template_num=temp_num1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1



stop
;;;;;
current_tag1='Mar18_gaussian'
n_bases_lsf1=1L
temp_num1=26L

ircsrv_rvfit, first_pix=first_pix, npix_select=npix_select, visualize=1, run=0, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1
stop
ircsrv_rvfit, first_pix=first_pix, npix_select=npix_select, visualize=0, run=1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1

ircsrv_templatefit, first_pix=first_pix, npix_select=npix_select, input_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_256_383_Mar18_gaussian_r1.fits', s_template_numout=temp_num1

ircsrv_rvfit, first_pix=first_pix, npix_select=npix_select, visualize=1, run=0, s_iter=1, s_template_num=temp_num1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1

ircsrv_rvfit, first_pix=first_pix, npix_select=npix_select, visualize=0, run=1, s_iter=1, s_template_num=temp_num1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1

stop
;4bases

ircsrv_rvfit, first_pix=first_pix, npix_select=npix_select, visualize=0, run=0, current_tag='Mar18_4lsf', n_bases_lsf=4L, initial_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_256_383_Mar18_3lsf_r3.fits'

ircsrv_rvfit, first_pix=first_pix, npix_select=npix_select, visualize=0, run=1, current_tag='Mar18_4lsf', n_bases_lsf=4L

ircsrv_templatefit, first_pix=first_pix, npix_select=npix_select, input_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_256_383_Mar18_4lsf_r1.fits', s_template_numout=27L

ircsrv_rvfit, first_pix=first_pix, npix_select=npix_select, visualize=0, run=0, s_iter=1, s_template_num=27L, current_tag='Mar18_4lsf', n_bases_lsf=4L

ircsrv_rvfit, first_pix=first_pix, npix_select=npix_select, visualize=0, run=1, s_iter=1, s_template_num=27L, current_tag='Mar18_4lsf', n_bases_lsf=4L

;5bases

ircsrv_rvfit, first_pix=first_pix, npix_select=npix_select, visualize=0, run=0, current_tag='Mar18_5lsf', n_bases_lsf=5L, initial_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_256_383_Mar18_4lsf_r3.fits'

ircsrv_rvfit, first_pix=first_pix, npix_select=npix_select, visualize=0, run=1, current_tag='Mar18_5lsf', n_bases_lsf=5L

ircsrv_templatefit, first_pix=first_pix, npix_select=npix_select, input_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_256_383_Mar18_5lsf_r1.fits', s_template_numout=28L

ircsrv_rvfit, first_pix=first_pix, npix_select=npix_select, visualize=0, run=0, s_iter=1, s_template_num=28L, current_tag='Mar18_5lsf', n_bases_lsf=5L

ircsrv_rvfit, first_pix=first_pix, npix_select=npix_select, visualize=0, run=1, s_iter=1, s_template_num=28L, current_tag='Mar18_5lsf', n_bases_lsf=5L

end
