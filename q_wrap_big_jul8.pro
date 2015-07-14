pro q_wrap_big_jul8, make_template=make_template

if n_elements(make_template) eq 0 then make_template=0

;if n_elements(model_tag) eq 0 then model_tag=strjoin((strsplit(systime(),' ',/extract))[1:2])+'temporary'

description=''

;;;EXTRA FOR BOTH
model_tag='Jul10ksquicktestsmooth'
model_out_tag='Jul10ksfit2' ;FOR FITTING ONLY, ACTUALLY
telluric_option=3 ;0-original telluric 3-TAPAS

;;;EXTRA for template caller
min_type='amoeba'
lab_depth=0
norm=0
temp_file='../data/epoch/18Jan2011/temp_results/phoenix_template.fits'
chi2_tol=0.03
smooth_opt=1
visualize=0

;;;EXTRA for q_wrap
run=1
npix_select=128L
n_bases_lsf=1L
n_other=3L
s_template_num=9 ;template corresponding to tag
fmode='amoeba'
visualize=1

;;;KEYWORDS for q_wrap
first_first_pix=0L
pix_step=128L ;NO OVERLAP

;;;MAKE EXTRA STRS
ex={min_type:min_type, $
    lab_depth:lab_depth, $
    norm:norm, $
    temp_file:temp_file, $
    chi2_tol:chi2_tol, $
    smooth_opt:smooth_opt, $
    telluric_option:telluric_option, $
    visualize:visualize, $
    model_tag:model_tag $
    }



ex2={run:run, $
     npix_select:npix_select, $
     n_bases_lsf:n_bases_lsf, $
     n_other:n_other, $
     s_template_num:s_template_num, $
     telluric_option:telluric_option, $
     fmode:fmode, $
     visualize:visualize, $
     model_tag:model_tag $
     }



extn=tag_names(ex)
ex2tn=tag_names(ex2)

openw, lun, '/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/'+model_out_tag+'_readme.txt', /get_lun
printf, lun, model_tag + ' Description:'
printf, lun, description
printf, lun, ''
for i=0,n_elements(extn)-1 do begin
    printf, lun, extn[i],',', ex.(i)
endfor
printf, lun, ''
for i=0,n_elements(ex2tn)-1 do begin
    printf, lun, ex2tn[i], ',',ex2.(i)
endfor
printf, lun, ''
printf, lun, 'first_first_pix',',', first_first_pix
printf, lun, 'pix_step',',', pix_step

close, lun
free_lun, lun


if make_template eq 0 and file_test( '../data/epoch/18Jan2011/temp_results/template_'+model_tag+'.fits') eq 0 then message, "Template corresponding to model_tag does not exist.  Check model_tag, or turn on /make_template option"

if make_template then template_caller_jul8, _EXTRA=ex


;q_wrap_jul1, first_first_pix=first_first_pix, pix_step=pix_step, model_tag_out=model_out_tag, _EXTRA=ex2



q_wrap_jul10, first_first_pix=first_first_pix, pix_step=pix_step, model_out_tag=model_out_tag, _EXTRA=ex2

end
