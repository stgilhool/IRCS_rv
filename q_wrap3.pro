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
    mode_1='mpfit'
    run1=1L
n_lsf_ga=2L

npop=10L
ngen_max=10L
ex={npop:npop, $
    ngen_max:ngen_max $
   }

current_tag1='Jun12_ga2_'+strtrim(npop,2)+'_'+strtrim(ngen_max,2)
;Testing optimal npop and ngen on exposure 0 only

;ircsrv_ga, first_pix=first_pix, npix_select=npix_select, visualize=0, run=0, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, n_lsf_ga=n_lsf_ga, mode=mode_1, n_other=n_other1, _extra=ex

n_lsf_ga=2L
npop=100L
ngen_max=10L
ex={npop:npop, $
    ngen_max:ngen_max $
   }

current_tag1='Jun12_ga2_'+strtrim(npop,2)+'_'+strtrim(ngen_max,2)
;Testing optimal npop and ngen on exposure 0 only

;ircsrv_ga, first_pix=first_pix, npix_select=npix_select, visualize=0, run=0, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, n_lsf_ga=n_lsf_ga, mode=mode_1, n_other=n_other1, _extra=ex


n_lsf_ga=2L
npop=1000L
ngen_max=10L
ex={npop:npop, $
    ngen_max:ngen_max $
   }

current_tag1='Jun12_ga2_'+strtrim(npop,2)+'_'+strtrim(ngen_max,2)
;Testing optimal npop and ngen on exposure 0 only

;ircsrv_ga, first_pix=first_pix, npix_select=npix_select, visualize=0, run=0, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, n_lsf_ga=n_lsf_ga, mode=mode_1, n_other=n_other1, _extra=ex

n_lsf_ga=2L
npop=10L
ngen_max=1000L
ex={npop:npop, $
    ngen_max:ngen_max $
   }

current_tag1='Jun12_ga2_'+strtrim(npop,2)+'_'+strtrim(ngen_max,2)
;Testing optimal npop and ngen on exposure 0 only

;ircsrv_ga, first_pix=first_pix, npix_select=npix_select, visualize=0, run=0, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, n_lsf_ga=n_lsf_ga, mode=mode_1, n_other=n_other1, _extra=ex


n_lsf_ga=2L
npop=1000L
ngen_max=1000L
ex={npop:npop, $
    ngen_max:ngen_max $
   }

current_tag1='Jun12_ga2_'+strtrim(npop,2)+'_'+strtrim(ngen_max,2)
;Testing optimal npop and ngen on exposure 0 only

;ircsrv_ga, first_pix=first_pix, npix_select=npix_select, visualize=0, run=0, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, n_lsf_ga=n_lsf_ga, mode=mode_1, n_other=n_other1, _extra=ex

n_lsf_ga=2L
npop=500L
ngen_max=150L
ex={npop:npop, $
    ngen_max:ngen_max $
   }
n_other1=3L ;THIS IS THE BIG CHANGE (adding a linear term to normalization)

current_tag1='Jun14_ga2_'+strtrim(npop,2)+'_'+strtrim(ngen_max,2)
;Testing optimal npop and ngen on exposure 0 only

;ircsrv_ga, first_pix=first_pix, npix_select=npix_select, visualize=0, run=0, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, n_lsf_ga=n_lsf_ga, mode=mode_1, n_other=n_other1, _extra=ex


n_lsf_ga=2L
npop=500L
ngen_max=150L
ex={npop:npop, $
    ngen_max:ngen_max, $
    ga_stop:1 $
   }
n_other1=3L

n_bases_lsf1=2L 

current_tag1='Jun15_iter_'+strtrim(npop,2)+'_'+strtrim(ngen_max,2)
;ITERATING ON GA SOLS

ircsrv_ga, first_pix=first_pix, npix_select=npix_select, visualize=1, run=1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, n_lsf_ga=n_lsf_ga, mode=mode_1, n_other=n_other1, initial_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/i500150.fits', _extra=ex

npop=500L
ngen_max=20L
ex={npop:npop, $
    ngen_max:ngen_max, $
    ga_stop:1 $	
   }
n_other1=3L

n_bases_lsf1=2L 

current_tag1='Jun15_iter_'+strtrim(npop,2)+'_'+strtrim(ngen_max,2)
;ITERATING ON GA SOLS

ircsrv_ga, first_pix=first_pix, npix_select=npix_select, visualize=1, run=1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, n_lsf_ga=n_lsf_ga, mode=mode_1, n_other=n_other1, initial_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/i50020.fits', _extra=ex


n_lsf_ga=4L  ;HERE's the change for this one
npop=500L
ngen_max=500L
ex={npop:npop, $
    ngen_max:ngen_max $
   }

current_tag1='Jun14_ga4_'+strtrim(npop,2)+'_'+strtrim(ngen_max,2)
;Testing optimal npop and ngen on exposure 0 only

;ircsrv_ga, first_pix=first_pix, npix_select=npix_select, visualize=0, run=0, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, n_lsf_ga=n_lsf_ga, mode=mode_1, n_other=n_other1, _extra=ex

end
