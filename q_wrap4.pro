pro q_wrap4, first_pix=first_pix, npix_select=npix_select, visualize=visualize, smooth_first=smooth_first

npix=1024L

if n_elements(first_pix) eq 0 then first_pix = 0L
if n_elements(npix_select) eq 0 then npix_select=128L < (npix-first_pix)
if n_elements(visualize) eq 0 then visualize=1
if n_elements(smooth_first) eq 0 then smooth_first = 0


n_bases_lsf1=2L
run1=1L
mode_1='mpfit'
n_other1=3L
npix_select=900L
current_tag1='Jun15allchunks'
initial_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_310_398_Jun15_iter_500_150_mpfit.fits'
first_first_pix=43L
output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun15_allchunks_temp.fits'
first_first_pix=0L
output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun15_allchunks2_temp.fits'
first_first_pix=43L
output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun15_ac_bigstep_temp.fits'


;for i=0,10 do begin
;    first_pix=first_first_pix + (i*89L)
;
;    ircsrv_allchunks, first_pix=first_pix, npix_select=npix_select, visualize=1, run=run1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1, n_other=n_other1, initial_file=initial_file, output_file=output_file
;
;endfor
;
;new_output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun15_allchunks_temp.fits'
;new_output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun15_allchunks2.fits'
;new_output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun15_ac_bigstep.fits'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;FITTING FOR 89 pixel chunks all across spectrum

; first_first_pix=10L
; output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun15_ac_manychunk_temp.fits'


; for i=0,910L do begin
;     first_pix=first_first_pix + i

;     ircsrv_allchunks, first_pix=first_pix, npix_select=npix_select, visualize=1, run=run1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1, n_other=n_other1, initial_file=initial_file, output_file=output_file

; endfor

; new_output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun15_ac_manychunk.fits'

; fiducial_str=mrdfits(output_file,1)
; r=replicate(fiducial_str, 911L)

; for i = 0, 910L do begin
;     tmp=mrdfits(output_file, i+1)
;     r[i]=tmp
; endfor

; mwrfits, r, new_output_file

; plot, r.delta_rv, ps=6


; endfor

; mwrfits, r, new_output_file

; plot, r.delta_rv, ps=6


; endfor

; mwrfits, r, new_output_file

; plot, r.delta_rv, ps=6


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;REPEATING THE ABOVE BUT WITH THE NEW TEMPLATE/TELLURIC AND ALL EXPOSURES


 first_first_pix=10L
initial_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_1023_Jun19_tapaspar_mpfit.fits'




npix_select=89L
n_bases_lsf1=1L
pix_step=89L ;NO OVERLAP
n_other1=2L
s_template_num=8 ;Phoenix TEMPLATE no smooth
mode_1='amoeba'
visualize=0
new_output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun30_ac_manychunk_phoenixnosmooth.fits'


for visit=0,13 do begin

    ex={visit:visit}
    output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun30_ac_manychunk_phoenixnosmooth_temp_' +strtrim(visit,2)+'.fits'
    for i=0,900L/pix_step do begin
    
        first_pix=first_first_pix + i*pix_step
        
        ircsrv_allchunks, first_pix=first_pix, npix_select=npix_select, visualize=visualize, run=run1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1, n_other=n_other1, initial_file=initial_file, output_file=output_file, s_template_num=s_template_num, _extra=ex

    endfor
endfor



for visit=0,13 do begin

    output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun30_ac_manychunk_phoenixnosmooth_temp_' +strtrim(visit,2)+'.fits'

    if visit eq 0 then begin
        fits_info, output_file, /silent, n_ext=n_ext
        fiducial_str=mrdfits(output_file,1)

        r=replicate(fiducial_str, n_ext, 14)

    endif

    for i = 0, n_ext-1 do begin
        tmp=mrdfits(output_file, i+1)
        r[i,visit]=tmp
    endfor
    
    if visit eq 0 then plot, r[*,0].delta_rv, ps=6 else oplot, r[*,visit].delta_rv, ps=6, color=200*visit
endfor


 
 mwrfits, r, new_output_file





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;REPEATING THE ABOVE BUT WITH THE NEW TEMPLATE/TELLURIC AND ALL EXPOSURES


 first_first_pix=10L
initial_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_1023_Jun19_tapaspar_mpfit.fits'




npix_select=89L
n_bases_lsf1=1L
pix_step=89L ;NO OVERLAP
n_other1=2L
s_template_num=7 ;Phoenix TEMPLATE smoothed
mode_1='amoeba'
visualize=0
new_output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun30_ac_manychunk_phoenixsmooth.fits'


for visit=0,13 do begin

    ex={visit:visit}
    output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun30_ac_manychunk_phoenixsmooth_temp_' +strtrim(visit,2)+'.fits'
    for i=0,900L/pix_step do begin
    
        first_pix=first_first_pix + i*pix_step
        
        ircsrv_allchunks, first_pix=first_pix, npix_select=npix_select, visualize=visualize, run=run1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1, n_other=n_other1, initial_file=initial_file, output_file=output_file, s_template_num=s_template_num, _extra=ex

    endfor
endfor



for visit=0,13 do begin

    output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun30_ac_manychunk_phoenixsmooth_temp_' +strtrim(visit,2)+'.fits'

    if visit eq 0 then begin
        fits_info, output_file, /silent, n_ext=n_ext
        fiducial_str=mrdfits(output_file,1)

        r=replicate(fiducial_str, n_ext, 14)

    endif

    for i = 0, n_ext-1 do begin
        tmp=mrdfits(output_file, i+1)
        r[i,visit]=tmp
    endfor
    
    if visit eq 0 then plot, r[*,0].delta_rv, ps=6 else oplot, r[*,visit].delta_rv, ps=6, color=200*visit
endfor


 
 mwrfits, r, new_output_file

stop





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;REPEATING THE ABOVE BUT WITH THE NEW TEMPLATE/TELLURIC AND ALL EXPOSURES


 first_first_pix=10L
initial_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_1023_Jun19_tapaspar_mpfit.fits'




npix_select=89L
n_bases_lsf1=1L
pix_step=89L ;NO OVERLAP
n_other1=2L
s_template_num=6 ;Phoenix TEMPLATE
mode_1='amoeba'
visualize=0
new_output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun30_ac_manychunk_phoenixold.fits'


for visit=0,13 do begin

    ex={visit:visit}
    output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun30_ac_manychunk_phoenixold_temp_' +strtrim(visit,2)+'.fits'
    for i=0,900L/pix_step do begin
    
        first_pix=first_first_pix + i*pix_step
        
        ircsrv_allchunks, first_pix=first_pix, npix_select=npix_select, visualize=visualize, run=run1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1, n_other=n_other1, initial_file=initial_file, output_file=output_file, s_template_num=s_template_num, _extra=ex

    endfor
endfor



for visit=0,13 do begin

    output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun30_ac_manychunk_phoenixold_temp_' +strtrim(visit,2)+'.fits'

    if visit eq 0 then begin
        fits_info, output_file, /silent, n_ext=n_ext
        fiducial_str=mrdfits(output_file,1)

        r=replicate(fiducial_str, n_ext, 14)

    endif

    for i = 0, n_ext-1 do begin
        tmp=mrdfits(output_file, i+1)
        r[i,visit]=tmp
    endfor
    
    if visit eq 0 then plot, r[*,0].delta_rv, ps=6 else oplot, r[*,visit].delta_rv, ps=6, color=200*visit
endfor


 
 mwrfits, r, new_output_file

stop



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;REPEATING THE ABOVE BUT WITH THE NEW TEMPLATE/TELLURIC AND ALL EXPOSURES


 first_first_pix=10L
initial_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_1023_Jun19_tapaspar_mpfit.fits'




npix_select=89L
n_bases_lsf1=1L
pix_step=89L ;NO OVERLAP
n_other1=2L
s_template_num=6 ;Phoenix TEMPLATE
mode_1='amoeba'
visualize=1
new_output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun30_ac_manychunk_smooth.fits'


for visit=0,13 do begin

    ex={visit:visit}
    output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun30_ac_manychunk_smooth_temp_' +strtrim(visit,2)+'.fits'
    for i=0,900L/pix_step do begin
    
        first_pix=first_first_pix + i*pix_step
        
        ircsrv_allchunks, first_pix=first_pix, npix_select=npix_select, visualize=visualize, run=run1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1, n_other=n_other1, initial_file=initial_file, output_file=output_file, s_template_num=s_template_num, _extra=ex

    endfor
endfor



for visit=0,13 do begin

    output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun30_ac_manychunk_smooth_temp_' +strtrim(visit,2)+'.fits'

    if visit eq 0 then begin
        fits_info, output_file, /silent, n_ext=n_ext
        fiducial_str=mrdfits(output_file,1)

        r=replicate(fiducial_str, n_ext, 14)

    endif

    for i = 0, n_ext-1 do begin
        tmp=mrdfits(output_file, i+1)
        r[i,visit]=tmp
    endfor
    
    if visit eq 0 then plot, r[*,0].delta_rv, ps=6 else oplot, r[*,visit].delta_rv, ps=6, color=200*visit
endfor


 
 mwrfits, r, new_output_file

stop





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;REPEATING THE ABOVE BUT WITH THE NEW TEMPLATE/TELLURIC AND ALL EXPOSURES


 first_first_pix=10L
initial_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_1023_Jun19_tapaspar_mpfit.fits'




npix_select=89L
n_bases_lsf1=1L
pix_step=89L ;NO OVERLAP
n_other1=2L
s_template_num=5 ;AVG TEMPLATE
mode_1='amoeba'
visualize=0
new_output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun29_ac_manychunk_smooth.fits'


for visit=0,13 do begin

    ex={visit:visit}
    output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun29_ac_manychunk_smooth_temp_' +strtrim(visit,2)+'.fits'
    for i=0,400L/pix_step do begin
    
        first_pix=first_first_pix + i*pix_step
        
        ircsrv_allchunks, first_pix=first_pix, npix_select=npix_select, visualize=visualize, run=run1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1, n_other=n_other1, initial_file=initial_file, output_file=output_file, s_template_num=s_template_num, _extra=ex

    endfor
endfor



for visit=0,13 do begin

    output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun29_ac_manychunk_smooth_temp_' +strtrim(visit,2)+'.fits'

    if visit eq 0 then begin
        fits_info, output_file, /silent, n_ext=n_ext
        fiducial_str=mrdfits(output_file,1)

        r=replicate(fiducial_str, n_ext, 14)

    endif

    for i = 0, n_ext-1 do begin
        tmp=mrdfits(output_file, i+1)
        r[i,visit]=tmp
    endfor
    
    if visit eq 0 then plot, r[*,0].delta_rv, ps=6 else oplot, r[*,visit].delta_rv, ps=6, color=200*visit
endfor


 
 mwrfits, r, new_output_file

stop






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;REPEATING THE ABOVE BUT WITH THE NEW TEMPLATE/TELLURIC AND ALL EXPOSURES


 first_first_pix=10L
initial_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_1023_Jun19_tapaspar_mpfit.fits'




npix_select=89L
n_bases_lsf1=1L
pix_step=89L ;NO OVERLAP
n_other1=2L
s_template_num=4 ;AVG TEMPLATE
mode_1='amoeba'
visualize=0
new_output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun29_ac_manychunk_allexp.fits'


;for visit=1,13 do begin
for visit=0,0 do begin
    ex={visit:visit}
    output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun29_ac_manychunk_allexp_temp_' +strtrim(visit,2)+'.fits'
    for i=0,400L/pix_step do begin
    
        first_pix=first_first_pix + i*pix_step
        
        ircsrv_allchunks, first_pix=first_pix, npix_select=npix_select, visualize=visualize, run=run1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1, n_other=n_other1, initial_file=initial_file, output_file=output_file, s_template_num=s_template_num, _extra=ex

    endfor
endfor



for visit=0,13 do begin

    output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun29_ac_manychunk_allexp_temp_' +strtrim(visit,2)+'.fits'

    if visit eq 0 then begin
        fits_info, output_file, /silent, n_ext=n_ext
        fiducial_str=mrdfits(output_file,1)

        r=replicate(fiducial_str, n_ext, 14)

    endif

    for i = 0, n_ext-1 do begin
        tmp=mrdfits(output_file, i+1)
        r[i,visit]=tmp
    endfor
    
    if visit eq 0 then plot, r[*,0].delta_rv, ps=6 else oplot, r[*,visit].delta_rv, ps=6, color=200*visit
endfor


 
 mwrfits, r, new_output_file

stop

;;;;
;;;;
;;;;


 first_first_pix=10L
initial_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_1023_Jun19_tapaspar_mpfit.fits'




npix_select=89L
n_bases_lsf1=1L
pix_step=89L/2
n_other1=2L
s_template_num=3 ;NEW TEMPLATE
mode_1='amoeba'
visualize=0
new_output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun23_ac_manychunk_allexp.fits'


for visit=1,13 do begin
    ex={visit:visit}
    output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun23_ac_manychunk_allexp_temp_' +strtrim(visit,2)+'.fits'
    for i=0,910L/pix_step do begin
    
        first_pix=first_first_pix + i*pix_step
        
        ircsrv_allchunks, first_pix=first_pix, npix_select=npix_select, visualize=visualize, run=run1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1, n_other=n_other1, initial_file=initial_file, output_file=output_file, s_template_num=s_template_num, _extra=ex

    endfor
endfor


 
for visit=0,13 do begin

    output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun23_ac_manychunk_allexp_temp_' +strtrim(visit,2)+'.fits'

    if visit eq 0 then begin
        fits_info, output_file, /silent, n_ext=n_ext
        fiducial_str=mrdfits(output_file,1)

        r=replicate(fiducial_str, n_ext, 14)

    endif

    for i = 0, n_ext-1 do begin
        tmp=mrdfits(output_file, i+1)
        r[i,visit]=tmp
    endfor
    
    if visit eq 0 then plot, r[*,0].delta_rv, ps=6 else oplot, r[*,visit].delta_rv, ps=6, color=200*visit
endfor


 
 mwrfits, r, new_output_file

stop



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;REPEATING THE ABOVE BUT WITH THE NEW TEMPLATE/TELLURIC 


 first_first_pix=10L
initial_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_1023_Jun19_tapaspar_mpfit.fits'

  output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun22_ac_manychunk_gauss_temp.fits'
npix_select=89L
n_bases_lsf1=1L
pix_step=3L
n_other1=2L
s_template_num=3 ;NEW TEMPLATE
mode_1='amoeba'
visualize=0

 for i=0,910L/pix_step do begin
     first_pix=first_first_pix + i*pix_step

     ircsrv_allchunks, first_pix=first_pix, npix_select=npix_select, visualize=visualize, run=run1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1, n_other=n_other1, initial_file=initial_file, output_file=output_file, s_template_num=s_template_num

 endfor

 new_output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun22_ac_manychunk_gaussall.fits'


fits_info, output_file, /silent, n_ext=n_ext
 fiducial_str=mrdfits(output_file,1)

 r=replicate(fiducial_str, n_ext)

 for i = 0, n_ext-1 do begin
     tmp=mrdfits(output_file, i+1)
     r[i]=tmp
 endfor

 plot, r.delta_rv, ps=6
 
 mwrfits, r, new_output_file

stop






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;REPEATING THE ABOVE BUT WITH THE NEW TEMPLATE/TELLURIC 


 first_first_pix=10L
initial_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_1023_Jun19_tapaspar_mpfit.fits'

; output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun17_ac_manychunk_temp.fits'
  output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun19_ac_manychunk_gauss_temp.fits'
npix_select=89L
n_bases_lsf1=1L
pix_step=10L
n_other1=2L

 for i=910L-92L,910L do begin
     first_pix=first_first_pix + i+pix_step
mode_1='amoeba'
     ircsrv_allchunks, first_pix=first_pix, npix_select=npix_select, visualize=0, run=run1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1, n_other=n_other1, initial_file=initial_file, output_file=output_file, s_template_num=2

 endfor

 new_output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun19_ac_manychunk_gaussall.fits'


fits_info, output_file, /silent, n_ext=n_ext
 fiducial_str=mrdfits(output_file,1)

 r=replicate(fiducial_str, n_ext)

 for i = 0, n_ext-1 do begin
     tmp=mrdfits(output_file, i+1)
     r[i]=tmp
 endfor

 plot, r.delta_rv, ps=6
 
 mwrfits, r, new_output_file

stop

 first_first_pix=10L
initial_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_1023_Jun19_tapaspar_mpfit.fits'

; output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun17_ac_manychunk_temp.fits'
  output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun19_ac_manychunk_shorttest.fits'
npix_select=89L

pix_step=10L

 for i=11,910L/pix_step do begin
     first_pix=first_first_pix + i*pix_step
mode_1='amoeba'
     ircsrv_allchunks, first_pix=first_pix, npix_select=npix_select, visualize=1, run=run1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1, n_other=n_other1, initial_file=initial_file, output_file=output_file, s_template_num=2

 endfor

 new_output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun19_ac_manychunk.fits'

stop

 fiducial_str=mrdfits(output_file,1)
 r=replicate(fiducial_str, 911L)

 for i = 0, 910L do begin
     tmp=mrdfits(output_file, i+1)
     r[i]=tmp
 endfor

 plot, r.delta_rv, ps=6
 stop
 mwrfits, r, new_output_file


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;GETTING STARTING PARAMS FOR TAPAS
; n_other1=13L
; first_pix=0
; npix_select=1024L
; current_tag1='Jun19_tapas'

; output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_1023_Jun19_tapaspar_mpfit.fits'
;run1=0

; ircsrv_allchunks, first_pix=first_pix, npix_select=npix_select, visualize=1, run=run1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1, n_other=n_other1, initial_file=initial_file, output_file=output_file, s_template_num=2

;stop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;GETTING STARTING PARAMS FOR TAPAS
; n_other1=13L
; first_pix=0
; npix_select=1024L
; current_tag1='Jun18_tapas'

; output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_310_398_Jun18_tapaspar_mpfit.fits'

; ircsrv_allchunks, first_pix=first_pix, npix_select=npix_select, visualize=1, run=run1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1, n_other=n_other1, initial_file=initial_file, output_file=output_file, s_template_num=1



stop
end
