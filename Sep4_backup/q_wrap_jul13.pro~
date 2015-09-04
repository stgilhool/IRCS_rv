pro q_wrap_jul10, first_first_pix=first_first_pix, pix_step=pix_step, model_out_tag=model_out_tag, _EXTRA=ex2


npix=1024L

model_tag=ex2.model_tag
;if n_elements(model_out_tag) eq 0 then model_out_tag=model_tag

initial_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_1023_Jun19_tapaspar_mpfit.fits'
new_output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/rvfit_'+model_out_tag+'_allexp.fits'
output_base='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/rvfit_'+model_out_tag


for visit=0,13 do begin


    output_file=output_base+'_' +strtrim(visit,2)+'.fits'
    for i=0,4 do begin
    
        first_pix=first_first_pix + i*pix_step
        
        ircsrv_rvjul10, first_pix=first_pix, visit=visit, output_file=output_file, initial_file=initial_file, _EXTRA=ex2

    endfor
endfor



for visit=0,13 do begin
    output_file=output_base+'_' +strtrim(visit,2)+'.fits'


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
print, new_output_file + " written to disk"

end
