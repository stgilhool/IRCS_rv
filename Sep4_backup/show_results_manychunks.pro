pro show_results_manychunks

first_first_pix=10L
output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun15_ac_manychunk_temp.fits'
npix_select=89L

;for i=0,910L do begin
;    first_pix=first_first_pix + i

;    ircsrv_allchunks, first_pix=first_pix, npix_select=npix_select, visualize=1, run=run1, current_tag=current_tag1, n_bases_lsf=n_bases_lsf1, mode=mode_1, n_other=n_other1, initial_file=initial_file, output_file=output_file, s_template_num=1

;endfor

new_output_file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/Jun15_ac_manychunk_partial.fits'

fits_info, output_file, n_ext=n_ext, /silent
fiducial_str=mrdfits(output_file,1)

r=replicate(fiducial_str, n_ext)

for i = 0, n_ext-1 do begin
    tmp=mrdfits(output_file, i+1)
    r[i]=tmp
endfor
window, 1
plot, r.delta_rv, ps=6
stop

mwrfits, r, new_output_file



stop

end
