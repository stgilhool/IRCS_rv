
pro lsfcompare

filename='~/RV_projects/IRCS_rv/cal_results/Oct31/model_4_4_4_13.fits'

fits_info, filename, /silent, N_ext=n_ext

a=ptrarr(n_ext, /allocate_heap)

print, "total num of ext : ", n_ext

for extension=0, n_ext-1 do begin
	*(a[extension])=mrdfits(filename, extension+1, /silent)
	status=tag_exist((*(a[extension])), 'min_type')
	if status eq 1 then begin
		print, "Ext: ", extension+1, " Min_type: " + (*(a[extension])).min_type, " chi2: ", (*(a[extension])).chi2
	endif else begin
		print, "Ext: ", extension+1, " Min_type: Unknown", " chi2: ", (*(a[extension])).chi2
	endelse
endfor

end
