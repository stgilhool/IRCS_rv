pro degen_test, file1, file2

a=mrdfits(file1, 1)
b=mrdfits(file2, 1)

npix_over=n_elements(a.lsf[0,*])
npix_lsf=n_elements(a.lsf[*,0])

wl_diff=a.trial_soln_over-b.trial_soln_over

awl_der=dblarr(npix_over)
bwl_der=dblarr(npix_over)

print, "starting first for loop"

for i=0, npix_over-1 do begin
	if i eq 0 then begin
		awl_der[i]=a.trial_soln_over[i+1]-a.trial_soln_over[i]
		bwl_der[i]=b.trial_soln_over[i+1]-b.trial_soln_over[i]
	endif else if i eq npix_over-1 then begin
		awl_der[i]=a.trial_soln_over[i]-a.trial_soln_over[i-1]
		bwl_der[i]=b.trial_soln_over[i]-b.trial_soln_over[i-1]
	endif else begin
		awl_der[i]=(a.trial_soln_over[i+1]-a.trial_soln_over[i-1])/2d0
		bwl_der[i]=(b.trial_soln_over[i+1]-b.trial_soln_over[i-1])/2d0
	endelse
endfor

print, "starting second for loop"

apeak_int=dblarr(npix_over)
bpeak_int=dblarr(npix_over)

for i=0, npix_over-1 do begin
	if (i lt 40) or (i gt npix_over-40) then begin
		apeak_int[i]=0
		bpeak_int[i]=0
	endif else begin
		amem=0d0
		bmem=0d0
		astop=0
		bstop=0
		print, i
		for j=28, npix_lsf-29 do begin
			if amem eq 0 and bmem eq 0 then begin
				atot=total(a.lsf[0:j,i])	
				btot=total(b.lsf[0:j,i])
			;	print, atot, amem
			endif else begin
				atot=amem+a.lsf[j,i]
				btot=btot+b.lsf[j,i]
			;	print, atot, amem
			endelse

			if atot lt 0.5 then amem=atot $
				else if astop eq 0 then begin
					alsfdiff=atot-amem
					alsfdiffmem=0.5-amem
					alsffracdiff=alsfdiffmem/alsfdiff
					apeak_int[i]=j + alsffracdiff - 35L
					astop=1
				endif

			if btot lt 0.5 then bmem=btot $
                                else if bstop eq 0 then begin
                                        blsfdiff=btot-bmem
                                        blsfdiffmem=0.5-bmem
                                        blsffracdiff=blsfdiffmem/blsfdiff
                                        bpeak_int[i]=j + blsffracdiff - 35L
					bstop=1           
                                endif
			if astop eq 1 and bstop eq 1 then break
		endfor
		
	endelse
endfor

print, "plotting"

lsf_diff=apeak_int-bpeak_int

;wl_diff_pix=(a.trial_soln_over/awl_der)-(b.trial_soln_over/bwl_der)
wl_diff_pix=wl_diff/awl_der


plot, lsf_diff+wl_diff_pix

stop

end		
