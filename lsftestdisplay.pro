pro lsftestdisplay, nbases, fixed_sigma=fixed_sigma

if (n_elements(fixed_sigma) eq 0) then fixed_sigma=0
if (fixed_sigma eq 0) then begin
	if nbases eq 4 then $
		filename='./cal_results/lsf_test/model_5_4_4_13.fits' else $
	if nbases eq 3 then $
		filename='./cal_results/lsf_test/model_5_3_3_13.fits'
endif else if fixed_sigma ne 0 then begin
	if nbases eq 4 then $
        	filename='./cal_results/lsf_test_nosigchange/model_5_4_4_13.fits' else $
	if nbases eq 3 then $
        	filename='./cal_results/lsf_test_nosigchange/model_5_3_3_13.fits'
endif

fits_info, filename, n_ext=n_ext, /silent

sigma=dblarr(n_ext)
chi2perdof=dblarr(n_ext)
delta_sig=dblarr(n_ext)

!p.multi=0
npix=1024L*7L
x=dindgen(npix)
window, 0
for ext=1, n_ext do begin

	f=mrdfits(filename, ext)
	sigma[ext-1]=f.gh0_guess[0]
	chi2perdof[ext-1]=f.chi2_per_dof
	delta_sig[ext-1]=f.gh1_result[0]
	sig_line=f.gh0_guess[0]+(x*f.gh1_result[0])
	
	if ext eq 1 then plot, x, sig_line, /xs, yr=[3.5,6.5], /ys, xtitle='oversampled pixel', ytitle='sigma (x)' else oplot, x, sig_line
endfor

window, 1
plot, sigma, chi2perdof, ps=6, xr=[3.5, 6.5], yr=[0.9*min(chi2perdof),1.1*max(chi2perdof)], xtitle="sigma", ytitle="chi2/DOF"

stop


end
