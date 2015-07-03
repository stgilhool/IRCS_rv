pro give_to_cullen
;arranging wl_soln for cullen




npix=1024L
oversamp=7l

epoch='18Jan2011'
rootpath='/home/stgilhool/RV_projects/IRCS_rv/data/'
epochpath=rootpath+'epoch/'+epoch+'/'
calibpath=epochpath+'calib_results/'

calib_file='GJ273_18Jan2011_AB1_5_7_7_13.fits'
calib_ext=13
model_file=calibpath+calib_file
model_par=mrdfits(model_file, calib_ext)


wl_coeff=model_par.wl_result
;;;Contstruct wl array and oversampled wl array from coeffs
x=dindgen(npix)
xx=(dindgen(npix*oversamp)-(oversamp/2))/oversamp

wl_soln=poly(x, wl_coeff)
;samp_index=(lindgen(npix)*oversamp)+(oversamp/2)

wl_soln_over=poly(xx, wl_coeff)

outstr={wl_soln:wl_soln, wl_soln_over:wl_soln_over, wl_coeff:wl_coeff}

mwrfits, outstr, '../data/wl_soln.fits'

end


