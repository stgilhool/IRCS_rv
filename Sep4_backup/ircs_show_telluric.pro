pro ircs_show_telluric

npix=1024L
oversamp=3L
npix_select=128L
first_pix=0L
npix_model=npix_select*oversamp


rootpath='/home/stgilhool/RV_projects/IRCS_rv/data/'
calib_file='GJ273_18Jan2011_AB1_5_7_7_13.fits'
calib_ext=13
;model_file=calibpath+calib_file
model_file=rootpath+'epoch/18Jan2011/calib_results/'+calib_file
model_par=mrdfits(model_file, calib_ext)



modelpath=rootpath+'supplemental/'
;READ IN TELLURIC
h2ostr=mrdfits(modelpath+'sky_h2o.fits',1)
co2ch4str=mrdfits(modelpath+'sky_co2andch4.fits', 1)

wl_telluric_long=h2ostr.wave
h2o_long=h2ostr.trans
co2ch4_long=co2ch4str.trans

m24index=where(wl_telluric_long gt 2.275 and wl_telluric_long lt 3.365)
wl_telluric=wl_telluric_long[m24index]
h2o=h2o_long[m24index]
co2ch4=co2ch4_long[m24index]


x=dindgen(npix)
xx=(dindgen(npix*oversamp)-(oversamp/2))/oversamp

wl_coeff=model_par.wl_result

wl_soln=poly(x, wl_coeff)
wl_soln_over=poly(xx, wl_coeff)

;Redefine those vectors to reflect only the first 128 pixels
x_select=dindgen(npix_select)+first_pix
xx_select=(dindgen(npix_model)-(oversamp/2))/oversamp+first_pix

wl_soln_select=poly(x_select, wl_coeff)
wl_soln_over_select=poly(xx_select, wl_coeff)


!p.multi=[0,1,2]

h2o_plot=interpol(h2o, wl_telluric, wl_soln_over_select)
co2ch4_plot=interpol(co2ch4, wl_telluric, wl_soln_over_select)

plot, h2o_plot, /xs
plot, co2ch4_plot, /xs

stop

h2o_plot=interpol(h2o, wl_telluric, wl_soln_over)
co2ch4_plot=interpol(co2ch4, wl_telluric, wl_soln_over)

plot, h2o_plot, /xs
plot, co2ch4_plot, /xs

stop





end
