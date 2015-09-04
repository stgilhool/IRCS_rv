pro sky_comp

rootpath='/home/stgilhool/RV_projects/IRCS_rv/data/'
;epochpath=rootpath+'epoch/'+epoch+'/'
;objectpath=epochpath+'final_spectra/'
;flatpath=epochpath+'final_spectra/'
;calibpath=epochpath+'calib_results/'
;outpath=epochpath+'temp_results/'
modelpath=rootpath+'supplemental/'

;Read in old telluric
h2ostr=mrdfits(modelpath+'sky_h2o.fits',1)
co2ch4str=mrdfits(modelpath+'sky_co2andch4.fits', 1)

wl_telluric_long=h2ostr.wave
h2o_long=h2ostr.trans
co2ch4_long=co2ch4str.trans

m24index=where(wl_telluric_long gt 2.275 and wl_telluric_long lt 2.365)
wl_telluric=wl_telluric_long[m24index]
h2o=h2o_long[m24index]
co2ch4=co2ch4_long[m24index]

;Read in new telluric
h2otfile=modelpath+'tapas_h2o.dat'
co2tfile=modelpath+'tapas_co2.dat'
ch4tfile=modelpath+'tapas_ch4.dat'

h2otresfile=modelpath+'tapas_h2orestest.dat'
h2otvacfile=modelpath+'tapas_h2ovac.dat'

;readcol, h2otfile, h2otwl, h2o_tapas, format='D,D'
;readcol, co2tfile, co2twl, co2_tapas, format='D,D'
;readcol, ch4tfile, ch4twl, ch4_tapas, format='D,D'

;readcol, h2otresfile, h2otreswl, h2ores_tapas, format='D,D'
;readcol, h2otvacfile, h2otvacwl, h2ovac_tapas, format='D,D'


;window, 0
;plot, wl_telluric, h2o^0.3, /xs, title='H2O with Tapas oplotted'
;oplot, h2otwl/1000d0, h2o_tapas, color=200


;ch4_tapas=interpol(ch4_tapas, ch4twl, co2twl)
;window, 1
;plot, wl_telluric, co2ch4, /xs, title='CO2CH4 with Tapas oplotted'
;oplot, co2twl/1000d0, co2_tapas*ch4_tapas, color=200

;window, 2
;h2ot=interpol(h2o_tapas, h2otwl/1000d0, wl_telluric)
;plot, wl_telluric, h2o^0.3, /xs, title='Tapas H2O interpolated onto old H2O wls'
;oplot, wl_telluric, h2ot, color=200

;window, 3
;h2otvac=interpol(h2ovac_tapas, h2otvacwl/1000d0, wl_telluric)
;plot, wl_telluric, h2o^0.3, /xs, title='Tapas H2O (vac) interpolated onto old H2O wls'
;oplot, wl_telluric, h2otvac, color=200

c2=mrdfits(modelpath+'tapas_co2vac.fits', 1)
c4=mrdfits(modelpath+'tapas_ch4vac.fits', 1)
cc=mrdfits(modelpath+'tapas_co2ch4vac.fits', 1)

;plot, cc.wavelength, cc.transmittance, /xs, title='co2ch4'
;wait,1
;oplot, c2.wavelength, c2.transmittance*c4.transmittance, color=200

print, total(c2.wavelength-c4.wavelength)
print, total(cc.wavelength-c2.wavelength)

window,1
cnew=interpol(cc.transmittance, cc.wavelength/1000d0, wl_telluric)

plot, wl_telluric, co2ch4, /xs, xr=[2.3,2.31]
;wait, 1
oplot, wl_telluric, cnew^0.5, color=200


stop
end

