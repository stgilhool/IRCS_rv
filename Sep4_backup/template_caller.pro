pro template_caller
;Just calls the template-making function in a loop in order to fit
;each piece
npix=1024L
oversamp=7L

; for i=0,7 do begin
;     first_pix=0 > (i*128L-20L)
;     npix_select=(1024L-first_pix) < 148L
;     ircsrv_maketemplate, visualize=1, first_pix=first_pix, npix_select=npix_select
; endfor


;Stitch together the template
stellar_template_over=dblarr(1024L*7L)


for i=0,7 do begin
    first_pix_file=0 > (i*128L-20L)
    npix_select_file=(1024L-first_pix_file) < 148L   

    first_pix=i*128L
    npix_select=128L
    file='/home/stgilhool/RV_projects/IRCS_rv/data/epoch/18Jan2011/temp_results/GJ273_18Jan2011_AB0_'+strtrim(first_pix_file,2)+'_'+strtrim(first_pix_file+npix_select_file-1,2)+'_tapas2.fits'
    t=mrdfits(file, 1)
    if i eq 0 then str=replicate(t, 8)
    str[i]=t
    stellar_chunk=t.template
    stellar_template_over[first_pix*7L:(first_pix+npix_select)*7L-1]=stellar_chunk[first_pix*7L:(first_pix+npix_select)*7L-1]
endfor

mwrfits, str, 'tapas_str.fits'
stop


low_reject=0.53
high_reject=3.0





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

plot, stellar_template_over

snorm=continuum_fit(wl_soln_over, stellar_template_over, low_rej=low_reject, high_rej=high_reject)
stemp=stellar_template_over/snorm


window, 1
plot, stemp, title='normalized one'

print, 'normalize? if yes, set norm=1 and cont, else just cont'

norm=0

stop

if norm eq 1 then stellar_template_over=stemp


openw, tlun, 'template.dat', /get_lun
for line=0, n_elements(stellar_template_over)-1 do begin
    printf, tlun, wl_soln_over[line], stellar_template_over[line]
endfor
close, tlun
free_lun, tlun

stop



end
