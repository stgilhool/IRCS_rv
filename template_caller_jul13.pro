pro template_caller_jul13, _EXTRA=ex

;THIS ONE SET UP TO JUST MAKE FIRST 512 PIXELS AND THEN FILL IN 0's
;for the test

model_tag=ex.model_tag




ircsrv_maketemplate_ksjul13, first_pix=0L, npix_select=512L, _EXTRA=ex



;ircsrv_maketemplate_ksjul13, first_pix=512L, npix_select=512L, _EXTRA=ex

npix=1024L
oversamp=7L
nexp=14L
path='../data/epoch/18Jan2011/temp_results/'
fbase='GJ273_18Jan2011_AB'
alltemplates=dblarr(npix*oversamp, nexp)

for exp=0,nexp-1 do begin
    fhf=path+fbase+strtrim(exp,2)+'_0_511_'+model_tag+'.fits'
 ;   shf=path+fbase+strtrim(exp,2)+'_512_1023_'+model_tag+'.fits'
    fh=mrdfits(fhf, 1)
 ;   sh=mrdfits(shf, 1)
    
    if exp eq 0 then begin
        fhs=replicate(fh, nexp)
;        shs=replicate(sh, nexp)
        plot, fh.template[0:512*oversamp-1], /xs
    endif
    
    fhs[exp]=fh
 ;   shs[exp]=sh
    oplot, fh.template, color=200, ps=3
    ;;;STITCH TEMPLATE BACK TOGETHER
    alltemplates[0:512*oversamp-1,exp]=fh.template[0:512*oversamp-1]
    alltemplates[512*oversamp:npix*oversamp-1,exp]=1d0 ;sh.template[512*oversamp:npix*oversamp-1]
;    stop
endfor

low_reject=0.53
high_reject=3.0

;;;MAKE WL SOLN
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
;;;END MAKE WL SOLN


;;;SHIFT TEMPLATES ONTO GRID assuming exp 0 is grid
;Readin bcvs
readcol, 'GJ273bcv.dat', bcv, format='D'
delta_bcv = bcv-bcv[0]
bcv_rv=-1d0*delta_bcv
bigc=299792458d0 ;m/s 

wl_soln_arr=rebin(wl_soln_over, n_elements(wl_soln_over), nexp)
bcv_rv_arr=rebin(reform((bcv_rv*1000d0),1,nexp),n_elements(wl_soln_over),nexp)
wl_shifted=wl_soln_arr*(1d0+(bcv_rv_arr/bigc))
wl_shifted2=wl_soln_arr*(1d0-(bcv_rv_arr/bigc))
;Interpolate
help, alltemplates
help, wl_shifted
help, wl_soln_arr
;stop
template_shifted=dblarr(n_elements(wl_soln_over), nexp)
template_shifted2=dblarr(n_elements(wl_soln_over), nexp)
for exp=0,nexp-1 do begin
    
    template_shifted[*,exp]=interpol(alltemplates[*,exp], wl_shifted[*,exp], wl_soln_over)
    template_shifted2[*,exp]=interpol(alltemplates[*,exp], wl_shifted2[*,exp], wl_soln_over)
endfor

final_template_mean=mean(template_shifted, dimension=2)
final_template_med=median(template_shifted, dimension=2)

;CHECK SIGN



;Interpolate


final_template_mean2=mean(template_shifted2, dimension=2)
final_template_med2=median(template_shifted2, dimension=2)

plot, wl_soln_over, final_template_mean, /xs, title='Mean template'
oplot, wl_soln_over, final_template_mean2, color=200

window, 1
plot, wl_soln_over, final_template_med, /xs, title='Median template'
oplot, wl_soln_over, final_template_med2, color=200

;CHECK SIGN AGAIN
!p.multi=[0,1,2]
pix1=200L
pix2=400L

explist=[0,13]
;for exp=0,nexp-1 do begin
foreach exp, explist do begin
    
;    plot, wl_soln_over[pix1*oversamp:pix2*oversamp-1], template_shifted[pix1*oversamp:pix2*oversamp-1,0], /xs
plot, wl_soln_over[pix1*oversamp:pix2*oversamp-1], alltemplates[pix1*oversamp:pix2*oversamp-1,0], /xs 
   if exp gt 0 then oplot, wl_soln_over[pix1*oversamp:pix2*oversamp-1], template_shifted[pix1*oversamp:pix2*oversamp-1,exp], color=200
;    plot, wl_soln_over[pix1*oversamp:pix2*oversamp-1],
;    template_shifted2[pix1*oversamp:pix2*oversamp-1,0],/xs
    plot, wl_soln_over[pix1*oversamp:pix2*oversamp-1], alltemplates[pix1*oversamp:pix2*oversamp-1,0],/xs
    if exp gt 0 then oplot, wl_soln_over[pix1*oversamp:pix2*oversamp-1], template_shifted2[pix1*oversamp:pix2*oversamp-1,exp], color=200

    wait,1
;endfor
endforeach
!p.multi=0

;;;LOOKS LIKE SIGN 2 is better

outstr={wl_soln:wl_soln_over, $
        temp_mean:final_template_mean2, $
        temp_med:final_template_med2 $
        }
;stop, 'Write the data to disk?'
mwrfits, outstr, '../data/epoch/18Jan2011/temp_results/template_'+model_tag+'.fits', /create

;stop

end
