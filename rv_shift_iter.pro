function rvfit, p

common rvinfo, wl_soln, wl_soln_full, samp_index, wl_soln_over_full, int_lab_over_full, int_obs, lsf, npix_lsf, npix_model, first_pix, npix_select, oversamp, norm_factor, err, visualize, npix_trim_start, npix_trim_end, min_type, nparam_other, last_guess, nlines, COwl_shift, h2o, co2ch4, tellwl, visit

last_guess=p

delta_rv=p[0]
co2ch4_depth=p[1]
h2o_depth=p[2]

stellar_depth=p[3:3+nlines-1]
stellar_width=p[(3+nlines):-1]



bigc=299792458d0 ;m/s



;;;;;;;;;;;;;;;
wl_shifted=wl_soln_over_full*(1d0+(delta_rv/bigc))


;Make stellar template on oversampled wavelength grid
line_strength=dblarr(nlines, n_elements(wl_soln_over_full))

for lin=0, nlines-1 do begin
    line_strength[lin, *]=stellar_depth[lin]*exp(-0.5d0*((wl_soln_over_full-COwl_shift[lin])/stellar_width[lin])^2)
    ;stop
endfor

stellar_trans_decomp=1d0-line_strength
stellar_spec=product(stellar_trans_decomp, 1)

;stop
star_model=interpol(stellar_spec, wl_shifted, wl_soln_over_full)

;stop

;Deal with telluric
;scale telluric depth
h2o_depth=h2o^h2o_depth
co2ch4_depth=co2ch4^co2ch4_depth
;shift lines and spline them onto soln
;tellwl_shift=tellwl*(1d0+(tell_rv/bigc))
h2o_shift=interpol(h2o_depth, tellwl, wl_soln_over_full)
co2ch4_shift=interpol(co2ch4_depth, tellwl, wl_soln_over_full)

;;;;;;;
;multiply and convolve
product_model=star_model*int_lab_over_full*h2o_shift* co2ch4_shift
obs_noconv=product_model
;convolve

;;;;
;Convolve model with lsf
;;;;

int_conv=dblarr(npix_model)
int_conv_mtx=dblarr(npix_lsf, npix_model)


for index=0,npix_model-1 do begin
    model_index=index+first_pix*oversamp
    if model_index lt (npix_lsf/2) or model_index gt (n_elements(obs_noconv)-1-(npix_lsf/2)) then $
      int_conv_mtx[*,index]=replicate(1d0, npix_lsf) $
    else int_conv_mtx[*,index]=obs_noconv[model_index-(npix_lsf/2):model_index+(npix_lsf/2)]
endfor
int_conv_temp=int_conv_mtx*lsf
int_conv=total(int_conv_temp, 1, /double)


;;;;
;Down-sample the model to IRCS resolution
;;;;

npix_tophat=oversamp
tophat=replicate(1d0/npix_tophat,npix_tophat)
int_avg=convol(int_conv, tophat)

star_model_avg=convol(star_model, tophat)

int_lab_over_avg=convol(int_lab_over_full, tophat)
samp_index_full=lindgen(n_elements(wl_soln_full))*oversamp+(oversamp/2)
star_model_down=star_model_avg[samp_index_full]
int_lab_down=int_lab_over_avg[samp_index_full]
h2o_shift_avg=convol(h2o_shift, tophat)
h2o_shift_down=h2o_shift_avg[samp_index_full]
h2o_shift_dim=h2o_shift_down[first_pix:first_pix+npix_select-1]
co2ch4_shift_avg=convol(co2ch4_shift, tophat)
co2ch4_shift_down=co2ch4_shift_avg[samp_index_full]
co2ch4_shift_dim=co2ch4_shift_down[first_pix:first_pix+npix_select-1]
star_model_dim=star_model_down[first_pix:first_pix+npix_select-1]
int_lab_dim=int_lab_down[first_pix:first_pix+npix_select-1]


int_model=int_avg[samp_index] 


;;;;
;Correct subtle variations in continuum
;;;;

;if nparam_other gt 0 then int_model=int_model*norm_factor



;;;;
;;;Get chi^2
;;;;

int_obs_dim=int_obs[first_pix:first_pix+npix_select-1]
err_dim=err[first_pix:first_pix+npix_select-1]
nan=where(finite(int_model) eq 0, nancount)
if nancount gt 0 then begin
    err_dim[nan]=1d-5
    int_model[nan]=1d0
endif

res=(int_obs_dim-int_model)
res2=res^2
dev=res/err_dim
dev2=dev^2


chi1_vec=dev[npix_trim_start:-1*(npix_trim_end+1)]
chi2_vec=dev2[npix_trim_start:-1*(npix_trim_end+1)]
chi2=total(chi2_vec, /double, /nan)



;;;;
;If visualize is set, do animation of fitting spectrum for full spectrum
;;;;

if (visualize eq 1) then begin
    
    ;title_str="RV shift: "+strtrim(delta_rv,2)+" | tell_rv: " + strtrim(tell_rv,2)+" | visit: " + strtrim(visit,2) + " | chi2: "+strtrim(chi2,2)
;    title_str="RV shift: "+strtrim(delta_rv,2)+" | visit: " + strtrim(visit,2)+" | depth: " + strtrim(temp_depth_coeff,2) + " | chi2: "+strtrim(chi2,2)

;    plot, wl_soln, int_obs_dim,yr=[0.1,1.1], /xs,
;    title=title_str,xtitle="Wavelength (microns)", charsize=1.0

;    plot, wl_soln, int_lab_dim + 2d0, yr=[0,3.5], /xs,
;    title=title_str,xtitle="Wavelength (microns)", charsize=1.2 
;window, 0, xsize=500, ysize=500
phi=findgen(32)*(!PI*2/32.)
phi = [ phi, phi(0) ]
usersym, cos(phi), sin(phi), /fill



;set_plot, 'ps'
set_plot, 'x'
;device, /encapsulated, filename = 'rv_fit.eps'
;!p.font=0
;device, /color, bits=8
;loadct, 12


    plot, wl_soln, int_lab_dim + 1.6d0, yr=[0,2.8], /ys, /xs, ytickname=replicate(' ', 5), yticks=1,yminor=1, xtitle="Wavelength (microns)", charsize=1.2 
    oplot, wl_soln, star_model_dim+1.0d0
    oplot, wl_soln, h2o_shift_dim+0.6d0
    oplot, wl_soln, co2ch4_shift_dim+0.6d0, color=200
    oplot, wl_soln, int_obs_dim
    oplot, wl_soln, int_model, color=200, ps=8, symsize=0.6 
    oplot, wl_soln, res + 0.25, ps=3
 ;device,/close   
;    plot, wl_soln, res, title="Residuals | RMS : ", xtitle="Wavelength (microns)", yr=[-0.1, 0.1], ps=3, /xs
    ;plot, wl_soln, res, title="Residuals | RMS : ", xtitle="Wavelength (microns)", ps=3, /xs

 ;  plot, wl_soln, star_model_dim, /xs

  ; plot, wl_soln, telluric_shift_dim, /xs

  ;  oplot, wl_soln*(1D0+(delta_rv/bigc)), temp_select, color=200
    
endif

;stop

;;;;
;Report results
;;;;



print, chi2
print, "----------"


if min_type eq 'amoeba' then return, chi2 $
else if min_type eq 'mpfit' then return, chi1_vec $
else message, 'incorrect min_type... cannot return a chi2 value'


end


pro rv_shift_iter, epoch, object, trace, visualize=visualize, first_pix=first_pix, npix_select=npix_select, min_type=min_type

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;MAIN BODY ;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


npix=1024L
if n_elements(visualize) eq 0 then visualize=0
;Define keywords for fitting just a range
if n_elements(first_pix) eq 0 then first_pix=0
if n_elements(npix_select) eq 0 then npix_select=npix-first_pix
if n_elements(min_type) eq 0 then min_type='mpfit'


;set file paths
rootpath='/home/stgilhool/RV_projects/IRCS_rv/data/'
epochpath=rootpath+'epoch/'+epoch+'/'
objectpath=epochpath+'final_spectra/'
flatpath=epochpath+'final_spectra/'
calibpath=epochpath+'calib_results/'
outpath=epochpath+'temp_results/'
modelpath=rootpath+'supplemental/'
outputpath=epochpath+'rv_results/'

;p=paths(epoch=epoch, object=object, trace=trace)
;rootpath=p.dataroot
;epochpath=p.epochfold
;objectpath=p.specfold
;flatpath=p.specfold
;outpath=p.tempfold
;modelpath=p.dataroot+'supplemental/'


;Define flat file names and full path file names
ABflat_filename=strjoin(['NH3', object, epoch, 'AB1'], '_')+'.fits'
BAflat_filename=strjoin(['NH3', object, epoch, 'BA1'], '_')+'.fits'
ABflat_file=flatpath+strjoin(['NH3', object, epoch, 'AB1'], '_')+'.fits'
BAflat_file=flatpath+strjoin(['NH3', object, epoch, 'BA1'], '_')+'.fits'


;Define object file names and full path file names
ABobj_listname=strjoin([object, epoch, 'AB'], '_')+'.list'
BAobj_listname=strjoin([object, epoch, 'BA'], '_')+'.list'
ABobj_list=objectpath+strjoin([object, epoch, 'AB'], '_')+'.list'
BAobj_list=objectpath+strjoin([object, epoch, 'BA'], '_')+'.list'


;Read in object filenames and then spectra
readcol, ABobj_list, ABobj_filename, format='A'
readcol, BAobj_list, BAobj_filename, format='A'

ABobj_file=objectpath+ABobj_filename
BAobj_file=objectpath+BAobj_filename

n_ABobj=n_elements(ABobj_file)
n_BAobj=n_elements(BAobj_file)


;Initialize arrays of spectra
temp_str=mrdfits(ABobj_file[0], 1)
temp_spectrum=temp_str.spectrum
npix=n_elements(temp_spectrum)
ABspec_arr=dblarr(n_ABobj, npix)
BAspec_arr=dblarr(n_BAobj, npix)
ABflat_arr=dblarr(n_ABobj, npix)
BAflat_arr=dblarr(n_BAobj, npix)


ABmjd_arr=dblarr(n_ABobj)
BAmjd_arr=dblarr(n_BAobj)

;FIX (ADD ERROR ARRAYS)



;Set continuum normalization params
low_reject=0.53
high_reject=3.0

;Loop through each object file and flip and normalize it
for f=0, n_ABobj-1 do begin
    ABobj_spec=mrdfits(ABobj_file[f], 1)
    BAobj_spec=mrdfits(BAobj_file[f], 1)

    ABspec_backwards=ABobj_spec.spectrum
    BAspec_backwards=BAobj_spec.spectrum

    ;Flip the spectra
    ABspec=reverse(ABspec_backwards)
    BAspec=reverse(BAspec_backwards)

    ;FIX remove this later but output SNR
    ABerr=reverse(ABobj_spec.sigma)
    BAerr=reverse(BAobj_spec.sigma)
    ABSNR=mean(ABspec/ABerr)
    BASNR=mean(BAspec/BAerr)
    print, '--------------------'
    print, 'AB SNR: '
    print, ABSNR
    print, '--------------------'
    print, 'BA SNR: '
    print, BASNR
    print, '--------------------'

    ;Normalize the model spectrum and the observed spectrum
    ABnorm=continuum_fit(dindgen(npix), ABspec, low_rej=low_reject, high_rej=high_reject)
    BAnorm=continuum_fit(dindgen(npix), BAspec, low_rej=low_reject, high_rej=high_reject)

    ABspec_arr[f,*]=ABspec/ABnorm
    BAspec_arr[f,*]=BAspec/BAnorm

    ;Get MJD
    ABhead=ABobj_spec.header
    BAhead=BAobj_spec.header
    
    ABmjd_arr[f]=sxpar(ABhead, 'MJD')
    BAmjd_arr[f]=sxpar(BAhead, 'MJD')

endfor

;Read flat files and flip and normalize them too    

ABflat_spec=mrdfits(ABflat_file, 1)
BAflat_spec=mrdfits(BAflat_file, 1)

ABflat_backwards=ABflat_spec.NH3_spectrum
BAflat_backwards=BAflat_spec.NH3_spectrum

;Flip
ABflat=reverse(ABflat_backwards)
BAflat=reverse(BAflat_backwards)

;Normalize
ABflatnorm=continuum_fit(dindgen(npix), ABflat, low_rej=low_reject, high_rej=high_reject)
BAflatnorm=continuum_fit(dindgen(npix), BAflat, low_rej=low_reject, high_rej=high_reject)

;Make array
ABflat_arr=rebin(reform(ABflat/ABflatnorm, 1, npix), n_ABobj, npix)
BAflat_arr=rebin(reform(BAflat/BAflatnorm, 1, npix), n_BAobj, npix)


;;;;;;;;
;Divide spec by flat
ABtemp1=ABspec_arr/ABflat_arr
BAtemp1=BAspec_arr/BAflat_arr

;Take median
ABtemp_rough=median(ABtemp1, dimension=1, /double)
BAtemp_rough=median(BAtemp1, dimension=1, /double)

;;;;;;;

;Read in calibration results
calib_file='GJ273_18Jan2011_AB1_5_7_7_13.fits'
calib_ext=13
model_file=calibpath+calib_file
model_par=mrdfits(model_file, calib_ext)


common rvinfo, wl_soln, wl_soln_full, samp_index, wl_soln_over_full, int_lab_over_full, int_obs, lsf, npix_lsf, npix_model, firstpix, npixselect, oversamp, norm_factor, err, visual, npix_trim_start, npix_trim_end, mintype, nparam_other, last_guess, nlines, COwl_shift, h2o, co2ch4, tellwl, visit

visual=visualize
npixselect=npix_select
firstpix=first_pix
mintype=min_type

;Assume lin_switch is on
lin_switch=1

;;;Parameters from amoeba
wl_coeff=model_par.wl_result

gh_coeff=model_par.gh0_result
gh_coeff_out=gh_coeff
sigma=gh_coeff[0]
gh_coeff[0]=1d0


if lin_switch eq 1 then begin
    gh_lin_coeff=model_par.gh1_result
    gh_lin_coeff_out=gh_lin_coeff
    sigma_lin=gh_lin_coeff[0]
    gh_lin_coeff[0]=0d0
    other=model_par.other_result
endif else begin
    sigma_lin=0d0
endelse

tau_scale=other[0]
if n_elements(other) gt 1 then norm_pts=other[1:n_elements(other)-1]

nparam=model_par.nparam
nparam_wl=n_elements(wl_coeff)
nparam_gh=n_elements(gh_coeff)
nparam_other=n_elements(other)
nparam_gh1=nparam-nparam_wl-nparam_gh-nparam_other


if n_elements(other) gt 1 then begin
    nparam_norm=nparam_other-1
    x_x=lindgen(npix_select)+first_pix
    x_norm=dblarr(nparam_norm)
    x_norm[0]=first_pix
    x_norm[nparam_norm-1]=first_pix+npix_select-1
    x_increment=npix_select/(nparam_norm-2)
    for i=1, nparam_norm-2 do x_norm[i]=x_norm[i-1]+x_increment

    norm_factor=interpol(norm_pts, x_norm, x_x, /spline)    
    
endif
;;;

;Other parameters and constants
oversamp=model_par.oversamp
npix_lsf=(oversamp*10L)+1L
c=299792458.D

;Grab other structure variables just in case
;npix_trim_start=model_par.npix_trim_start
;npix_trim_end=model_par.npix_trim_end
npix_trim_start=5L
npix_trim_end=5L

npix_model=npix_select*oversamp






;read in the crazy template I made
;tempstr=mrdfits(outpath+'GJ273_18Jan2011_AB2_0_1023_template.fits', 1)
;temp=tempstr.template

;temp_select=temp[first_pix:first_pix+npix_select-1]



;Read in LAB SPECTRUM
modelfile= modelpath+'NH3_model.dat'
readcol, modelfile, wl_lab, int_lab, format='D,D'

norm_lab=continuum_fit(wl_lab, int_lab, low_rej=low_reject, high_rej=high_reject)
int_lab_copy=int_lab
int_lab=int_lab_copy/norm_lab


;Adjust the optical depth of the lab spectrum (I think this is okay to
;do for the high-resolution spectrum even though the wl scale isn't linear)
int_lab_depth=int_lab^tau_scale


;read in telluric spectrum
;tell=mrdfits(modelpath+'skytable.fits', 1)
;tellwl=tell.lam
;telluric=tell.trans
h2ostr=mrdfits(modelpath+'sky_h2o.fits',1)
co2ch4str=mrdfits(modelpath+'sky_co2andch4.fits', 1)

tell_wl_long=h2ostr.wave
m24index=where(tell_wl_long gt 2.275 and tell_wl_long lt 3.365)
h2o_long=h2ostr.trans
co2ch4_long=co2ch4str.trans
tellwl=tell_wl_long[m24index]
h2o=h2o_long[m24index]
co2ch4=co2ch4_long[m24index]


;;;Contstruct wl array and oversampled wl array from coeffs
x=dindgen(npix)
xx=(dindgen(npix*oversamp)-(oversamp/2))/oversamp

wl_soln=poly(x, wl_coeff)
;samp_index=(lindgen(npix)*oversamp)+(oversamp/2)

wl_soln_over=poly(xx, wl_coeff)

;Interpolate the lab spectrum onto the oversampled trial grid
int_lab_over_full=interpol(int_lab_depth, wl_lab, wl_soln_over)

;Save full length vectors
wl_soln_full=wl_soln
wl_soln_over_full=wl_soln_over


;Redefine those vectors to reflect only the specified range


x_select=dindgen(npix_select)+first_pix
xx_select=(dindgen(npix_model)-(oversamp/2))/oversamp+first_pix

wl_soln=poly(x_select, wl_coeff)
samp_index=(lindgen(npix_select)*oversamp)+(oversamp/2)
wl_soln_over=poly(xx_select, wl_coeff)








;;;;
;Make LSF
;;;;

lsf_tic=tic()

;populate lsf coordinates
x_lsf=rebin((dindgen(npix_lsf)-(npix_lsf/2))/oversamp, npix_lsf, npix_model)
pixel_2d=rebin(reform((((dindgen(npix_model)-(oversamp/2))/oversamp)+first_pix), 1, npix_model), npix_lsf, npix_model)
sig_arr=sigma+(sigma_lin*pixel_2d)
x_sig_2d=x_lsf/sig_arr
x_sig_3d=rebin(x_sig_2d, npix_lsf, npix_model, nparam_gh)

;Make non-normalized gh-polynomial cube (npix_lsf x npix_model x
;basis)

H_coeff=sg_hermite_coeff()
gh_cube_nonorm=dblarr(npix_lsf, npix_model, nparam_gh)
gh_cube_norm=dblarr(npix_lsf, npix_model, nparam_gh)
gh_cube=dblarr(npix_lsf, npix_model, nparam_gh)
basis_coeff=dblarr(npix_lsf, npix_model, nparam_gh)
for basis=0, nparam_gh-1 do begin
    ;make gh_polynomial for given basis
    gh_cube_nonorm[*,*,basis]=poly(x_sig_2d, H_coeff[*,basis])*exp(-0.5*(x_sig_2d^2))
    ;normalize each basis
    norm=rebin(reform(total(abs(gh_cube_nonorm[*,*,basis]), 1, /double), 1, npix_model), npix_lsf, npix_model)
    gh_cube_norm[*,*,basis]=gh_cube_nonorm[*,*,basis]/norm
    ;Adjust amplitudes according to gh parameters
    if lin_switch eq 1 then basis_coeff[*,*,basis]=gh_coeff[basis]+(gh_lin_coeff[basis]*pixel_2d) $
      else if lin_switch eq 0 then basis_coeff[*,*,basis]=replicate(gh_coeff[basis], npix_lsf, npix_model)
    gh_cube[*,*,basis]=basis_coeff[*,*,basis]*gh_cube_norm[*,*,basis]
endfor

;Sum the bases
if nparam_gh gt 1 then lsf_nonorm=total(gh_cube, 3, /double) $
  else lsf_nonorm=gh_cube[*,*,0]

;account for unrealistically negative lsfs

min_vec=min(lsf_nonorm, dimension=1)
neg_ind=where(min_vec lt 0, neg_count)
if neg_count gt 0 then begin
    neg_vec=dblarr(npix_model)
    neg_vec[neg_ind]=min_vec[neg_ind]

    neg_arr=rebin(reform(neg_vec, 1, npix_model), npix_lsf, npix_model)

    lsf_nonorm_copy=lsf_nonorm
    lsf_nonorm=lsf_nonorm-neg_arr
endif

;normalize lsf
lsf_norm=rebin(reform(total(lsf_nonorm, 1, /double), 1, npix_model), npix_lsf, npix_model)
lsf=lsf_nonorm/lsf_norm
if lin_switch eq 0 then lsf=lsf[*,0]


;;;;
;Fix LSF center of mass
;;;;
lsf_com=rebin(reform(total(lsf*x_sig_2d, 1, /double), 1, npix_model), npix_lsf, npix_model)
x_sig_2d_com=x_sig_2d+lsf_com
x_sig_3d_com=rebin(x_sig_2d_com, npix_lsf, npix_model, nparam_gh)

;;;Remake LSF
;Make non-normalized gh-polynomial cube (npix_lsf x npix_model x
;basis)

gh_cube_nonorm_com=dblarr(npix_lsf, npix_model, nparam_gh)
gh_cube_norm_com=dblarr(npix_lsf, npix_model, nparam_gh)
gh_cube_com=dblarr(npix_lsf, npix_model, nparam_gh)
for basis=0, nparam_gh-1 do begin
    ;make gh_polynomial for given basis
    gh_cube_nonorm_com[*,*,basis]=poly(x_sig_2d_com, H_coeff[*,basis])*exp(-0.5*(x_sig_2d_com^2))
    ;normalize each basis
    norm_com=rebin(reform(total(abs(gh_cube_nonorm_com[*,*,basis]), 1, /double), 1, npix_model), npix_lsf, npix_model)
    gh_cube_norm_com[*,*,basis]=gh_cube_nonorm_com[*,*,basis]/norm_com
    
    gh_cube_com[*,*,basis]=basis_coeff[*,*,basis]*gh_cube_norm_com[*,*,basis]
endfor

;Sum the bases
if nparam_gh gt 1 then lsf_nonorm_com=total(gh_cube_com, 3, /double) $
  else lsf_nonorm_com=gh_cube_com[*,*,0]

;account for unrealistically negative lsfs

min_vec_com=min(lsf_nonorm_com, dimension=1)
neg_ind_com=where(min_vec_com lt 0, neg_count_com)
if neg_count_com gt 0 then begin
    neg_vec_com=dblarr(npix_model)
    neg_vec_com[neg_ind_com]=min_vec_com[neg_ind_com]

    neg_arr_com=rebin(reform(neg_vec_com, 1, npix_model), npix_lsf, npix_model)

    lsf_nonorm_com_copy=lsf_nonorm_com
    lsf_nonorm_com=lsf_nonorm_com-neg_arr_com
endif

;normalize lsf
lsf_norm_com=rebin(reform(total(lsf_nonorm_com, 1, /double), 1, npix_model), npix_lsf, npix_model)
lsf_com=lsf_nonorm_com/lsf_norm_com
if lin_switch eq 0 then lsf_com=lsf_com[*,0]

;make copy of old lsf 
lsf_copy=lsf
lsf=lsf_com


lsf_t=toc(lsf_tic)
print, "LSF construction took : ", lsf_t, " seconds."


;Now fit each observation

;for visit=0, n_ABobj-1 do begin
for visit=7, 7 do begin

;Define Template
;temp=ABtemp_rough
bigc=299792458d0 ;m/s

;get header
struct=mrdfits(ABobj_file[visit], 1)
head=struct.header
;get bcv correction
bcvcorr_IRCS, head, params
bcv=params[0]

;Try using the 12CO line list
linefilename='12CO_1st_overtone_rest_wl.dat'
linefile=modelpath+linefilename

;read in the central wls
readcol, linefile, COwl, format='D'
COwl_full=COwl
;COwl=COwl_full[3:5]
nlines=n_elements(COwl)

;shift co lines to star frame
;also gj273 apparently has 18.22 km/s rv
gjrv_corr=-18220d0
COwl_shift=COwl*(1+((gjrv_corr+(bcv*1000d0))/bigc))


;read in observation
    int_obs=ABspec_arr[visit,*]
    mjd=ABmjd_arr[visit]

;FIX read this stuff into an array along with the observations and
;flats above
    tempstr=mrdfits(ABobj_file[visit], 1)
    error_nonorm=reverse(tempstr.sigma)
    tempspec=reverse(tempstr.spectrum)
    tempnorm=continuum_fit(dindgen(npix), tempspec, low_rej=low_reject, high_rej=high_reject)
    err=error_nonorm/tempnorm
;don't count first and last few
    err[0:20]=1d8
    err[-21:-1]=1d8




;directly compare them using AMOEBA
;;;;;;;;;;;;;;;;;;;;;;
;;;  Run Modeling Function   ;;;;
;;;;;;;;;;;;;;;;;;;;;;

clock=tic()

;Define inputs to modeling function
ftol=1d-11

rv_guess=[0d0]
;tellrv_guess=[0d0]
co2ch4_depth_guess=[0.3d0]
h2o_depth_guess=[0.3d0]
stell_depth_guess=replicate(0.3d0, nlines)
stell_width_guess=replicate(6.7d-5,nlines) 


rv_scale=[100d0]
;tellrv_scale=[1d0]
co2ch4_depth_scale=[0.2d0]
h2o_depth_scale=[0.2d0]
stell_depth_scale=replicate(0.2d0, nlines)
stell_width_scale=replicate(1d-6, nlines)


guess=[rv_guess, co2ch4_depth_guess, h2o_depth_guess, stell_depth_guess, stell_width_guess]
scale=[rv_scale, co2ch4_depth_scale, h2o_depth_scale, stell_depth_scale, stell_width_scale]
;scale=[rv_scale, depth_scale]


if visualize eq 1 then begin
    set_plot, 'x'
    window, 0, xsize=1000, ysize=400
endif
;!p.multi=[0,1,4]
!p.multi=0
;;;Run Minimization scheme

;read in last guess
;guess=[-2799.7870,-1081.8196,      0.99741826,      0.20462968,      0.12492732,      0.40048046,   8.1290099d-05,   0.00011054021,   2.2794023d-05]
;r=rvfit(guess)
;stop

if min_type eq 'amoeba' then begin
    r=amoeba3(ftol, scale=scale, p0=guess,function_name='rvfit', $
              function_value=fval, nmax=150000L)
    chi2=fval[0]
endif else if min_type eq 'mpfit' then begin
    r=mpfit('rvfit', guess, bestnorm=chi2, ftol=ftol, /quiet)
endif


!p.multi=0

;stop clock
process_time=toc(clock)


;OUTPUT RESULTS

model_id=strjoin([object, epoch,'AB'+strtrim(visit,2), strtrim(first_pix,2), strtrim(first_pix+npix_select-1,2), 'lines'],'_')

if n_elements(r) eq 1 then begin
    print, "Minimization scheme failed to converge"
    result_vec=last_guess
    rvfinal=last_guess[0]
    chi2=-1
    
endif else begin
    
    print, "---------------------------------"
    print, "Model ID: ", model_id
    print, "---------------------------------"
    
    print, chi2
    
    result_vec=r
    rvfinal=r[0]
endelse

print, "---------------------------------"
print, "Processing took ", process_time, " seconds."



outfile=model_id+".fits"

outfinal=outputpath+outfile

if file_test(outfinal) then begin
    fits_info, outfinal, n_ext=prev_ext, /silent
    extension=prev_ext+1
endif else extension=1


output_str={OBJ:object, $
            EPOCH:epoch, $
            TRACE:trace, $
            CALIB_FILE:calib_file, $
            CALIB_EXT:calib_ext, $
            WL_COEFF:wl_coeff, $
            GH0_COEFF:gh_coeff_out, $
            GH1_COEFF:gh_lin_coeff_out, $
            OTHER_COEFF:other, $
            OVERSAMP:oversamp, $
            FTOL:ftol, $
            NPIX_TRIM_START:npix_trim_start, $
            NPIX_TRIM_END:npix_trim_end, $
            FIRST_PIX:first_pix, $
            NPIX_SELECT:npix_select, $
            TIME:process_time, $
            MIN_TYPE:min_type, $
            EXT:extension, $
            CHI2:chi2, $
            MJD:mjd, $
            DELTA_RV:rvfinal, $
            RESULT_VEC:result_vec $
           }

stop
mwrfits, output_str, outfinal
 

endfor

end
