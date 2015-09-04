function templatefunc, p

common func_info, wl_soln, samp_index, wl_soln_over_full, int_lab_over, int_obs, lsf, npix_lsf, npix_model, first_pix, npix_select, oversamp, norm_factor, err, visualize, npix_trim_start, npix_trim_end, min_type, nparam_other, last_guess


;Make pre-conv observation
;Oversample the rough template
ABtemp=p
last_guess=p
ABtemp_over_select=interpol(ABtemp, wl_soln, wl_soln_over_full)

ABtemp_over=ABtemp_over_select
;ABtemp_over=replicate(1d0, n_elements(int_lab_over))
;ABtemp_over[first_pix*oversamp:(first_pix+npix_select)*oversamp-1]=ABtemp_over_select

obs_noconv=ABtemp_over*int_lab_over


;;;;
;Convolve model with lsf
;;;;

int_conv=dblarr(npix_model)
int_conv_mtx=dblarr(npix_lsf, npix_model)


for index=0,npix_model-1 do begin
    model_index=index+first_pix*oversamp
    if model_index lt (npix_lsf/2) or model_index gt (n_elements(int_lab_over)-1-(npix_lsf/2)) then $
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
int_model=int_avg[samp_index] 


;;;;
;Correct subtle variations in continuum
;;;;

if nparam_other gt 0 then int_model=int_model*norm_factor



;;;;
;;;Get chi^2
;;;;

int_obs_dim=int_obs[first_pix:first_pix+npix_select-1]
err_dim=err[first_pix:first_pix+npix_select-1]

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
    
    title_str="Template Making Process | chi2: "+strtrim(chi2,2)

    plot, wl_soln, int_obs_dim,yr=[0.1,1.1], /xs, title=title_str,xtitle="Wavelength (microns)", ytitle="Relative Instensity", charsize=2.0 
    oplot, wl_soln, int_model, color=200 

    plot, wl_soln, res, title="Residuals | RMS : ", xtitle="Wavelength (microns)", yr=[-0.1, 0.1], ps=3, /xs

    plot, wl_soln, ABtemp, /xs
    
endif



;;;;
;Report results
;;;;

amoebafunc_t=toc(amoebafunc_tic)
;print, "One full iteration of amoebafunc took : ",amoebafunc_t, " s"
;print, "LSF construction took ", (lsf_t/amoebafunc_t)*100d0, " percent of the time"
;print, "Convolution took ", (convol_t/amoebafunc_t)*100d0, " percent of the time"
;print, "Plotting took ", plot_t, " seconds and ", (plot_t/amoebafunc_t)*100d0, " percent of the time"


print, chi2
print, "----------"

; if model_func ne 0 then begin
;     print, "Stopping for user analysis of model"
;     print, "Define variable for outfilename"
;     outfilename='tempout.fits'
;     outstr = {lsf:lsf, $
;               lsf_decomp: lsf_decomp, $
;               wl_soln: wl_soln, $
;               wl_soln_over: wl_soln_over, $
;               oversamp: oversamp, $
;               int_model: int_model, $
;               norm_factor: norm_factor, $
;               int_conv: int_conv $
;               }
;     stop
;     ;IRCS_make_figs, chi2perdegree, wl_soln, lsf, oversamp, npix_lsf, int_obs, int_model, npix_trim_start, npix_trim_end, lin_switch
;     ;stop
;     mwrfits, outstr, outfilename, /create
;     print, "outfile written"
;     stop

; endif else begin

    if min_type eq 'amoeba' then return, chi2 $
    else if min_type eq 'mpfit' then return, chi1_vec $
    else message, 'incorrect min_type... cannot return a chi2 value'
;endelse



end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;MAIN BODY ;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro make_template, epoch, object, trace, visualize=visualize, first_pix=first_pix, npix_select=npix_select, min_type=min_type

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

;Now pick a fiducial stellar observation and make the temp the guess
;and all free parameters

;We're going to use AB2

;Read in calibration results
calib_file='GJ273_18Jan2011_AB1_5_7_7_13.fits'
calib_ext=13
model_file=calibpath+calib_file
model_par=mrdfits(model_file, calib_ext)


common func_info, wl_soln, samp_index, wl_soln_over_full, int_lab_over, int_obs, lsf, npix_lsf, npix_model, firstpix, npixselect, oversamp, norm_factor, err, visual, npix_trim_start, npix_trim_end, mintype, nparam_other, last_guess

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
npix_trim_start=0L
npix_trim_end=0L


;Read in LAB SPECTRUM
modelfile= modelpath+'NH3_model.dat'
readcol, modelfile, wl_lab, int_lab, format='D,D'

norm_lab=continuum_fit(wl_lab, int_lab, low_rej=low_reject, high_rej=high_reject)
int_lab_copy=int_lab
int_lab=int_lab_copy/norm_lab


;Adjust the optical depth of the lab spectrum (I think this is okay to
;do for the high-resolution spectrum even though the wl scale isn't linear)
int_lab_depth=int_lab^tau_scale

;;;Contstruct wl array and oversampled wl array from coeffs
x=dindgen(npix)
xx=(dindgen(npix*oversamp)-(oversamp/2))/oversamp

wl_soln=poly(x, wl_coeff)
;samp_index=(lindgen(npix)*oversamp)+(oversamp/2)

wl_soln_over=poly(xx, wl_coeff)

;Interpolate the lab spectrum onto the oversampled trial grid
int_lab_over=interpol(int_lab_depth, wl_lab, wl_soln_over)

;Save full length vectors
wl_soln_full=wl_soln
wl_soln_over_full=wl_soln_over
int_lab_over_full=int_lab_over

;Redefine those vectors to reflect only the specified range
npix_model=npix_select*oversamp

x_select=dindgen(npix_select)+first_pix
xx_select=(dindgen(npix_model)-(oversamp/2))/oversamp+first_pix

wl_soln=poly(x_select, wl_coeff)
samp_index=(lindgen(npix_select)*oversamp)+(oversamp/2)
wl_soln_over=poly(xx_select, wl_coeff)

ABtemp_rough_select=ABtemp_rough[first_pix:first_pix+npix_select-1]

;FIX remove this when no longer needed
;write out the rough temps
roughtempname=strjoin([object, epoch, trace, 'roughtemp'],'_')+'.dat'
roughtempselectname=strjoin([object, epoch, trace, strtrim(first_pix,2), strtrim(first_pix+npix_select-1,2), 'roughtemp'],'_')+'.dat'
roughtempout=outpath+roughtempname
roughtempselectout=outpath+roughtempselectname
openw, lun, roughtempout, /get_lun
openw, lunsel, roughtempselectout, /get_lun
for line=0, n_elements(ABtemp_rough)-1 do printf, lun, ABtemp_rough[line]
for line=0, n_elements(ABtemp_rough_select)-1 do printf, lunsel, ABtemp_rough_select[line]
close, lun
close, lunsel
free_lun, lun
free_lun, lunsel

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



;OBSERVATION
int_obs=ABspec_arr[1,*]
mjd=ABmjd_arr[1]

tempstr=mrdfits(ABobj_file[1], 1)
error_nonorm=reverse(tempstr.sigma)
tempspec=reverse(tempstr.spectrum)
tempnorm=continuum_fit(dindgen(npix), tempspec, low_rej=low_reject, high_rej=high_reject)
err=error_nonorm/tempnorm
;don't count first and last few
err[0:5]=1d8
err[-6:-1]=1d8




;directly compare them using AMOEBA
;;;;;;;;;;;;;;;;;;;;;;
;;;  Run Modeling Function   ;;;;
;;;;;;;;;;;;;;;;;;;;;;

clock=tic()

;Define inputs to modeling function
ftol=5d-8
guess=ABtemp_rough_select
scale=0.05d0*ABtemp_rough_select

if visualize eq 1 then window, 0, xsize=1000, ysize=400
!p.multi=[0,1,3]
;;;Run Minimization scheme


if min_type eq 'amoeba' then begin
    r=amoeba3(ftol, scale=scale, p0=guess,function_name='templatefunc', $
              function_value=fval, nmax=150000L)
    chi2=fval[0]
endif else if min_type eq 'mpfit' then begin
    r=mpfit('templatefunc', guess, bestnorm=chi2, ftol=ftol, /quiet)
endif


!p.multi=0

;stop clock
process_time=toc(clock)


;OUTPUT RESULTS

model_id=strjoin([object, epoch, 'AB2', strtrim(first_pix,2), strtrim(first_pix+npix_select-1,2)],'_')

if n_elements(r) eq 1 then begin
    print, "Minimization scheme failed to converge"
    r=last_guess
    chi2=-1
    
endif else begin
    
    print, "---------------------------------"
    print, "Model ID: ", model_id
    print, "---------------------------------"
    
    print, chi2
    
    
endelse

print, "---------------------------------"
print, "Processing took ", process_time, " seconds."



outfile=model_id+"_template.fits"

outfinal=outpath+outfile

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
            TEMPLATE:r $
           }


mwrfits, output_str, outfinal
 

end
