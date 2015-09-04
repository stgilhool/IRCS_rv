pro ircsrv_readobs, object=object, epoch=epoch, trace=trace, obs_list=obs_list, spectra=spectra, errors=errors, mjds=mjds, headers=headers

if n_elements(object) eq 0 then object='GJ273'
if n_elements(epoch) eq 0 then epoch='18Jan2011'
if n_elements(trace) eq 0 then trace='AB'



;set file paths
rootpath='/home/stgilhool/RV_projects/IRCS_rv/data/'
epochpath=rootpath+'epoch/'+epoch+'/'
objectpath=epochpath+'final_spectra/'
flatpath=epochpath+'final_spectra/'
calibpath=epochpath+'calib_results/'
outpath=epochpath+'temp_results/'
modelpath=rootpath+'supplemental/'
outputpath=epochpath+'rv_results/'


;Define object file names and full path file names

if n_elements(obs_list) eq 0 then begin
    obs_listname=strjoin([object, epoch, trace], '_')+'.list'
    
    obs_list=objectpath+obs_listname
endif



;Read in object filenames and then spectra
readcol, obs_list, obs_filename, format='A'


obs_file=objectpath+obs_filename


n_obs=n_elements(obs_file)



;Initialize arrays of spectra
temp_str=mrdfits(obs_file[0], 1)
temp_spectrum=temp_str.spectrum
npix=n_elements(temp_spectrum)
spec_arr=dblarr(n_obs, npix)





mjd_arr=dblarr(n_obs)


;define error arrays
err_arr=dblarr(n_obs, npix)




;Set continuum normalization params
low_reject=0.53
high_reject=3.0


;Loop through each object file and flip and normalize it
for f=0, n_obs-1 do begin
    obs_spec=mrdfits(obs_file[f], 1)


    spec_backwards=obs_spec.spectrum


    ;Flip the spectra
    spec=reverse(spec_backwards)


        ;Read in and flip errors
    err=reverse(obs_spec.sigma)


    ;Get MJD
    
    head=obs_spec.header
    if f eq 0 then begin
        head_save=replicate({head:head}, n_obs)
    endif else head_save[f].head=head
    

    ;Normalize the model spectrum and the observed spectrum
    norm=continuum_fit(dindgen(npix), spec, low_rej=low_reject, high_rej=high_reject)

    
    ;Populate arrays
    ;spectra
    spec_arr[f,*]=spec/norm

    ;errors
    err_arr[f,*]=err/norm

    ;MJD
    mjd_arr[f]=sxpar(head, 'MJD')


endfor

spectra=spec_arr
errors=err_arr
mjds=mjd_arr
headers=head_save

end
