pro make_1Dspec_test, epoch, object, mask=mask, no_flats=no_flats


;Error related info
gain=3.8d0 ;electrons/ADU
read_noise=68d0 ;electrons/NDR
dark_current=0.05d0 ;electrons/second

;set file paths
rootpath='/home/stgilhool/RV_projects/IRCS_rv/data/'
epochpath=rootpath+'epoch/'+epoch+'/'
objectpath=rootpath+object+'/'+epoch+'/obs/'
flatpath=rootpath+'cal_files/'+epoch+'/'
outpath=epochpath+'final_spectra/'

if file_test(outpath, /directory) eq 0 then file_mkdir, outpath

;No mask for now, so
if n_elements(mask) eq 0 then mask=replicate(0, 1024)

;Start with reading in object frames
obj_list_file=epochpath+object+'_'+epoch+'_'+'1Dspec.list'

readcol, obj_list_file, obj_file, format='A'
n_obj=n_elements(obj_file)

;Loop through each object file and process it
for f=0, n_obj-1 do begin

    outstr=strsplit(obj_file[f], '.', /extract)
    outstr[1]='1D'
    ;outfile=outpath+strjoin(outstr, '.')
outfile=outpath+strjoin([strjoin([outstr[0],outstr[1]], '_'),outstr[2]], '.')
        

    obj_spec=mrdfits(obj_file[f], 0, obj_head, /dscale)*gain
    obj_spec_1D=total(obj_spec, 2, /double)

    exptime=sxpar(obj_head,'EXPTIME')
    ndr=sxpar(obj_head,'NDR')
    
    obj_var=obj_spec + (ndr*read_noise^2) + (exptime*dark_current)
    obj_var_1D=total(obj_var, 2, /double)
    obj_err=sqrt(obj_var_1D)

    outstruct={header:obj_head, spectrum:obj_spec_1D, sigma:obj_err, mask:mask}

    mwrfits, outstruct, outfile, /create

endfor

end

