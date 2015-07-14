pro make_1Dspec_jul9, epoch, object, mask=mask, no_flats=no_flats

;keyword
if n_elements(no_flats) eq 0 then no_flats=1

;Error related info
gain=3.8d0 ;electrons/ADU
read_noise=68d0 ;electrons/NDR
dark_current=0.05d0 ;electrons/second
flat_time=0.5 ;seconds

;set file paths
rootpath='/home/stgilhool/RV_projects/IRCS_rv/data/'
epochpath=rootpath+'epoch/'+epoch+'/'
objectpath=rootpath+object+'/'+epoch+'/obs/'
flatpath=rootpath+'cal_files/'+epoch+'/'
outpath=objectpath+'final_spectra/'

if file_test(outpath, /directory) eq 0 then file_mkdir, outpath

;No mask for now, so
if n_elements(mask) eq 0 then mask=replicate(0, 1024)

;maskf=objectpath+'bpmask.0001.fits'
;ma=mrdfits(maskf,0)
;mask=total(ma,2)
;if n_elements(mask) ne 1024 then message, "summed over wrong dimension in mask"
;flag=where(mask gt 0, fcount)
;if fcount gt 0 then mask[flag]=1

;Start with reading in object frames
obj_list_file=objectpath+object+'_'+epoch+'_'+'1Dspec.list'

readcol, obj_list_file, obj_file_nameonly, format='A'
obj_file=objectpath+obj_file_nameonly

sig_list_file=objectpath+object+'_'+epoch+'_'+'1Dsig.list'
readcol, sig_list_file, sig_file_nameonly, format='A'
sig_file=objectpath+sig_file_nameonly

bp_list_file=objectpath+object+'_'+epoch+'_'+'1Dbp.list'
readcol, bp_list_file, bp_file_nameonly, format='A'
bp_file=objectpath+bp_file_nameonly


n_obj=n_elements(obj_file)

;Loop through each object file and process it
for f=0, n_obj-1 do begin

                                ;THIS PART RELIES ON THE ENDING BEING
                                ;.0001.fits WITH NO OTHER
                                ;PERIODS. PROBABLY SHOULD FIND A MORE
                                ;ELEGANT WAY
    outstr=strsplit(obj_file_nameonly[f], '.', /extract)
    outstr[1]='1D'
    obj_outfile=outpath+strjoin([strjoin([outstr[0],outstr[1]], '_'),outstr[2]], '.')
    
    print,'reading objs for file '+strtrim(f,2)
    obj_spec=mrdfits(obj_file[f], 0, obj_head, /dscale)*gain
    obj_spec_1D=total(obj_spec, 2, /double)

    exptime=sxpar(obj_head,'EXPTIME')
    ndr=sxpar(obj_head,'NDR')
    
    ;THIS IS CHANGED AS OF JUL8,2015.  READNOISE FOR A GIVEN FRAME IS
    ;GIVEN BY RN/SQRT(NDR)

;     obj_var=obj_spec + ((read_noise^2)/ndr) + (exptime*dark_current)
;     obj_var_1D=total(obj_var, 2, /double)
;     obj_err=sqrt(obj_var_1D)
    print, 'readin err for file '+strtrim(f,2)
    obj_err_2d=mrdfits(sig_file[f], 0, /dscale)*sqrt(gain)  ;/sqrt(gain) this was to convert wrong sigma
    var_2d=obj_err_2d^2
    obj_var=total(var_2d, 2, /double)
    obj_err=sqrt(obj_var)


    ;Mask bp from bpmask
    print,'reading bpmask for file '+strtrim(f,2)
    bp_2d=mrdfits(bp_file[f], 0)
    bp_mask=total(bp_2d,2)
    mask=bp_mask
    bp_flag=where(bp_mask gt 0, bpfcount)

    ;CHECK DIMS
    help, obj_spec
    help, obj_err_2d
    help, bp_2d


    if bpfcount gt 0 then begin
        ;mask adjacent pixels NEVERMIND, let's do it after cr masking
        ; bp_flag_up=bp_flag+1
;         bp_flag_down=bp_flag-1
        
;         bp_flag_cat=[bp_flag, bp_flag_up, bp_flag_down]
;         bp_flag_cat_sort=bp_flag_cat[sort(bp_flag_cat)]
;         bp_flag_cat_uniq=bp_flag_cat_sort[uniq(bp_flag_cat_sort)]
;         bp_flag_all=bp_flag_cat_uniq[where(bp_flag_cat_uniq ge 0 and bp_flag_cat_uniq lt 1024)]

        mask[bp_flag]=1
    endif

    


    
    ;Mask cosmic rays
    spec_med=median(obj_spec, dimension=2)
    spec_mean=mean(obj_spec, dimension=2)
    spec_var=mean((obj_spec-rebin(spec_mean, n_elements(obj_spec[*,0]), n_elements(obj_spec[0,*])))^2, dimension=2)
    spec_sig=sqrt(spec_var)



    nerrtol=3
    for col=0,n_elements(spec_sig)-1 do begin
        ;spec_fit_coeff=robust_poly_fit(lindgen(n_elements(obj_spec[0,*])), obj_spec[col,*], 2, spec_fit, fit_sig)
;        spec_fit_coeff=gaussfit(lindgen(n_elements(obj_spec[0,*]))-(n_elements(obj_spec[0,*)/2), obj_spec[col,*], 2, spec_fit, fit_sig)
        ;spec_fit=gaussfit(lindgen(n_elements(obj_spec[0,*])),
        ;obj_spec[col,*], coeff, nterms=4,yerror=yerror)
        spec_fit=mpfitpeak(lindgen(n_elements(obj_spec[0,*])), obj_spec[col,*], coeff, nterms=4,yerror=yerror, error=obj_err_2d, /lorentzian)
        err_tol=nerrtol*yerror
        for row=0,n_elements(obj_spec[0,*])-1 do begin
;            if abs(obj_spec[col,row]-spec_med[col]) gt spec_sig[col] then begin
;            if abs(obj_spec[col,row]-spec_fit[col]) gt 3*fit_sig then
;            begin 
            
            ;if f eq 7 then err_tol=3*yerror
            if abs(spec_fit[row]-obj_spec[col,row]) gt err_tol then begin 
               plot, obj_spec[col,*], title=strtrim(f,2)+'_'+strtrim(col,2)+', '+strtrim(row,2)
               oplot, spec_fit, color=200
                rep=1
                stop
                if rep eq 1 then begin
                    if mask[col] ne 1 then mask[col]=1
                    xx=lindgen(n_elements(obj_spec[0,*]))
                    x=xx[where(xx ne row)]
                    new_spec=obj_spec[col,x]
                    new_spec_fit=gaussfit(x, new_spec, nterms=4)
                    newnew_spec_fit=interpol(new_spec_fit, x, xx)
                    obj_spec_1d[col]=total(newnew_spec_fit)
                endif
                
            endif
        endfor
    endfor

    mask_ind=where(mask ne 0, mmcount)
    if mmcount gt 0 then begin
        ;mask adjacent pixels to masked ones
         bp_flag_up=mask_ind+1
         bp_flag_down=mask_ind-1
        
         bp_flag_cat=[bp_flag, bp_flag_up, bp_flag_down]
         bp_flag_cat_sort=bp_flag_cat[sort(bp_flag_cat)]
         bp_flag_cat_uniq=bp_flag_cat_sort[uniq(bp_flag_cat_sort)]
         bp_flag_all=bp_flag_cat_uniq[where(bp_flag_cat_uniq ge 0 and bp_flag_cat_uniq lt 1024)]
         mask[bp_flag_all]=1
         ;dial up error in masked pixels
        obj_err[bp_flag_all]=1d10
    endif

    obj_outstruct={header:obj_head, spectrum:obj_spec_1D, sigma:obj_err, mask:mask}

    print, 'writing structure for file '+strtrim(f,2)
    mwrfits, obj_outstruct, obj_outfile, /create

endfor


;Only do this part if no_flats is set to 0

for exp=1,1 do begin ;this is for when there are more extracted apertures for the flats

if no_flats eq 1 then break

    for pattern=0, 1 do begin
        
        if pattern eq 0 then pat='_AB' else if pattern eq 1 then pat='_BA'
        
        list_ext=pat+strtrim(exp,2)+'.list'

        calONcellIN_list=epochpath+'calONcellIN_'+object+list_ext
        calOFFcellIN_list=epochpath+'calOFFcellIN_'+object+list_ext
        calONcellOUT_list=epochpath+'calONcellOUT_'+object+list_ext
        calOFFcellOUT_list=epochpath+'calOFFcellOUT_'+object+list_ext

;read in file names to string arrays
        readcol, calONcellIN_list, calONcellIN_file, format='A'
        readcol, calOFFcellIN_list, calOFFcellIN_file, format='A'
        readcol, calONcellOUT_list, calONcellOUT_file, format='A'
        readcol, calOFFcellOUT_list, calOFFcellOUT_file, format='A'
;add path to epochpath
        calONcellIN_file=epochpath+calONcellIN_file
        calOFFcellIN_file=epochpath+calOFFcellIN_file
        calONcellOUT_file=epochpath+calONcellOUT_file
        calOFFcellOUT_file=epochpath+calOFFcellOUT_file

;initialize arrays to store spectra and error data
        nflats=n_elements(calONcellIN_file)
        temp=mrdfits(calONcellIN_file[0], 0)
        ncols=n_elements(temp[0,*])
        nrows=n_elements(temp[*,0])

        calONcellIN_arr=dblarr(nflats, nrows, ncols)
        calOFFcellIN_arr=dblarr(nflats, nrows, ncols)
        calONcellOUT_arr=dblarr(nflats, nrows, ncols)
        calOFFcellOUT_arr=dblarr(nflats, nrows, ncols)

        calONcellIN_var=dblarr(nflats, nrows, ncols)
        calOFFcellIN_var=dblarr(nflats, nrows, ncols)
        calONcellOUT_var=dblarr(nflats, nrows, ncols)
        calOFFcellOUT_var=dblarr(nflats, nrows, ncols)


;Read in spectra in number of photo-electrons
        for file=0,nflats-1 do begin
    
            calONcellIN_arr[file, *, *]=mrdfits(calONcellIN_file[file], 0, /dscale)*gain
            calOFFcellIN_arr[file, *, *]=mrdfits(calOFFcellIN_file[file], 0, /dscale)*gain
            calONcellOUT_arr[file, *, *]=mrdfits(calONcellOUT_file[file], 0, /dscale)*gain
            calOFFcellOUT_arr[file, *, *]=mrdfits(calOFFcellOUT_file[file], 0, /dscale)*gain
            
        endfor

;variance at all pixels in each pixel in each flat
        calONcellIN_var=calONcellIN_arr + (read_noise^2) + (dark_current*flat_time)
        calOFFcellIN_var=calOFFcellIN_arr + (read_noise^2) + (dark_current*flat_time)
        calONcellOUT_var=calONcellOUT_arr + (read_noise^2) + (dark_current*flat_time)
        calOFFcellOUT_var=calOFFcellOUT_arr + (read_noise^2) + (dark_current*flat_time)
        
;Sum the columns and spectra to get total counts per pixel (DISABLED
;BECAUSE OF COSMIC RAYS)
        ;calONcellIN_tot=total(total(calONcellIN_arr, 3), 1)
        ;calOFFcellIN_tot=total(total(calOFFcellIN_arr, 3), 1)
        ;calONcellOUT_tot=total(total(calONcellOUT_arr, 3), 1)
        ;calOFFcellOUT_tot=total(total(calOFFcellOUT_arr, 3), 1)

        ;calONcellIN_tot_var=total(total(calONcellIN_var, 3),1)
        ;calOFFcellIN_tot_var=total(total(calOFFcellIN_var, 3),1)
        ;calONcellOUT_tot_var=total(total(calONcellOUT_var, 3),1)
        ;calOFFcellOUT_tot_var=total(total(calOFFcellOUT_var, 3),1)

;;;Try making 1D spectrum by throwing out the highest and lowest
;;;files at each pixel to counter cosmic rays, and then summing
        ;Make each file a 1D spectrum

         calONcellIN_1D_arr=total(calONcellIN_arr, 3)
         calOFFcellIN_1D_arr=total(calOFFcellIN_arr, 3)
         calONcellOUT_1D_arr=total(calONcellOUT_arr, 3)
         calOFFcellOUT_1D_arr=total(calOFFcellOUT_arr, 3)
        
         calONcellIN_1D_var=total(calONcellIN_var, 3)
         calOFFcellIN_1D_var=total(calOFFcellIN_var, 3)
         calONcellOUT_1D_var=total(calONcellOUT_var, 3)
         calOFFcellOUT_1D_var=total(calOFFcellOUT_var, 3)

                                ;TAKE OUT STARTING HERE TO REVERT TO
                                ;THE CODE THAT REMOVES JUST 1 MIN/MAX
         dims=size(calONcellIN_1D_arr, /dimensions)
         ONINzero=replicate(1d0, dims)
         OFFINzero=replicate(1d0, dims)
         ONOUTzero=replicate(1d0, dims)
         OFFOUTzero=replicate(1d0, dims)

         for iter=0,1 do begin
         
         ;Find which files are max and min at each pixel
             maxONIN=max(calONcellIN_1D_arr, maxONINind, dimension=1, min=minONIN, subscript_min=minONINind, /nan)
             maxOFFIN=max(calOFFcellIN_1D_arr, maxOFFINind, dimension=1, min=minOFFIN, subscript_min=minOFFINind, /nan)
             maxOFFOUT=max(calOFFcellOUT_1D_arr, maxOFFOUTind, dimension=1, min=minOFFOUT, subscript_min=minOFFOUTind, /nan)
             maxONOUT=max(calONcellOUT_1D_arr, maxONOUTind, dimension=1, min=minONOUT, subscript_min=minONOUTind, /nan)
         

             ;NAN the corresponding entries
             calONcellIN_1D_arr[[maxONINind,minONINind]]=!values.f_nan
             calOFFcellIN_1D_arr[[maxOFFINind,minOFFINind]]=!values.f_nan
             calONcellOUT_1D_arr[[maxONOUTind,minONOUTind]]=!values.f_nan
             calOFFcellOUT_1D_arr[[maxOFFOUTind,minOFFOUTind]]=!values.f_nan
        
             calONcellIN_1D_var[[maxONINind,minONINind]]=!values.f_nan
             calOFFcellIN_1D_var[[maxOFFINind,minOFFINind]]=!values.f_nan
             calONcellOUT_1D_var[[maxONOUTind,minONOUTind]]=!values.f_nan
             calOFFcellOUT_1D_var[[maxOFFOUTind,minOFFOUTind]]=!values.f_nan

         endfor

         ;Now zero the actual NANed elements
         ONINnan=where(finite(calONcellIN_1D_arr) eq 0)	
         OFFINnan=where(finite(calOFFcellIN_1D_arr) eq 0)	
         ONOUTnan=where(finite(calONcellOUT_1D_arr) eq 0)	
         OFFOUTnan=where(finite(calOFFcellOUT_1D_arr) eq 0)	

         calONcellIN_1D_arr[ONINnan]=0d0
         calOFFcellIN_1D_arr[OFFINnan]=0d0
         calONcellOUT_1D_arr[ONOUTnan]=0d0
         calOFFcellOUT_1D_arr[OFFOUTnan]=0d0
         
         calONcellIN_1D_var[ONINnan]=0d0
         calOFFcellIN_1D_var[OFFINnan]=0d0
         calONcellOUT_1D_var[ONOUTnan]=0d0
         calOFFcellOUT_1D_var[OFFOUTnan]=0d0


                                ;AND END HERE TO REPLACE CODE WITH BELOW

        
         ;Find which files are max and min at each pixel
         ;THIS METHOD WORKED FOR REMOVING JUST 1 MIN AND MAX PIXEL
;          maxONIN=max(calONcellIN_1D_arr, maxONINind, dimension=1, min=minONIN, subscript_min=minONINind)
;          maxOFFIN=max(calOFFcellIN_1D_arr, maxOFFINind, dimension=1, min=minOFFIN, subscript_min=minOFFINind)
;          maxOFFOUT=max(calOFFcellOUT_1D_arr, maxOFFOUTind, dimension=1, min=minOFFOUT, subscript_min=minOFFOUTind)
;          maxONOUT=max(calONcellOUT_1D_arr, maxONOUTind, dimension=1, min=minONOUT, subscript_min=minONOUTind)
         
                                ;Zero the max and min entries in both
                                ;the spectrum and variance arrays

;          calONcellIN_1D_arr[[maxONINind,minONINind]]=0d0
;          calOFFcellIN_1D_arr[[maxOFFINind,minOFFINind]]=0d0
;          calONcellOUT_1D_arr[[maxONOUTind,minONOUTind]]=0d0
;          calOFFcellOUT_1D_arr[[maxOFFOUTind,minOFFOUTind]]=0d0
        
;          calONcellIN_1D_var[[maxONINind,minONINind]]=0d0
;          calOFFcellIN_1D_var[[maxOFFINind,minOFFINind]]=0d0
;          calONcellOUT_1D_var[[maxONOUTind,minONOUTind]]=0d0
;          calOFFcellOUT_1D_var[[maxOFFOUTind,minOFFOUTind]]=0d0

         ;Consolodate all to 1D
        
         calONcellIN_tot=total(calONcellIN_1D_arr, 1)
         calOFFcellIN_tot=total(calOFFcellIN_1D_arr, 1)
         calONcellOUT_tot=total(calONcellOUT_1D_arr, 1)
         calOFFcellOUT_tot=total(calOFFcellOUT_1D_arr, 1)

         calONcellIN_tot_var=total(calONcellIN_1D_var,1)
         calOFFcellIN_tot_var=total(calOFFcellIN_1D_var,1)
         calONcellOUT_tot_var=total(calONcellOUT_1D_var,1)
         calOFFcellOUT_tot_var=total(calOFFcellOUT_1D_var,1)


;Now take on-off
        cellIN_arr=calONcellIN_tot-calOFFcellIN_tot
        cellOUT_arr=calONcellOUT_tot-calOFFcellOUT_tot

;errors
        cellIN_arr_var=calONcellIN_tot_var+calOFFcellIN_tot_var
        cellOUT_arr_var=calONcellOUT_tot_var+calOFFcellOUT_tot_var

;Finally, do gc/nc
        observed_spectrum=cellIN_arr/cellOUT_arr
;errors
        var_final=(observed_spectrum^2)*((cellIN_arr_var/(cellIN_arr^2))+(cellOUT_arr_var/(cellOUT_arr^2)))
        error=sqrt(var_final)

        ;make final structure
        NH3_outstruct={nh3_spectrum:observed_spectrum,nh3_sigma:error,nh3_mask:mask}

        ;output structure
        NH3_outfile=outpath+'NH3_'+object+'_'+epoch+pat+strtrim(exp,2)+'.fits'

        mwrfits, NH3_outstruct, NH3_outfile, /create
        
    endfor
endfor



end



