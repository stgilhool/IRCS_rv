pro make_1Dspec, epoch, object, mask=mask, no_flats=no_flats

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
outpath=epochpath+'final_spectra/'

if file_test(outpath, /directory) eq 0 then file_mkdir, outpath

;No mask for now, so
if n_elements(mask) eq 0 then mask=replicate(0, 1024)

;Start with reading in object frames
obj_list_file=epochpath+object+'_'+epoch+'_'+'1Dspec.list'

readcol, obj_list_file, obj_file_nameonly, format='A'
obj_file=epochpath+obj_file_nameonly
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
    
    obj_spec=mrdfits(obj_file[f], 0, obj_head, /dscale)*gain
    obj_spec_1D=total(obj_spec, 2, /double)

    exptime=sxpar(obj_head,'EXPTIME')
    ndr=sxpar(obj_head,'NDR')
    
    obj_var=obj_spec + (ndr*read_noise^2) + (exptime*dark_current)
    obj_var_1D=total(obj_var, 2, /double)
    obj_err=sqrt(obj_var_1D)

    obj_outstruct={header:obj_head, spectrum:obj_spec_1D, sigma:obj_err, mask:mask}

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



