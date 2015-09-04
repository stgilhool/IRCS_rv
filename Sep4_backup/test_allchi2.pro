; IDL Version 8.2.2 (linux x86_64 m64)
; Journal File for stgilhool@iroquois.physics.upenn.edu
; Working directory: /RAID/home/stgilhool/RV_projects/IRCS_rv/IDL_pros
; Date: Thu Mar 12 07:57:03 2015
pro test_allchi2, extension, test=test


;;;SET TEST KEYWORD 
if n_elements(test) eq 0 then begin
    file='/home/stgilhool/RV_projects/IRCS_rv/data/h2o_test/test1/h2otest.fits'
endif else begin
    file='/home/stgilhool/RV_projects/IRCS_rv/data/h2o_test/test'+ $
      strtrim(test,2)+$
      '/h2otest.fits'
endelse

;;;TEST AND SET EXTENSION ARGUMENT
fits_info, file, n_ext=n_ext, /silent
if n_params() eq 0 then extension = 1
if extension gt n_ext then message, "Extension does not exist"
if extension lt 1 then message, "Extension 1 is lowest number"


;;;READ IN CORRECT EXTENSION
a=mrdfits(file, extension)

;;;CALCULATE CHI2 FOR EACH EXPOSURE
n_exp=n_elements(a.h2o_depth_guess)
chi2_arr=dblarr(n_exp)
for i=0, n_exp-1 do begin
  chi2_arr[i]=total(((a.model_arr[*,i]-a.obs_arr[*,i])/a.err_arr[*,i])^2, /double)
endfor

;;;READ IN RV RUN INPUT AS WELL
varray=[0,1,2,3,4,6,10,11,12,13]

if n_elements(varray) ne n_exp then message, "Varray must match"

chi2rv_arr=dblarr(n_exp)

rfile='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_127_t7newpen_r2free_spec.fits'

foreach visit, varray, index do begin

    r=mrdfits(rfile, visit+1)
    
    chi2rv_arr[index]=total(((r.model-r.obs)/r.error)^2, /double)

endforeach

print, chi2_arr
print, chi2rv_arr

print, "The following should be positive if all is well:"
print, chi2rv_arr-chi2_arr


stop



end

