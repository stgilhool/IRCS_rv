pro NH3_maketable, nparam_wl, nparam_gh, nparam_gh_lin, nparam_other

nwl=strtrim(nparam_wl, 2)
ngh0=strtrim(nparam_gh, 2)
ngh1=strtrim(nparam_gh_lin, 2)
nother=strtrim(nparam_other, 2)

resultpath='/home/stgilhool/RV_projects/IRCS_rv/cal_results/Oct28/'
nameroot='model_'
nameext='.fits'

filename=resultpath+nameroot+nwl+'_'+ngh0+'_'+ngh1+'_'+nother+nameext

fits_info, filename, /silent, n_ext=last_ext
r=mrdfits(filename, last_ext)

;stop

;wl_res=strtrim(r.wl_result, 2)
;gh0_res=strtrim(r.gh0_result,2)
;gh1_res=strtrim(r.gh1_result,2)
;other_res=strtrim(r.other_result,2)

wl_res=string(r.wl_result, format='(g10.4)')
gh0_res=string(r.gh0_result, format='(g10.4)')
gh1_res=string(r.gh1_result, format='(g10.4)')
other_res=string(r.other_result, format='(g10.4)')


d='$'
q='&'
b=' '

print, "Order"+q+"Wavelength coefficients"+q+"G-H basis coefficients"+q+"G-H linear shift"+q+"Scaling parameters\\"
for i=0, max([nparam_wl, nparam_gh, nparam_gh_lin, nparam_other])-1 do begin
    num=strtrim(i,2)
    if i lt nparam_wl and i lt nparam_gh and i lt nparam_gh_lin and i lt nparam_other then begin
      print, d+num+d+q+d+wl_res[i]+d+q+d+gh0_res[i]+d+q+d+gh1_res[i]+d+q+d+other_res[i]+d+'\\'
  endif else if i ge nparam_wl and i lt nparam_gh then begin
      print, d+num+d+q+b+q+d+gh0_res[i]+d+q+d+gh1_res[i]+d+q+d+other_res[i]+d+'\\'
  endif else if i ge nparam_wl and i ge nparam_gh then begin
      print, b+q+b+q+b+q+b+q+d+other_res[i]+d+'\\'
  endif

endfor

end
