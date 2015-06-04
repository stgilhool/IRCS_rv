pro show, exp


                                ;Set up big dots
phi=findgen(32)*(!PI*2/32.)
phi = [ phi, phi(0) ]
usersym, cos(phi), sin(phi), /fill




ttestfile2='tellurictest.fits'
t_ext=1

stestfile='penaltytest.fits'
s_ext=30


;;----------------------


trv1='GJ273_18Jan2011_AB_0_127_ttest2_r1.fits'
trv2='GJ273_18Jan2011_AB_0_127_ttest2_r2.fits'

srv1='GJ273_18Jan2011_AB_0_127_smooth30.fits'
srv2='GJ273_18Jan2011_AB_0_127_smooth30_r2.fits'

orv2='GJ273_18Jan2011_AB_0_127_run2sign.fits'


;;--------------------


trv1spec='GJ273_18Jan2011_AB_0_127_ttest2_r1_spec.fits'
trv2spec='GJ273_18Jan2011_AB_0_127_ttest2_r2_spec.fits'

srv1spec='GJ273_18Jan2011_AB_0_127_smooth30_spec.fits'
srv2spec='GJ273_18Jan2011_AB_0_127_smooth30_r2_spec.fits'

orv2spec='GJ273_18Jan2011_AB_0_127_run2sign_spec.fits'


;;-------------------


print, '0 Result of first template, 2 iterations RV fit'
print, '1 Result of fitting for a better stellar template'
print, '2 Result of 1 iteration with new stellar template (RV fixed)'
print, '3 Result of 2nd iteration with new stellar template (all parameters free)'
print, '4 Result of fitting for a better telluric template'
print, '5 Result of 1 iteration with new stellar and telluric (RV fixed)'
print, '6 Result of 2nd iteration with new stellar and telluric (all parameters free)'

;exp=1 ;Just for exposure 0

guess=dblarr(22, 7)
result=dblarr(22, 7)
model=dblarr(128, 7)
obs=dblarr(128, 7)
err=dblarr(128, 7)
wl_grid=dblarr(128, 7)
chi2=dblarr(7)


label=['RV from Template 1', 'Shifting the Template', 'Run 1 with new template', $
       'Run 2 with new template', 'Shifting Telluric', 'Run 1 with new telluric', 'Run 2 with new telluric']
list=[orv2, stestfile, srv1, srv2, ttestfile2, trv1, trv2]
list2=[orv2spec, srv1spec, srv1spec, srv2spec, trv1spec, trv1spec, trv2spec]
combined_index=[0,1,2,3,4,-1,6,-1,-1,-1,10,11,12,13]
test_index=[0,1,2,3,4,6,10,11,12,13]

foreach entry, list, index do begin
    if index ne 1 and index ne 4 then begin
        str=mrdfits(entry, exp+1)
        guess[*,index]=str.guess
        result[*,index]=str.result
        ;chi2[index]=str.chi2
        strspec=mrdfits(list2[index], exp+1)
        print, strspec.params eq str.result
        ;if strspec.params ne str.result then begin
        ;    print, strspec.params
        ;     print, str.result
        ;    stop
        ;endif
        ;if strspec.chi2 ne str.chi2 then begin
        ;    print, strspec.chi2
        ;    print, str.chi2
        ;    stop
        ;endif
        model[*,index]=strspec.model
        obs[*,index]=strspec.obs
        err[*,index]=strspec.error
        chi2[index]=total(((obs[*,index]-model[*,index])/err[*,index])^2, /double)
    endif else if index eq 1 then begin
        if combined_index[exp] lt 0 then begin
            label[index]='No data'
        endif else begin
            num=where(test_index eq exp)
            str=mrdfits(entry, s_ext)
            obs[*,index]=str.obs_arr[*,num]
            model[*,index]=str.model_arr[*,num]
            err[*,index]=str.err_arr[*,num]
            chi2[index]=total(((obs[*,index]-model[*,index])/err[*,index])^2, /double)   
        endelse
    endif else if index eq 4 then begin
        if combined_index[exp] lt 0 then label[index]='No data' else begin
            num=where(test_index eq exp)
            str=mrdfits(entry, t_ext)
            obs[*,index]=str.obs_arr[*,num]
            model[*,index]=str.model_arr[*,num]
            err[*,index]=str.err_arr[*,num]
            guess[0,index]=str.delta_rv_guess[num]
            result[0,index]=str.delta_rv_final[num]
            guess[1,index]=str.h2o_depth_guess[num]
            result[1,index]=str.h2o_depth_final[num]
            guess[2,index]=str.co2ch4_depth_guess[num]
            result[2,index]=str.co2ch4_depth_final[num]
            chi2[index]=total(((obs[*,index]-model[*,index])/err[*,index])^2, /double)
        endelse
    endif
endforeach


!p.multi=[0,2,7]


for run=0,6 do begin
    titlestr=label[run] + "  | chi2: " + strtrim(chi2[run], 2)
    plot, obs[*,run], /xs, yr=[0,1.2], title=titlestr, charsize=2.0
    oplot, model[*,run], ps=8, color=200, symsize=0.3
    plot, obs[*,run]-model[*,run], ps=3, yr=[-0.1,0.1], /xs
endfor

!p.multi=0
stop


end
