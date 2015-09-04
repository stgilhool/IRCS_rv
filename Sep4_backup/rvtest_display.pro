pro rvtest_display, first_read=first_read


;READIN FILES
file_tag='wlgh_test3'
file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_310_398_Mar24_'+file_tag+'_r0.fits'

fits_info, file, n_ext=n_trials, /silent

;INITIALIZE VECTORS

if n_elements(first_read) eq 0 then first_read=1 else first_read=0

first_read=0

if first_read eq 1 then begin 
    visit=dblarr(n_trials)
    chi2=dblarr(n_trials)
    chi2_nopen=dblarr(n_trials)
    wl_start=dblarr(n_trials)
    wl_result=dblarr(n_trials)
    sig_start=dblarr(n_trials)
    sig_result=dblarr(n_trials)
    mjd=dblarr(n_trials)
    bcv=dblarr(n_trials)
    delta_bcv=dblarr(n_trials)
    delta_rv=dblarr(n_trials)
    final_rv=dblarr(n_trials)
    status_flag=dblarr(n_trials)
    penalty_flag=dblarr(n_trials)
    chi2_flag=dblarr(n_trials)
    
    for i=0,n_trials-1 do begin
        
        f=mrdfits(file,i+1)
        
        visit[i]=f.visit
        chi2[i]=f.chi2
        chi2_nopen[i]=f.chi2_nopen
        wl_start[i]=f.wl_start
        wl_result[i]=f.wl_result
        sig_start[i]=f.sig_start
        sig_result[i]=f.sig_result
        mjd[i]=f.mjd
        bcv[i]=f.bcv
        delta_bcv[i]=f.delta_bcv
        delta_rv[i]=f.delta_rv
        final_rv[i]=f.final_rv
        status_flag[i]=f.status_flag
        penalty_flag[i]=f.penalty_flag
        chi2_flag[i]=f.chi2_flag
        
    endfor
    
endif else begin

restore, filename=file_tag+'.sav'

n_exp=max(visit)
;get n_wl
trials_per_exp=where(visit eq 1, n_trials_per_exp)
n_wl=n_elements(uniq(wl_start[trials_per_exp]))
n_sig=n_trials_per_exp/n_wl

chi2_arr=dblarr(n_wl, n_sig, n_exp)
chi2_nopen_arr=dblarr(n_wl, n_sig, n_exp)
wl_start_arr=dblarr(n_wl, n_sig, n_exp)
wl_result_arr=dblarr(n_wl, n_sig, n_exp)
sig_start_arr=dblarr(n_wl, n_sig, n_exp)
sig_result_arr=dblarr(n_wl, n_sig, n_exp)
final_rv_arr=dblarr(n_wl, n_sig, n_exp)
status_flag_arr=lonarr(n_wl, n_sig, n_exp)
penalty_flag_arr=lonarr(n_wl, n_sig, n_exp)
chi2_flag_arr=lonarr(n_wl, n_sig, n_exp)

visit_arr=lonarr(n_wl, n_sig, n_exp)
mjd_arr=dblarr(n_wl, n_sig, n_exp)


for exp=0, n_exp-1 do begin
    for iwl=0, n_wl-1 do begin
        i1=(iwl*n_sig)+(exp*n_trials_per_exp)
        i2=i1+n_sig-1
        chi2_arr[iwl,0,exp]=reform(chi2[i1:i2], 1, n_sig)
        chi2_nopen_arr[iwl,0,exp]=reform(chi2_nopen[i1:i2], 1, n_sig)
        wl_start_arr[iwl,0,exp]=reform(wl_start[i1:i2], 1, n_sig)
        wl_result_arr[iwl,0,exp]=reform(wl_result[i1:i2], 1, n_sig)
        sig_start_arr[iwl,0,exp]=reform(sig_start[i1:i2], 1, n_sig)
        sig_result_arr[iwl,0,exp]=reform(sig_result[i1:i2], 1, n_sig)
        final_rv_arr[iwl,0,exp]=reform(final_rv[i1:i2], 1, n_sig)
        status_flag_arr[iwl,0,exp]=reform(status_flag[i1:i2], 1, n_sig)
        penalty_flag_arr[iwl,0,exp]=reform(penalty_flag[i1:i2], 1, n_sig)
        chi2_flag_arr[iwl,0,exp]=reform(chi2_flag[i1:i2], 1, n_sig)
        visit_arr[iwl,0,exp]=reform(visit[i1:i2], 1, n_sig)
        mjd_arr[iwl,0,exp]=reform(mjd[i1:i2], 1, n_sig)
    endfor
endfor

;print, wl_start_arr[*,0,0]
;print, sig_start_arr[*,0,0]
;stop

endelse    
best_rv=dblarr(n_exp)
best_wlrv=dblarr(n_exp)
mjd=mjd_arr[0,0,*]
for exp=0,n_exp-1 do begin
    c2_2d=chi2_arr[*,*,exp]
    rv_2d=final_rv_arr[*,*,exp]
    wl_2d=wl_result_arr[*,*,exp]
    chi2min=min(c2_2d, mindex)
    rvmin=min(abs(rv_2d), rvmindex)
    print, "Exposure: ", exp
    print, "Best wl: ", wl_2d[mindex]    
    print, "Corresponding chi2: ", c2_2d[mindex]
    print, "Corresponding RV: ", rv_2d[mindex]
    print, "---"
    print, "Best RV: ", rv_2d[rvmindex]
    print, "Corresponding chi2: ", c2_2d[rvmindex]
    print, "Corresponding wl: ", wl_2d[rvmindex]
    print, "---"
    if wl_2d[rvmindex] eq wl_2d[mindex] then agreed='YES' else agreed='NO'
    print, "Agreed: " + agreed
    print, "-------------------------------------------"
    
    best_rv[exp]=rv_2d[rvmindex]
    best_wlrv[exp]=rv_2d[mindex]
endfor

plot, mjd, best_wlrv, ps=6
oplot, mjd, best_rv, ps=6, symsize=1.5, color=200
stop
end
