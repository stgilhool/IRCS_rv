pro rvtest_display2


;READIN FILES
file_tag='Mar26_wl_allexp3'
file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_310_398_'+file_tag+'_r0.fits'

fits_info, file, n_ext=n_exp, /silent

;INITIALIZE VECTORS

    
sarr=ptrarr(n_exp, /allocate_heap)

for e=0,n_exp-1 do begin
    s=mrdfits(file, e+1)
    ;if e eq 0 then sarr=ptrarr(n_exp) else sarr[e]=s
    *sarr[e]=s
endfor

ssize=size(*sarr[0])
n_wl=ssize[1]
if ssize[0] eq 1 then n_sig=1L else n_wl=ssize[2]

trials_per_exp=n_sig*n_wl

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

model_arr=dblarr(n_elements((*sarr[0])[0,0].model_select), n_wl, n_sig, n_exp)
obs_arr=dblarr(n_elements((*sarr[0])[0,0].model_select), n_exp)
err_arr=dblarr(n_elements((*sarr[0])[0,0].model_select), n_exp)



for exp=0, n_exp-1 do begin
    chi2_arr[0,0,exp]=(*sarr[exp])[*,*].chi2
    chi2_nopen_arr[0,0,exp]=(*sarr[exp])[*,*].chi2_nopen
    wl_start_arr[0,0,exp]=(*sarr[exp])[*,*].wl_start
    wl_result_arr[0,0,exp]=(*sarr[exp])[*,*].wl_result
    sig_start_arr[0,0,exp]=(*sarr[exp])[*,*].sig_start
    sig_result_arr[0,0,exp]=(*sarr[exp])[*,*].sig_result
    final_rv_arr[0,0,exp]=(*sarr[exp])[*,*].final_rv
    status_flag_arr[0,0,exp]=(*sarr[exp])[*,*].status_flag
    penalty_flag_arr[0,0,exp]=(*sarr[exp])[*,*].penalty_flag
    chi2_flag_arr[0,0,exp]=(*sarr[exp])[*,*].chi2_flag
    visit_arr[0,0,exp]=(*sarr[exp])[*,*].visit
    mjd_arr[0,0,exp]=(*sarr[exp])[*,*].mjd
    
    model_arr[0,0,0,exp]=(*sarr[exp])[*,*].model_select
    obs_arr[0,exp]=(*sarr[exp])[0,0].obs_select
    err_arr[0,exp]=(*sarr[exp])[0,0].err_select

endfor


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
