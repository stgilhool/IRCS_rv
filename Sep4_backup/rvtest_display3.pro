function readin, file_tag, index=in, npix_select=ns

if n_elements(ns) eq 0 then message, "pass npix_select, please"
                                ;;Set file tag and name
if n_elements(file_tag) eq 0 then file_tag='310_398_Mar26_wl_allexp3'

file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_'+file_tag+'_r0.fits'

                                ;;Get appropriate indices
fits_info, file, n_ext=n_exp, /silent
if n_elements(in) eq 0 then in=lindgen(n_exp) else n_exp=n_elements(in)


;INITIALIZE VECTORS

sarr=ptrarr(n_exp, /allocate_heap)
                                ;;Populate pointer array with structures
foreach visit, in, index do begin
    
    s=mrdfits(file, visit+1)
    *sarr[index]=s
endforeach

fid_exp=0
                                ;;Get dimensions of arrays in str to initialize arrays
ssize=size(*sarr[0])
n_wl=ssize[1]

n_bases=n_elements((*sarr[0])[0,0].sig_start)

if ssize[0] eq 1 then n_sig=1L else n_sig=ssize[2]

;print, ssize

trials_per_exp=n_sig*n_wl

;help, file_tag
;help, n_exp
;help, in
;help, n_wl
;help, sarr
;help, s
;stop



                                ;;Initialize arrays
chi2_arr=dblarr(n_wl, n_sig, n_exp)
chi2_nopen_arr=dblarr(n_wl, n_sig, n_exp)
wl_start_arr=dblarr(n_wl, n_sig, n_exp)
wl_result_arr=dblarr(n_wl, n_sig, n_exp)

if n_bases eq 1 then begin
    sig_start_arr=dblarr(n_wl, n_sig, n_exp)
    sig_result_arr=dblarr(n_wl, n_sig, n_exp)
endif else begin
    sig_start_arr=dblarr(n_bases,n_wl, n_sig, n_exp)
    sig_result_arr=dblarr(n_bases,n_wl, n_sig, n_exp)
endelse
final_rv_arr=dblarr(n_wl, n_sig, n_exp)
status_flag_arr=lonarr(n_wl, n_sig, n_exp)
penalty_flag_arr=lonarr(n_wl, n_sig, n_exp)
chi2_flag_arr=lonarr(n_wl, n_sig, n_exp)
visit_arr=lonarr(n_wl, n_sig, n_exp)
mjd_arr=dblarr(n_wl, n_sig, n_exp)
                                ;;Spectral arrays
spectra_check = tag_exist((*sarr[fid_exp]),'model_select')

model_arr=dblarr(ns, n_wl, n_sig, n_exp)
obs_arr=dblarr(ns, n_exp)
err_arr=dblarr(ns, n_exp)

                                
                                
                                ;;Populate arrays		
foreach visit, in, exp do begin
    chi2_arr[0,0,exp]=(*sarr[exp])[*,*].chi2
    chi2_nopen_arr[0,0,exp]=(*sarr[exp])[*,*].chi2_nopen
    wl_start_arr[0,0,exp]=(*sarr[exp])[*,*].wl_start
    wl_result_arr[0,0,exp]=(*sarr[exp])[*,*].wl_result
    if n_bases eq 1 then begin
        sig_start_arr[0,0,exp]=(*sarr[exp])[*,*].sig_start
        sig_result_arr[0,0,exp]=(*sarr[exp])[*,*].sig_result
    endif else begin
        sig_start_arr[0,0,0,exp]=(*sarr[exp])[*,*].sig_start
        sig_result_arr[0,0,0,exp]=(*sarr[exp])[*,*].sig_result
    endelse
    final_rv_arr[0,0,exp]=(*sarr[exp])[*,*].final_rv
    status_flag_arr[0,0,exp]=(*sarr[exp])[*,*].status_flag
    penalty_flag_arr[0,0,exp]=(*sarr[exp])[*,*].penalty_flag
    chi2_flag_arr[0,0,exp]=(*sarr[exp])[*,*].chi2_flag
    visit_arr[0,0,exp]=(*sarr[exp])[*,*].visit
    mjd_arr[0,0,exp]=(*sarr[exp])[*,*].mjd
                                ;;Spectral Arrays
    if spectra_check eq 1 then begin
        model_arr[0,0,0,exp]=(*sarr[exp])[*,*].model_select
        obs_arr[0,exp]=(*sarr[exp])[0,0].obs_select
        err_arr[0,exp]=(*sarr[exp])[0,0].err_select
    endif else begin
        model_arr[*,*,*,exp]=make_array(type=5, value=-1,size=size(model_arr[*,*,*,exp]))
        obs_arr[*,exp]=make_array(type=5, value=-1,size=size(obs_arr[*,exp]))
        err_arr[*,exp]=make_array(type=5, value=-1,size=size(err_arr[*,exp]))
    endelse

endforeach


outstr={chi2_arr:chi2_arr, $
        chi2_nopen_arr:chi2_nopen_arr, $
        wl_start_arr:wl_start_arr, $
        wl_result_arr:wl_result_arr, $
        sig_start_arr:sig_start_arr, $
        sig_result_arr:sig_result_arr, $
        final_rv_arr:final_rv_arr, $
        status_flag_arr:status_flag_arr, $
        penalty_flag_arr:penalty_flag_arr, $
        chi2_flag_arr:chi2_flag_arr, $
        visit_arr:visit_arr, $
        mjd_arr:mjd_arr, $
        model_arr:model_arr, $
        obs_arr:obs_arr, $
        err_arr:err_arr  $
       }


return, outstr

end




pro rvtest_display3, ft1, ft2, ft3, ft4, ft5, ft6, ft7, ft8, ft9, ft10, first_pix=fp, npix_select=ns, ft0=ft0


if n_elements(ft0) eq 0 then ft0='Mar27_wlrv_allbrute'
if n_elements(fp) eq 0 then fp=310L
if n_elements(ns) eq 0 then ns=89L < (1024L-fp)
if fp+ns gt 1024L then message, "Not enough pixels exist to cover that range"

tag_head=strtrim(fp,2)+'_'+strtrim(fp+ns-1, 2)+'_'

tag_vname=['ft0','ft1','ft2','ft3','ft4','ft5','ft6','ft7','ft8','ft9','ft10'] 

vname=['rv1', 'rv2']

n_tags=n_params()+1

r=hash()                        ;;Holds output structures from runs

d=hash()                        ;;Holds variables and values from runs

brutelist=['Mar27_wlrv_allbrute','Mar27_wlrv_allbrute7','Apr01_lsf_brute','Apr03_lsf_brute','Apr03_lsf_brutemore','Apr03_lsf3_brute','Apr04_lsf_sig', 'Apr04_lsf_sigplus']

for tag=0, n_tags-1 do begin
    
    name=tag_vname[tag]
    tag_base=scope_varfetch(name)
    file_tag=tag_head+tag_base
    
    brute_check=where(brutelist eq tag_base, brute_switch)
    if brute_switch ge 1 then brute_switch=1 else brute_switch=0
    

    s=readin(file_tag, npix_select=ns) ;;Readin structure to s

    r[name]=s                   ;;Store struct in hash
    
    
    if brute_switch eq 1 then begin      ;;All brute run
        
        n_exp=n_elements(s.chi2_arr[0,0,*])
        
        for exp=0,n_exp-1 do begin
        

    
            if exp eq 7 then begin ;;patching in new run of outlier
                s=readin(tag_head+'Mar31_wlrv_allbrute7', npix_select=ns)
                
            endif

            
                        
            n_wl=n_elements(uniq(s.wl_start_arr[sort(s.wl_start_arr[*,0,0]),0,0]))
            n_rv=n_elements(uniq(s.final_rv_arr[sort(s.final_rv_arr[*,0,0]),0,0]))
            ch=dblarr(n_wl, n_rv)
            sig=dblarr(n_wl, n_rv)
            ;help, n_wl
            ;help, n_rv

            if exp eq 7 then begin
                wl=rebin(s.wl_result_arr[0:n_wl-1,0,0], n_wl, n_rv)
                rv=rebin(reform(s.final_rv_arr[lindgen(n_rv)*n_wl, 0,0], 1, n_rv), n_wl, n_rv)
                ch[*,*]=s.chi2_arr[*,0,0] 
                sig[*,*]=s.sig_result_arr[*,0,0]
                s=r[name]
                
            endif else begin
                
                
                wl=rebin(s.wl_result_arr[0:n_wl-1,0,exp], n_wl, n_rv)
                rv=rebin(reform(s.final_rv_arr[lindgen(n_rv)*n_wl, 0, exp], 1, n_rv), n_wl, n_rv)
                ch[*,*]=s.chi2_arr[*,0,exp]
                sig[*,*]=s.sig_result_arr[*,0,exp]
                            
            endelse
       

            dkeys=['wl'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'wls'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'sig'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'sigs'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'sf'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'pf'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'cf'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'ex'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'mjd'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'mod'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'obs'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'err'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'rv'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'ch'+strtrim(tag,2)+'_'+strtrim(exp,2) $
                   ]

            d[dkeys[0]]=wl
            d[dkeys[2]]=sig
            d[dkeys[8]]=s.mjd_arr[*,*,exp]
            d[dkeys[12]]=rv
            d[dkeys[13]]=ch
            
        endfor
    endif else begin            ;;All other runs
        
        n_exp=n_elements(s.chi2_arr[0,0,*])
        n_wl=n_elements(uniq(s.wl_start_arr[sort(s.wl_start_arr[*,*,0]), *,0]))
        n_sig=n_elements(uniq(s.sig_start_arr[sort(s.sig_start_arr[*,*,0]), *,0]))
        ;print, n_wl
        ;print, n_sig
        
        n_rv=n_wl

        for exp=0,n_exp-1 do begin
            
            wl=s.wl_result_arr[*,*,exp]
            wls=s.wl_start_arr[*,*,exp]
            sig=s.sig_result_arr[*,*,exp]
            sigs=s.sig_start_arr[*,*,exp]
            sf=s.status_flag_arr[*,*,exp]
            pf=s.penalty_flag_arr[*,*,exp]
            cf=s.chi2_flag_arr[*,*,exp]
            ex=s.visit_arr[*,*,exp]
            mjd=s.mjd_arr[*,*,exp]
            rv=s.final_rv_arr[*,*,exp]
            ch=s.chi2_arr[*,*,exp]
            

           
            dkeys=['wl'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'wls'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'sig'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'sigs'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'sf'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'pf'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'cf'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'ex'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'mjd'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'mod'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'obs'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'err'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'rv'+strtrim(tag,2)+'_'+strtrim(exp,2),$
                   'ch'+strtrim(tag,2)+'_'+strtrim(exp,2) $
                   ]         
            
            d[dkeys[0]]=wl
            d[dkeys[1]]=wls
            d[dkeys[2]]=sig
            d[dkeys[3]]=sigs
            d[dkeys[4]]=sf
            d[dkeys[5]]=pf
            d[dkeys[6]]=cf
            d[dkeys[7]]=ex
            d[dkeys[8]]=mjd
            d[dkeys[12]]=rv
            d[dkeys[13]]=ch
            
            
        endfor
    endelse

  
    
    
endfor

o=d.toStruct()


rvbest0=dblarr(n_exp)
rvbest1=dblarr(n_exp)
chbest0=dblarr(n_exp)
chbest1=dblarr(n_exp)
mjd0=dblarr(n_exp)
mjd1=dblarr(n_exp)
rvbestall=dblarr(n_exp)
wlbestall=dblarr(n_exp)
sigbestall=dblarr(n_exp)
chdiff=dblarr(n_exp)
chbestall=dblarr(n_exp)

for exp=0,n_exp-1 do begin
    ex=strtrim(exp,2)
    
    ;Best from contours
    chbv=d['ch0_'+ex]
    chb0rough=min(chbv, mind)
    chbv_nx=n_elements(chbv[*,0])
    chbv_ny=n_elements(chbv[0,*])

    mind2d=array_indices(d['ch0_'+ex], mind)
    pix_range=7L

    wlbv=d['wl0_'+ex]
    wlb0rough=wlbv[mind]
    wlbv_nx=n_elements(wlbv[*,0])
    wlbv_ny=n_elements(wlbv[0,*])
        

    frangex=[mind2d[0]-pix_range, mind2d[0]+pix_range]
    frangey=[mind2d[1]-pix_range, mind2d[1]+pix_range]
    
    
    
    rvbv= d['rv0_'+ex]
    rvb0rough=rvbv[mind]
    
    rvbv_nx=n_elements(rvbv[*,0])
    rvbv_ny=n_elements(rvbv[0,*])


    sigbv= d['sig0_'+ex]
    sigb0rough=sigbv[mind]
    
    sigbv_nx=n_elements(sigbv[*,0])
    sigbv_ny=n_elements(sigbv[0,*])
    
    wlb_arr=wlbv[(0L > frangex[0]):((wlbv_nx-1) < frangex[1]), (0L > frangey[0]):((wlbv_ny-1) < frangey[1])]
    rvb_arr=rvbv[(0L > frangex[0]):((rvbv_nx-1) < frangex[1]), (0L > frangey[0]):((rvbv_ny-1) < frangey[1])]
    chb_arr=chbv[(0L > frangex[0]):((chbv_nx-1) < frangex[1]), (0L > frangey[0]):((chbv_ny-1) < frangey[1])]
    sigb_arr=sigbv[(0L > frangex[0]):((sigbv_nx-1) < frangex[1]), (0L > frangey[0]):((sigbv_ny-1) < frangey[1])]

    
    ;surf_interp=tri_surf(chb_arr, xvalues=wlb_arr, yvalues=rvb_arr, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)
    surf_interp=tri_surf(chb_arr, wlb_arr, rvb_arr, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)
    surf_x=tri_surf(wlb_arr, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)
    surf_y=tri_surf(rvb_arr, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)

    surf_sig=tri_surf(sigb_arr, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)
    
    chb0=min(surf_interp, fmind)
    wlb0=surf_x[fmind]
    rvb0=surf_y[fmind]
    sigb0=surf_sig[fmind]

    ;print, "old ... interpolated"
    ;print, chb0rough, chb0
    ;print, wlb0rough, wlb0
    ;print, rvb0rough, rvb0
    
    ;stop
    
    if brute_switch eq 0 then begin ;;FIX THIS IS SEMI-HARD CODED
        ;Best from mpfit
        chb1=min(d['ch1_'+ex], mind1)
        wlbv1=d['wl1_'+ex]
        wlb1=wlbv1[mind1]
        rvbv1= d['rv1_'+ex]
        rvb1=rvbv1[mind1]
        sigbv1=d['sig1_'+ex]
        sigb1=sigbv1[mind1]
        
        ;contour, d['ch0_'+ex], d['wl0_'+ex], d['rv0_'+ex],
        ;levels=[chb0*1.00001+lindgen(30)], title="Contour of exp: "+ex+"
        ;| Best_fit: chi2="+strtrim(chb0,2)+", RV="+strtrim(rvb0,2)+",
        ;WL="+strtrim(wlb0,2), charsize=1.5, ytitle="RV (km/s)",
        ;xtitle="WL 0pt shift (um)"
        
    endif else begin
        chbv1=d['ch1_'+ex]
        chb1rough=min(chbv1, mind1)
        chbv1_nx=n_elements(chbv1[*,0])
        chbv1_ny=n_elements(chbv1[0,*])
        
        mind2d1=array_indices(d['ch1_'+ex], mind1)
        pix_range=2L
        
        wlbv1=d['wl1_'+ex]
        wlb1rough=wlbv1[mind1]
        wlbv1_nx=n_elements(wlbv1[*,0])
        wlbv1_ny=n_elements(wlbv1[0,*])
        
        frangex1=[mind2d1[0]-pix_range, mind2d1[0]+pix_range]
        frangey1=[mind2d1[1]-pix_range, mind2d1[1]+pix_range]
    
        
        wlb_arr1=wlbv1[(0L > frangex1[0]):((wlbv1_nx-1) < frangex1[1]), (0L > frangey1[0]):((wlbv1_ny-1) < frangey1[1])]
    
        rvbv1= d['rv1_'+ex]
        rvb1rough=rvbv1[mind1]
        rvbv1_nx=n_elements(rvbv1[*,0])
        rvbv1_ny=n_elements(rvbv1[0,*])
        rvb_arr1=rvbv1[(0L > frangex1[0]):((rvbv1_nx-1) < frangex1[1]), (0L > frangey1[0]):((rvbv1_ny-1) < frangey1[1])]
        
        chb_arr1=chbv1[(0L > frangex1[0]):((chbv1_nx-1) < frangex1[1]), (0L > frangey1[0]):((chbv1_ny-1) < frangey1[1])]
    
    
    ;surf_interp=tri_surf(chb_arr, xvalues=wlb_arr, yvalues=rvb_arr, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)
    surf_interp1=min_curve_surf(chb_arr1, wlb_arr1, rvb_arr1,/tps, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)
    surf_x1=tri_surf(wlb_arr1, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)
    surf_y1=tri_surf(rvb_arr1, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)
    
    
    chb1=min(surf_interp1, fmind1)
    wlb1=surf_x1[fmind1]
    rvb1=surf_y1[fmind1]
    endelse


;;;PLOT
contour, surf_interp, surf_x, surf_y, levels=[chb0*1.00001+lindgen(30)], title="Contour of exp: "+ex+" | Best_fit: chi2="+strtrim(chb0,2)+", RV="+strtrim(rvb0,2)+", WL="+strtrim(wlb0,2), charsize=1.5, ytitle="RV (km/s)", xtitle="WL 0pt shift (um)"
    

;contour, chbv1, wlbv1, rvbv1, levels=[chb1*1.00001+lindgen(30)], title="Contour of exp: "+ex+" | Best_fit: chi2="+strtrim(chb1,2)+", RV="+strtrim(rvb1,2)+", WL="+strtrim(wlb1,2), charsize=1.5, ytitle="RV (km/s)", xtitle="WL 0pt shift (um)"
    
;if brute_switch eq 0 then begin ;;FIX same as above    
    plot, [wlb1],[rvb1], ps=6, color=200, xs=4, ys=4, /noerase, xr=[!x.crange], yr=[!y.crange]
    xyouts,!x.crange[0]+1d-7, !y.crange[1]-0.05, 'Chi2: '+strtrim(chb1,2)+' | RV: '+strtrim(rvb1,2)+' | WL: '+strtrim(wlb1,2), color=200, charsize=2.0
;endif else begin
    


rvbest0[exp]=rvb0
rvbest1[exp]=rvb1
chbest0[exp]=chb0
chbest1[exp]=chb1


mjd0[exp]=(d['mjd0_'+ex])[0,0]
mjd1[exp]=(d['mjd0_'+ex])[0,0]

if chb0 lt chb1 then begin
    rvbestall[exp]=rvb0
    chdiff[exp]=chb0-chb1
    chbestall[exp]=chb0
    sigbestall[exp]=sigb0
    wlbestall[exp]=wlb0
endif else begin
    rvbestall[exp]=rvb1
    chdiff[exp]=chb1-chb0
    chbestall[exp]=chb1
    sigbestall[exp]=sigb1
    wlbestall[exp]=wlb1
endelse

stop

endfor

;;;Plot mjd vs rv
plot, mjd0, rvbest0, ps=6, yr=[(min(rvbest0)) < (min(rvbest1)),(max(rvbest0)) > (max(rvbest1))], /xs, title=RVs

oplot, mjd1, rvbest1, ps=6, color=200

xyouts, mjd0, rvbestall, strtrim(chdiff,2), color=99999

stop

openw,lun, 'startingparams.txt', /get_lun
for i=0,n_exp-1 do begin
    printf, lun, wlbestall[i], sigbestall[i], chbestall[i], rvbestall[i]
endfor
close, lun
free_lun, lun
stop

         ;;Set up wl vs rv

best_rv=dblarr(n_exp)
best_wlrv=dblarr(n_exp)
mjd=mjd_arr[0,0,*]
for exp=0,n_exp-1 do begin
    ch_2d=chi2_arr[*,*,exp]
    rv_2d=final_rv_arr[*,*,exp]
    wl_2d=wl_result_arr[*,*,exp]
    chi2min=min(ch_2d, mindex)
    rvmin=min(abs(rv_2d), rvmindex)
    print, "Exposure: ", exp
    print, "Best wl: ", wl_2d[mindex]    
    print, "Corresponding RV: ", rv_2d[mindex]    
    print, "Corresponding chi2: ", ch_2d[mindex]
    
    print, "---"

    print, "Best RV: ", rv_2d[rvmindex]
    print, "Corresponding wl: ", wl_2d[rvmindex]
    print, "Corresponding chi2: ", ch_2d[rvmindex]

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
