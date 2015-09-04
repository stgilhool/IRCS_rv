function readin, file_tag, first_pix=fp, npix_select=ns, dims=dims, expnum=in, field=field

;dims is output dimensions

if n_elements(ns) eq 0 then message, "pass npix_select, please"
                                ;;Set file tag and name
if n_elements(file_tag) eq 0 then file_tag1='310_398_Mar26_wl_allexp3' else $
  file_tag1=strtrim(fp,2)+'_'+strtrim(ns+fp-1,2)+'_'+file_tag

file='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_'+file_tag1+'_r0.fits'

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
trials_per_exp=n_elements(*(sarr)[fid_exp])


if n_elements(dims) eq 0 then dims=[trials_per_exp, 0] $
  else if n_elements(dims) ne 2 then message, "dims must have 2 dimensions" $
  else if dims[0]*dims[1] ne trials_per_exp then message, "input dimensions don't fit with input array"
    
nx=dims[0]
ny=dims[1]

if n_elements(field) eq 0 then message, "Please specify a field from the structure"

                                ;;check for existence of field
field_check = tag_exist((*sarr[fid_exp]), field, index=field_index)

if field_check eq 0 then message, field+" does not exist"


field_tot_len=n_elements((*sarr[fid_exp]).(field_index))
dims_tot_len=nx*ny
field_len=field_tot_len/dims_tot_len


                                ;;Initialize arrays
if field_len gt 1 then begin
    out_arr=dblarr(field_len, nx, ny, n_exp) 
    
    endif else out_arr=dblarr(nx,ny,n_exp)

                                ;;Populate arrays		
foreach visit, in, exp do begin
    if field_len gt 1 then begin
        out_arr[0,0,0,exp]=reform(((*sarr[exp]).(field_index)), field_len, nx, ny)
    endif else begin
        out_arr[0,0,exp]=reform(((*sarr[exp]).(field_index)), nx, ny)
    endelse

endforeach


return, out_arr

end




pro rvtest_display4, first_pix=fp, npix_select=ns, ft0=ft0, dims=ft0dims, expnum=expnum, copt=copt, bopt=bopt


;if n_elements(ft0) eq 0 then ft0='Mar27_wlrv_allbrute'
if n_elements(fp) eq 0 then fp=310L
if n_elements(ns) eq 0 then ns=89L < (1024L-fp)
if fp+ns gt 1024L then message, "Not enough pixels exist to cover that range"
if n_elements(copt) eq 0 then begin
    sopt=1
    copt=0
endif else sopt=0
if n_elements(bopt) eq 0 then bopt=0
if bopt then copt=1

tag_head=strtrim(fp,2)+'_'+strtrim(fp+ns-1, 2)+'_'

rvwllist=['Mar27_wlrv_allbrute','Mar31_wlrv_allbrute7','Apr01_lsf_brute','Apr03_lsf_brute','Apr03_lsf_brutemore','Apr03_lsf3_brute','Apr04_lsf_sig', 'Apr04_lsf_sigplus']

sigghlist=['Apr07_lsf2_brute','Apr20_lsf2_fine']

lsflist=['Apr21_lsf3_brute','Apr22_lsf3_wide','Apr22_lsf3_more','Apr23_lsf3_many']

if n_elements(ft0) eq 0 then ft0=sigghlist[1]

if n_elements(ft0dims) eq 0 then ft0dims=[41L,41L]


if n_elements(expnum) eq 0 then expnum=lindgen(14)
;expnum=[2,7,9]

nothing=where(lsflist eq ft0, newtest)

if newtest gt 0 then begin
    p0=readin(ft0, first_pix=fp, npix_select=ns, dims=ft0dims, field='par0_start', expnum=expnum)
    p1=readin(ft0, first_pix=fp, npix_select=ns, dims=ft0dims, field='par1_start', expnum=expnum)
endif else begin
    p0=readin(ft0, first_pix=fp, npix_select=ns, dims=ft0dims, field='sig_start', expnum=expnum)
    p1=readin(ft0, first_pix=fp, npix_select=ns, dims=ft0dims, field='gh_start', expnum=expnum)
endelse
c2=readin(ft0, first_pix=fp, npix_select=ns, dims=ft0dims, field='chi2', expnum=expnum)
c2no=readin(ft0, first_pix=fp, npix_select=ns, dims=ft0dims, field='chi2_nopen', expnum=expnum)
wls=readin(ft0, first_pix=fp, npix_select=ns, dims=ft0dims, field='wl_start', expnum=expnum)
wlr=readin(ft0, first_pix=fp, npix_select=ns, dims=ft0dims, field='wl_result', expnum=expnum)
rv=readin(ft0, first_pix=fp, npix_select=ns, dims=ft0dims, field='final_rv', expnum=expnum)
mjd=readin(ft0, first_pix=fp, npix_select=ns, dims=ft0dims, field='mjd', expnum=expnum)
sflag=readin(ft0, first_pix=fp, npix_select=ns, dims=ft0dims, field='status_flag', expnum=expnum)
pflag=readin(ft0, first_pix=fp, npix_select=ns, dims=ft0dims, field='penalty_flag', expnum=expnum)
cflag=readin(ft0, first_pix=fp, npix_select=ns, dims=ft0dims, field='chi2_flag', expnum=expnum)
vis=readin(ft0, first_pix=fp, npix_select=ns, dims=ft0dims, field='visit', expnum=expnum)
stat=readin(ft0, first_pix=fp, npix_select=ns, dims=ft0dims, field='status', expnum=expnum)



c2dims=size(c2, /dimensions)
;n_exp=c2dims[-1]
n_exp=n_elements(expnum)



if sopt eq 1 || bopt eq 1 then begin
    
    for exp=0,n_exp-1 do begin
        
        maxval=max(c2[*,*,exp])
        minval=min(c2[*,*,exp])
        diffval=maxval-(1.1d0*minval)
        num_levels=9L
        increm=diffval/(num_levels-1)
        help, maxval
        help, minval
        help, increm
        
        ;maxarr=[10000L, 5000L, 3000L, 2900L, 2800L, 2700L, 2600L, 2500L, 2400L]
        for w=0,num_levels-1 do begin
            maxvalue=(maxval+1)-(w*increm)
            help, maxvalue
            ;maxvalue=maxarr[w]
            for i=0,36 do begin 
                
                surface, c2[*,*,exp],p0[*,*,exp],p1[*,*,exp], az=i*10L,$
                  max_value=maxvalue, title="chi2 surface exp: "+ $
                  strtrim(exp,2)+" | bestchi2: "+ strtrim(minval,2)+ $
                  " | max_value: "+strtrim(maxvalue,2);, /lego
                
                wait, 0.1 
            endfor 
            wait, 1 
        endfor
        wait,2
    endfor
stop    
endif


if copt eq 1 || bopt eq 1 then begin
    
    exp=0
    minval=min(c2[*,*,exp])
    contour, c2[*,*,exp],p0[*,*,exp],p1[*,*,exp], levels=[minval*1.000001]
    
    for i=1,25 do begin
        contour, c2[*,*,exp],p0[*,*,exp],p1[*,*,exp], levels=[minval*1.000001+lindgen(i+1)*0.05], /overplot
        wait,0.5
        ;if i eq 20 then stop
        
        
    endfor
    
endif

stop

cs=sort(c2)
p0s=p0[cs]
p1s=p1[cs]
c2s=c2[cs]
ws=wlr[cs]
rs=rv[cs]

!p.multi=[0,1,3]
maxpoints=41L


for point=0, maxpoints do begin
    plot, p0s[0:point], p1s[0:point], ps=6, xr=[0.9*min(p0s[0:maxpoints]),1.1*max(p0s[0:maxpoints])], yr=[0.9*min(p1s[0:maxpoints]),1.1*max(p1s[0:maxpoints])]
    plot, lindgen(point+1),rs[0:point], xr=[0,maxpoints], yr=[0.9*min(rs[0:maxpoints]),1.1*max(rs[0:maxpoints])], ps=6
    plo!p.multi=0
plot, rs[0:50], ps=6
t, lindgen(point+1),ws[0:point], xr=[0,maxpoints], yr=[0.9*min(ws[0:maxpoints]),1.1*max(ws[0:maxpoints])],/ys, ps=6
    wait, 0.5
endfor

stop


end

; for tag=0, n_tags-1 do begin
    
;     name=tag_vname[tag]
;     tag_base=scope_varfetch(name)
;     file_tag=tag_head+tag_base
    
;     brute_check=where(brutelist eq tag_base, brute_switch)
;     if brute_switch ge 1 then brute_switch=1 else brute_switch=0
;     siggh_check=where(sigghlist eq tag_base, siggh_switch)
;     if siggh_switch ge 1 then siggh_switch=1 elst siggh_switch=0

;     s=readin(file_tag, npix_select=ns) ;;Readin structure to s

;     r[name]=s                   ;;Store struct in hash
    
    
;     if brute_switch eq 1 then begin      ;;All brute run
        
;         n_exp=n_elements(s.chi2_arr[0,0,*])
        
;         for exp=0,n_exp-1 do begin
        

    
;             if exp eq 7 and tag_base eq 'Mar27_wlrv_allbrute' then begin ;;patching in new run of outlier
;                 s=readin(tag_head+'Mar31_wlrv_allbrute7', npix_select=ns)
                
;             endif

            
                        
;             n_wl=n_elements(uniq(s.wl_start_arr[sort(s.wl_start_arr[*,0,0]),0,0]))
;             n_rv=n_elements(uniq(s.final_rv_arr[sort(s.final_rv_arr[*,0,0]),0,0]))
;             ch=dblarr(n_wl, n_rv)
;             sig=dblarr(n_wl, n_rv)
;             ;help, n_wl
;             ;help, n_rv

;             if exp eq 7 and tag_base eq 'Mar27_wlrv_allbrute' then begin
;                 wl=rebin(s.wl_result_arr[0:n_wl-1,0,0], n_wl, n_rv)
;                 rv=rebin(reform(s.final_rv_arr[lindgen(n_rv)*n_wl, 0,0], 1, n_rv), n_wl, n_rv)
;                 ch[*,*]=s.chi2_arr[*,0,0] 
;                 sig[*,*]=s.sig_result_arr[*,0,0]
;                 s=r[name]
                
;             endif else begin
                
                
;                 wl=rebin(s.wl_result_arr[0:n_wl-1,0,exp], n_wl, n_rv)
;                 rv=rebin(reform(s.final_rv_arr[lindgen(n_rv)*n_wl, 0, exp], 1, n_rv), n_wl, n_rv)
;                 ch[*,*]=s.chi2_arr[*,0,exp]
;                 sig[*,*]=s.sig_result_arr[*,0,exp]
                            
;             endelse
       

;             dkeys=['wl'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'wls'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'sig'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'sigs'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'sf'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'pf'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'cf'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'ex'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'mjd'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'mod'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'obs'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'err'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'rv'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'ch'+strtrim(tag,2)+'_'+strtrim(exp,2) $
;                    ]

;             d[dkeys[0]]=wl
;             d[dkeys[2]]=sig
;             d[dkeys[8]]=s.mjd_arr[*,*,exp]
;             d[dkeys[12]]=rv
;             d[dkeys[13]]=ch
            
;         endfor
;     endif else begin            ;;All other runs
        
;         n_exp=n_elements(s.chi2_arr[0,0,*])
;         n_wl=n_elements(uniq(s.wl_start_arr[sort(s.wl_start_arr[*,*,0]), *,0]))
;         n_sig=n_elements(uniq(s.sig_start_arr[sort(s.sig_start_arr[*,*,0]), *,0]))
;         ;print, n_wl
;         ;print, n_sig
        
;         n_rv=n_wl

;         for exp=0,n_exp-1 do begin
            
;             wl=s.wl_result_arr[*,*,exp]
;             wls=s.wl_start_arr[*,*,exp]
;             sig=s.sig_result_arr[*,*,exp]
;             sigs=s.sig_start_arr[*,*,exp]
;             sf=s.status_flag_arr[*,*,exp]
;             pf=s.penalty_flag_arr[*,*,exp]
;             cf=s.chi2_flag_arr[*,*,exp]
;             ex=s.visit_arr[*,*,exp]
;             mjd=s.mjd_arr[*,*,exp]
;             rv=s.final_rv_arr[*,*,exp]
;             ch=s.chi2_arr[*,*,exp]
            

           
;             dkeys=['wl'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'wls'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'sig'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'sigs'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'sf'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'pf'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'cf'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'ex'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'mjd'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'mod'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'obs'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'err'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'rv'+strtrim(tag,2)+'_'+strtrim(exp,2),$
;                    'ch'+strtrim(tag,2)+'_'+strtrim(exp,2) $
;                    ]         
            
;             d[dkeys[0]]=wl
;             d[dkeys[1]]=wls
;             d[dkeys[2]]=sig
;             d[dkeys[3]]=sigs
;             d[dkeys[4]]=sf
;             d[dkeys[5]]=pf
;             d[dkeys[6]]=cf
;             d[dkeys[7]]=ex
;             d[dkeys[8]]=mjd
;             d[dkeys[12]]=rv
;             d[dkeys[13]]=ch
            
            
;         endfor
;     endelse

  
    
    
; endfor

; o=d.toStruct()


; rvbest0=dblarr(n_exp)
; rvbest1=dblarr(n_exp)
; chbest0=dblarr(n_exp)
; chbest1=dblarr(n_exp)
; mjd0=dblarr(n_exp)
; mjd1=dblarr(n_exp)
; rvbestall=dblarr(n_exp)
; wlbestall=dblarr(n_exp)
; sigbestall=dblarr(n_exp)
; chdiff=dblarr(n_exp)
; chbestall=dblarr(n_exp)

; for exp=0,n_exp-1 do begin
;     ex=strtrim(exp,2)
    
;     ;Best from contours
;     chbv=d['ch0_'+ex]
;     chb0rough=min(chbv, mind)
;     chbv_nx=n_elements(chbv[*,0])
;     chbv_ny=n_elements(chbv[0,*])

;     mind2d=array_indices(d['ch0_'+ex], mind)
;     pix_range=7L

;     wlbv=d['wl0_'+ex]
;     wlb0rough=wlbv[mind]
;     wlbv_nx=n_elements(wlbv[*,0])
;     wlbv_ny=n_elements(wlbv[0,*])
        

;     frangex=[mind2d[0]-pix_range, mind2d[0]+pix_range]
;     frangey=[mind2d[1]-pix_range, mind2d[1]+pix_range]
    
    
    
;     rvbv= d['rv0_'+ex]
;     rvb0rough=rvbv[mind]
    
;     rvbv_nx=n_elements(rvbv[*,0])
;     rvbv_ny=n_elements(rvbv[0,*])


;     sigbv= d['sig0_'+ex]
;     sigb0rough=sigbv[mind]
    
;     sigbv_nx=n_elements(sigbv[*,0])
;     sigbv_ny=n_elements(sigbv[0,*])
    
;     wlb_arr=wlbv[(0L > frangex[0]):((wlbv_nx-1) < frangex[1]), (0L > frangey[0]):((wlbv_ny-1) < frangey[1])]
;     rvb_arr=rvbv[(0L > frangex[0]):((rvbv_nx-1) < frangex[1]), (0L > frangey[0]):((rvbv_ny-1) < frangey[1])]
;     chb_arr=chbv[(0L > frangex[0]):((chbv_nx-1) < frangex[1]), (0L > frangey[0]):((chbv_ny-1) < frangey[1])]
;     sigb_arr=sigbv[(0L > frangex[0]):((sigbv_nx-1) < frangex[1]), (0L > frangey[0]):((sigbv_ny-1) < frangey[1])]

    
;     ;surf_interp=tri_surf(chb_arr, xvalues=wlb_arr, yvalues=rvb_arr, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)
;     surf_interp=tri_surf(chb_arr, wlb_arr, rvb_arr, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)
;     surf_x=tri_surf(wlb_arr, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)
;     surf_y=tri_surf(rvb_arr, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)

;     surf_sig=tri_surf(sigb_arr, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)
    
;     chb0=min(surf_interp, fmind)
;     wlb0=surf_x[fmind]
;     rvb0=surf_y[fmind]
;     sigb0=surf_sig[fmind]

;     ;print, "old ... interpolated"
;     ;print, chb0rough, chb0
;     ;print, wlb0rough, wlb0
;     ;print, rvb0rough, rvb0
    
;     ;stop
    
;     if brute_switch eq 0 then begin ;;FIX THIS IS SEMI-HARD CODED
;         ;Best from mpfit
;         chb1=min(d['ch1_'+ex], mind1)
;         wlbv1=d['wl1_'+ex]
;         wlb1=wlbv1[mind1]
;         rvbv1= d['rv1_'+ex]
;         rvb1=rvbv1[mind1]
;         sigbv1=d['sig1_'+ex]
;         sigb1=sigbv1[mind1]
        
;         ;contour, d['ch0_'+ex], d['wl0_'+ex], d['rv0_'+ex],
;         ;levels=[chb0*1.00001+lindgen(30)], title="Contour of exp: "+ex+"
;         ;| Best_fit: chi2="+strtrim(chb0,2)+", RV="+strtrim(rvb0,2)+",
;         ;WL="+strtrim(wlb0,2), charsize=1.5, ytitle="RV (km/s)",
;         ;xtitle="WL 0pt shift (um)"
        
;     endif else begin
;         chbv1=d['ch1_'+ex]
;         chb1rough=min(chbv1, mind1)
;         chbv1_nx=n_elements(chbv1[*,0])
;         chbv1_ny=n_elements(chbv1[0,*])
        
;         mind2d1=array_indices(d['ch1_'+ex], mind1)
;         pix_range=2L
        
;         wlbv1=d['wl1_'+ex]
;         wlb1rough=wlbv1[mind1]
;         wlbv1_nx=n_elements(wlbv1[*,0])
;         wlbv1_ny=n_elements(wlbv1[0,*])
        
;         frangex1=[mind2d1[0]-pix_range, mind2d1[0]+pix_range]
;         frangey1=[mind2d1[1]-pix_range, mind2d1[1]+pix_range]
    
        
;         wlb_arr1=wlbv1[(0L > frangex1[0]):((wlbv1_nx-1) < frangex1[1]), (0L > frangey1[0]):((wlbv1_ny-1) < frangey1[1])]
    
;         rvbv1= d['rv1_'+ex]
;         rvb1rough=rvbv1[mind1]
;         rvbv1_nx=n_elements(rvbv1[*,0])
;         rvbv1_ny=n_elements(rvbv1[0,*])
;         rvb_arr1=rvbv1[(0L > frangex1[0]):((rvbv1_nx-1) < frangex1[1]), (0L > frangey1[0]):((rvbv1_ny-1) < frangey1[1])]
        
;         chb_arr1=chbv1[(0L > frangex1[0]):((chbv1_nx-1) < frangex1[1]), (0L > frangey1[0]):((chbv1_ny-1) < frangey1[1])]
    
    
;     ;surf_interp=tri_surf(chb_arr, xvalues=wlb_arr, yvalues=rvb_arr, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)
;     surf_interp1=min_curve_surf(chb_arr1, wlb_arr1, rvb_arr1,/tps, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)
;     surf_x1=tri_surf(wlb_arr1, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)
;     surf_y1=tri_surf(rvb_arr1, nx=(2*pix_range+1)*5, ny=(2*pix_range+1)*5)
    
    
;     chb1=min(surf_interp1, fmind1)
;     wlb1=surf_x1[fmind1]
;     rvb1=surf_y1[fmind1]
;     endelse


; ;;;PLOT
; contour, surf_interp, surf_x, surf_y, levels=[chb0*1.00001+lindgen(30)], title="Contour of exp: "+ex+" | Best_fit: chi2="+strtrim(chb0,2)+", RV="+strtrim(rvb0,2)+", WL="+strtrim(wlb0,2), charsize=1.5, ytitle="RV (km/s)", xtitle="WL 0pt shift (um)"
    

; ;contour, chbv1, wlbv1, rvbv1, levels=[chb1*1.00001+lindgen(30)], title="Contour of exp: "+ex+" | Best_fit: chi2="+strtrim(chb1,2)+", RV="+strtrim(rvb1,2)+", WL="+strtrim(wlb1,2), charsize=1.5, ytitle="RV (km/s)", xtitle="WL 0pt shift (um)"
    
; ;if brute_switch eq 0 then begin ;;FIX same as above    
;     plot, [wlb1],[rvb1], ps=6, color=200, xs=4, ys=4, /noerase, xr=[!x.crange], yr=[!y.crange]
;     xyouts,!x.crange[0]+1d-7, !y.crange[1]-0.05, 'Chi2: '+strtrim(chb1,2)+' | RV: '+strtrim(rvb1,2)+' | WL: '+strtrim(wlb1,2), color=200, charsize=2.0
; ;endif else begin
    


; rvbest0[exp]=rvb0
; rvbest1[exp]=rvb1
; chbest0[exp]=chb0
; chbest1[exp]=chb1


; mjd0[exp]=(d['mjd0_'+ex])[0,0]
; mjd1[exp]=(d['mjd0_'+ex])[0,0]

; if chb0 lt chb1 then begin
;     rvbestall[exp]=rvb0
;     chdiff[exp]=chb0-chb1
;     chbestall[exp]=chb0
;     sigbestall[exp]=sigb0
;     wlbestall[exp]=wlb0
; endif else begin
;     rvbestall[exp]=rvb1
;     chdiff[exp]=chb1-chb0
;     chbestall[exp]=chb1
;     sigbestall[exp]=sigb1
;     wlbestall[exp]=wlb1
; endelse

; stop

; endfor

; ;;;Plot mjd vs rv
; plot, mjd0, rvbest0, ps=6, yr=[(min(rvbest0)) < (min(rvbest1)),(max(rvbest0)) > (max(rvbest1))], /xs, title=RVs

; oplot, mjd1, rvbest1, ps=6, color=200

; xyouts, mjd0, rvbestall, strtrim(chdiff,2), color=99999

; stop

; openw,lun, 'startingparams.txt', /get_lun
; for i=0,n_exp-1 do begin
;     printf, lun, wlbestall[i], sigbestall[i], chbestall[i], rvbestall[i]
; endfor
; close, lun
; free_lun, lun
; stop

;          ;;Set up wl vs rv

; best_rv=dblarr(n_exp)
; best_wlrv=dblarr(n_exp)
; mjd=mjd_arr[0,0,*]
; for exp=0,n_exp-1 do begin
;     ch_2d=chi2_arr[*,*,exp]
;     rv_2d=final_rv_arr[*,*,exp]
;     wl_2d=wl_result_arr[*,*,exp]
;     chi2min=min(ch_2d, mindex)
;     rvmin=min(abs(rv_2d), rvmindex)
;     print, "Exposure: ", exp
;     print, "Best wl: ", wl_2d[mindex]    
;     print, "Corresponding RV: ", rv_2d[mindex]    
;     print, "Corresponding chi2: ", ch_2d[mindex]
    
;     print, "---"

;     print, "Best RV: ", rv_2d[rvmindex]
;     print, "Corresponding wl: ", wl_2d[rvmindex]
;     print, "Corresponding chi2: ", ch_2d[rvmindex]

;     print, "---"

;     if wl_2d[rvmindex] eq wl_2d[mindex] then agreed='YES' else agreed='NO'
;     print, "Agreed: " + agreed
;     print, "-------------------------------------------"
    
;     best_rv[exp]=rv_2d[rvmindex]
;     best_wlrv[exp]=rv_2d[mindex]
; endfor

; plot, mjd, best_wlrv, ps=6
; oplot, mjd, best_rv, ps=6, symsize=1.5, color=200
; stop


; end
