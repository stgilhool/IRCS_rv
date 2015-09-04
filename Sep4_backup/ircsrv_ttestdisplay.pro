pro ircsrv_ttestdisplay, test=test


;Set up big dots
phi=findgen(32)*(!PI*2/32.)
phi = [ phi, phi(0) ]
usersym, cos(phi), sin(phi), /fill



if n_elements(test) eq 0 then test=1

filename='/home/stgilhool/RV_projects/IRCS_rv/data/telluric_test/test'+strtrim(test,2)+'/tellurictest.fits'

fits_info, filename, n_ext=n_ext, /silent

;n_ext=8L

if n_ext eq 0 then message, "fits file has no extensions"

;Define arrays
t=mrdfits(filename, 1)
template=dblarr(n_elements(t.template), n_ext)
if test gt 8 then begin
    chi2_nopen=dblarr(n_ext)
    ncalls=dblarr(n_ext)
    model_size=size(t.model_arr)
    model_dims=model_size[lindgen(model_size[0])+1]
    model_arr=dblarr([model_dims, n_ext])
    obs_arr=t.obs_arr ;these are the same for all extensions
    err_arr=t.err_arr ;ditto
    n_exp=n_elements(err_arr[0,*])
endif
    
guess=t.guess
chi2=dblarr(n_ext)
pscale=dblarr(n_ext)
status=lonarr(n_ext)


;stop
for ext=1,n_ext do begin
    
    r=mrdfits(filename, ext)
    template[*,ext-1]=r.template
    chi2[ext-1]=r.chi2
    pscale[ext-1]=r.pscale
    status[ext-1]=r.status
    if test gt 8 then begin
        chi2_nopen[ext-1]=r.chi2_nopen
        ncalls[ext-1]=r.ncalls
        model_arr[*,*, ext-1]=r.model_arr
    endif
endfor

if test gt 8 then !p.multi=[0,1,n_ext+2] $
  else !p.multi=[0,1,n_ext+1]
       
 window, 0, ysize=1000

plot, pscale, chi2, ps=6, /xlog
oplot, pscale, chi2
if test gt 8 then begin
    plot, pscale, chi2_nopen, ps=6, /xlog
    oplot, pscale, chi2_nopen
endif
    
;plot, pscale, status, /xlog, ps=8, symsize=0.3, yr=[min(status)-1, max(status)+1]

;!p.multi=[0,1, n_ext]

;window, 1, ysize=750

for ext=1, n_ext do begin
    if test gt 8 then begin
        title_string='TrialNum: ' + strtrim(ext, 2) + $
          ' | pscale: ' + strtrim(pscale[ext-1],2) + $
          ' | chi2: ' + strtrim(chi2[ext-1], 2) + $
          ' | chi2_nopen: ' + strtrim(chi2_nopen[ext-1], 2) + $
          ' | ncalls: ' + strtrim(ncalls[ext-1], 2) + $
          ' | Status: ' + strtrim(status[ext-1],2)
    endif else begin
        title_string='TrialNum: ' + strtrim(ext, 2) + $
          ' | pscale: ' + strtrim(pscale[ext-1],2) + $
          ' | chi2: ' + strtrim(chi2[ext-1], 2) + $
          ' | Status: ' + strtrim(status[ext-1],2)
    endelse
    plot, template[*,ext-1], /xs, yr=[0.4,1.1], title=title_string, charsize=1.5
    oplot, guess, ps=8, symsize=0.3, color=200
endfor

;xnum=0
;while xnum ne 'q' do begin
    
;    read, xnum, prompt='Type exposure number'
    
;    file='/home/stgilhool/RV_projects/IRCS_rv/data/smooth_penalty_test/test3/penalty_'+strtrim(xnum,2)+'.sav'
if test le 8 then begin
    !p.multi=[0,1,2]

    window, 1

;    for exp=0,n-exp-1 do begin
    xnum=0
    plot, guess, /xs, yr=[0,1.1], title="Best one? ("+strtrim(xnum+1,2)+")"
    oplot, template[*,xnum], color=200 ;, ps=8, symsize=0.3
    plot, guess-template[*,xnum], ps=3, yr=[-0.1,0.1];yr=[-0.01,0.01]
    
    window, 2
    
    xnum=9
    plot, guess, /xs, yr=[0,1.1], title="Best one? ("+strtrim(xnum+1,2)+")"
    oplot, template[*,xnum], color=200
    plot, guess-template[*,xnum], ps=3, yr=[-0.1,0.1];yr=[-0.01,0.01]

endif else begin
    
    ;stop
    
    pnum=0
    while pnum ne 99 do begin
        
        read, pnum, prompt='Type pscale trial number from '+strtrim(1,2)+ $
          ' - ' + strtrim(n_ext,2)+' (99 to exit): '
        if pnum gt n_ext and pnum ne 99 then print, "Entered trial is too high" $
        else if pnum le 0 then print, "Entered trial is too low" $
        else if pnum eq 99 then print, "Exiting" $
        else begin
            !p.multi=0
            window, 1
            title_string='TrialNum: ' + strtrim(pnum, 2) + $
              ' | pscale: ' + strtrim(pscale[pnum-1],2) + $
              ' | chi2: ' + strtrim(chi2[pnum-1], 2) + $
              ' | chi2_nopen: ' + strtrim(chi2_nopen[pnum-1], 2) + $
              ' | ncalls: ' + strtrim(ncalls[pnum-1], 2) + $
              ' | Status: ' + strtrim(status[pnum-1],2)
            plot, template[*,pnum-1], /xs, yr=[0,1.1], title=title_string
            oplot, guess, color=200, ps=8, symsize=0.3
            
            expnum=0
            while expnum ne 99 do begin
                
                read, expnum, prompt='Type exposure number from '+strtrim(0,2)+ $
                  ' - ' +strtrim(n_exp-1,2)+' (-1 to display all exposures, 99 to return to trial selection): '
                if expnum ge n_exp and expnum ne 99 then print,  "Entered trial is too high" $
                else if expnum lt -1 then print, "Entered trial is too low" $
                else if expnum eq 99 then print, "Returning to trial selection" $
                else if expnum ge 0 then begin
                    
                    
                    chi2exp=total(((model_arr[*,expnum,pnum-1]-obs_arr[*,expnum]) / $
                                   err_arr[*,expnum])^2, /double)
                    penexp=chi2[pnum-1]-chi2_nopen[pnum-1]
                    title_string='Model and Obs for exposure: '+strtrim(expnum,2)+$
                      ' | Chi2: ' + strtrim(chi2exp,2) + $
                      ' | Smoothness Penalty: ' + strtrim(penexp,2)
                    !p.multi=[0,1,2]
                    
                    window, 2
                    plot, obs_arr[*,expnum], /xs, yr=[0,1.1], title=title_string
                    oplot, model_arr[*,expnum, pnum-1], color=200, ps=8, symsize=0.3
                    
                    plot, obs_arr[*,expnum]-model_arr[*,expnum,pnum-1], yr=[-0.1,0.1], ps=3, /xs
                endif else begin
                    !p.multi=[0,2,n_exp]
                    window, 2, xsize=1300, ysize=1000
                    for frame=0, n_exp-1 do begin
                        chi2exp=total(((model_arr[*,frame,pnum-1]-obs_arr[*,frame]) / $
                                       err_arr[*,frame])^2, /double)
                        penexp=chi2[pnum-1]-chi2_nopen[pnum-1]
                        title_string='Model and Obs for exposure: '+strtrim(frame,2)+$
                          ' | Chi2: ' + strtrim(chi2exp,2) + $
                          ' | Smoothness Penalty: ' + strtrim(penexp,2)
                        
                        plot, obs_arr[*,frame], /xs, yr=[0,1.1], title=title_string
                        oplot, model_arr[*,frame, pnum-1], color=200, ps=8, symsize=0.3
                        
                        plot, obs_arr[*,frame]-model_arr[*,frame,pnum-1], yr=[-0.1,0.1], ps=3, /xs
                    endfor
                    
                endelse
                
                
                
            endwhile
            
        endelse
    endwhile
endelse



!p.multi=0
window, 3

d_chi2=chi2[1:-1]-chi2[0:-2]
d_pscale=pscale[1:-1]-pscale[0:-2]
d_chi2_nopen=chi2_nopen[1:-1]-chi2_nopen[0:-2]
plot, d_chi2/d_pscale

d2_chi2_nopen=d_chi2_nopen[1:-1]-d_chi2_nopen[0:-1]

stop

plot, d_chi2_nopen/d_pscale

stop
                    
                
                         

end
        


