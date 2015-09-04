pro IRCS_cal_control_116, nparam_wl, nparam_gh, nparam_gh_lin, nparam_other, niter, random=random, wl_guess=wl_guess, gh_guess=gh_guess, gh_lin_guess=gh_lin_guess, other_guess=other_guess, min_type=min_type, visualize=visualize

;Set lin_switch
if nparam_gh_lin ne 0 then begin
    lin_switch=1 
    if nparam_gh_lin ne nparam_gh then message, "Unequal number of 0th and linear parameters for gh bases"
    endif else lin_switch=0

;Set default minimization scheme
if n_elements(min_type) eq 0 then min_type='amoeba'

;Set default visualization
if n_elements(visualize) eq 0 then visualize=0

;Turn random on/off
if n_elements(random) eq 0 then random=0

;;;Wavelength parameters for amoeba
    ;Potential entries in wl_guess array
wl_guess_all=dblarr(7)
wl_guess_all[0]=2.2912220d0 ;2.2912188d0                     
wl_guess_all[1]= 6.2510646d-05 ;6.2542585d-5                    
wl_guess_all[2]=-2.0597915d-09 ;-2.063680d-9                   
wl_guess_all[3]=-1.5431909d-12 ;-1.6401032d-12                  
wl_guess_all[4]=-8.4728472d-17 ;-8.863d-17
wl_guess_all[5]=0d0
wl_guess_all[6]=0d0
                                ;Same for scale
wl_scale_all=dblarr(7)
wl_scale_all[0]=5d-5
wl_scale_all[1]=5d-8
wl_scale_all[2]=1d-11
wl_scale_all[3]=1d-13
wl_scale_all[4]=5d-18
wl_scale_all[5]=1d-23
wl_scale_all[6]=1d-31
    
;;;Gauss-Hermite parameters for amoeba
gh_guess_all=dblarr(11)


                                ; gh_guess_all[0]=5.9240583d0
                                ; ;6.1341822d0		;guesses
                                ; before cb_hermite fix
;  gh_guess_all[1]=0.025492432d0 
;  gh_guess_all[2]=-0.19622154d0 
;  gh_guess_all[3]=0.020708712d0 
;  gh_guess_all[4]=-0.018781381d0
;  gh_guess_all[5]=0.016140888d0
;  gh_guess_all[6]=0.13626593d0
;  gh_guess_all[7]=-0.11817643d0
;  gh_guess_all[8]=0.079964228d0
;  gh_guess_all[9]=-0.097312890d0
;  gh_guess_all[10]=-0.12176638d0
 
gh_guess_all[0]=5.4801202d0	;guess for gaussian sigma
 gh_guess_all[1]=0.074284679d0
 gh_guess_all[2]=-0.12407160d0
 gh_guess_all[3]=-0.0014980603d0
 gh_guess_all[4]=-0.056740430d0
 gh_guess_all[5]=-0.0077393645d0
 gh_guess_all[6]=0.068103339d0
 gh_guess_all[7]=0.0048822754d0
 gh_guess_all[8]=-0.044190870d0
 gh_guess_all[9]=0.0079423653d0
 gh_guess_all[10]=-0.0021501101d0


gh_scale_all=dblarr(11)
gh_scale_all[0]=5d-1
gh_scale_all[1]=1d-2
gh_scale_all[2]=1d-2
gh_scale_all[3]=1d-2
gh_scale_all[4]=1d-2
gh_scale_all[5]=1d-2
gh_scale_all[6]=1d-2
gh_scale_all[7]=1d-2
gh_scale_all[8]=1d-2
gh_scale_all[9]=1d-2
gh_scale_all[10]=1d-2

;;;Gauss-Hermite parameters for amoeba
gh_lin_guess_all=dblarr(11)
; gh_lin_guess_all[0]=0.00030035155d0	;guesses before cb_hermite fix
; gh_lin_guess_all[1]=4.1945145d-05
; gh_lin_guess_all[2]=-4.6587297d-05
; gh_lin_guess_all[3]=-1.1758238d-05
; gh_lin_guess_all[4]=2.6893785d-06
; gh_lin_guess_all[5]=-5.5451940d-06
; gh_lin_guess_all[6]=6.9242824d-05
; gh_lin_guess_all[7]=-4.5808258d-05
; gh_lin_guess_all[8]=-1.8217831d-05
; gh_lin_guess_all[9]=-6.9966460d-05
; gh_lin_guess_all[10]=0.00010331374d0


gh_lin_guess_all[0]=-0.00020781873
gh_lin_guess_all[1]=-2.7751771e-06
gh_lin_guess_all[2]=6.6262665e-05
gh_lin_guess_all[3]=8.4012060e-06
gh_lin_guess_all[4]=-1.2985736e-06
gh_lin_guess_all[5]=-1.1491570e-06
gh_lin_guess_all[6]=-1.3797089e-06
gh_lin_guess_all[7]=1.6455046e-06
gh_lin_guess_all[8]=3.4848679e-06
gh_lin_guess_all[9]=-8.6695246e-06
gh_lin_guess_all[10]=6.1633282e-06

gh_lin_scale_all=dblarr(11)
gh_lin_scale_all[0]=1d-4
gh_lin_scale_all[1]=1d-6
gh_lin_scale_all[2]=1d-6
gh_lin_scale_all[3]=1d-6
gh_lin_scale_all[4]=1d-6
gh_lin_scale_all[5]=1d-6
gh_lin_scale_all[6]=1d-6
gh_lin_scale_all[7]=1d-6
gh_lin_scale_all[8]=1d-6
gh_lin_scale_all[9]=1d-6
gh_lin_scale_all[10]=1d-6

for run=0,niter-1 do begin
    
    ;;;Wavelength parameters for amoeba
    ;Now make arrays with proper number of parameters
    if n_elements(wl_guess) eq 0 then begin
        wl_guess=dblarr(nparam_wl)
        for i=0,nparam_wl-1 do wl_guess[i]=wl_guess_all[i]
    endif
    wl_scale=dblarr(nparam_wl)
    for i=0,nparam_wl-1 do wl_scale[i]=wl_scale_all[i]
    

    ;Shift guesses by a random amount
    if random eq 1 then begin
        random_shift_wl=randomn(seed, nparam_wl, /double)*(0.5d0)*wl_scale
        random_shift_wl[0]=0d0  ;don't shift 0th order guess
        wl_guess=wl_guess+random_shift_wl
        endif

    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;Gauss-Hermite parameters for amoeba

    if n_elements(gh_guess) eq 0 then begin
        gh_guess=dblarr(nparam_gh)
        for i=0,nparam_gh-1 do gh_guess[i]=gh_guess_all[i]
    endif
    gh_scale=dblarr(nparam_gh)
    for i=0,nparam_gh-1 do gh_scale[i]=gh_scale_all[i]
    

    ;Shift guesses by a random amount
    if random eq 1 then begin
        random_shift_gh=randomn(seed, nparam_gh, /double)*gh_scale*(0.5d0)
        gh_guess=gh_guess+random_shift_gh
    endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;Gauss-Hermite linear parameters for amoeba

if lin_switch eq 1 then begin
    if n_elements(gh_lin_guess) eq 0 then begin
        gh_lin_guess=dblarr(nparam_gh_lin)
        for i=0,nparam_gh_lin-1 do gh_lin_guess[i]=gh_lin_guess_all[i]
    endif
    gh_lin_scale=dblarr(nparam_gh)
    for i=0,nparam_gh-1 do gh_lin_scale[i]=gh_lin_scale_all[i]

    ;Shift guesses by a random amount
    if random eq 1 then begin
        random_shift_gh_lin=randomn(seed, nparam_gh, /double)*gh_lin_scale*(0.5d0)
        gh_lin_guess=gh_lin_guess+random_shift_gh_lin
    endif
endif else begin
    gh_lin_guess=[!values.f_nan]
    gh_lin_scale=[!values.f_nan]
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;Other parameters for amoeba
if n_elements(other_guess) eq 0 then begin
    if nparam_other le 13 then begin
    other_guess=dblarr(13)
    other_guess[0]=0.96073530d0
    other_guess[1]=0.99991994d0
    other_guess[2]=0.99942584d0
    other_guess[3]=0.99918029d0
    other_guess[4]=0.99746050d0
    other_guess[5]=0.99768481d0
    other_guess[6]=0.99627572d0
    other_guess[7]=0.99839309d0
    other_guess[8]=0.99862463d0
    other_guess[9]=0.99901706d0
    other_guess[10]=0.99859190d0
    other_guess[11]=0.99105320d0
    other_guess[12]=0.99016681d0

    other_guess=other_guess[0:nparam_other-1]
endif else begin
    other_guess=dblarr(nparam_other)
    for i=1,nparam_other-1 do other_guess[i]=1d0
    other_guess[0]=0.96073530d0  ;guess for tau
endelse

endif

other_scale=dblarr(nparam_other)
for i=0,nparam_other-1 do other_scale[i]=0.05d0
    

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



IRCS_cal_116, wl_guess, wl_scale, gh_guess, gh_scale, gh_lin_guess, gh_lin_scale, other_guess, other_scale, run, lin_switch=lin_switch, min_type=min_type, visualize=visualize

endfor

end
