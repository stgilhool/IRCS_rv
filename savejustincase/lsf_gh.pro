function lsf_gh, gh0_vec, gh1_vec=gh1_vec, oversamp=oversamp, first_pix=first_pix, npix_select=npix_select, verbose=verbose, decomp=decomp

;Idea for improvement: pass H_coeff in with function call?
npix=1024L ;CAUTION this needs to be tuned for non-IRCS data sets

;Deal with inputs
nparam_gh=n_elements(gh0_vec)
;Set verbose
if n_elements(verbose) eq 0 then verbose=0
;Set decomp
if n_elements(decomp) eq 0 then decomp=0
;Sigma/gh0
sigma=gh0_vec[0]
gh0_vec[0]=1d0
;sigma_lin/gh1
;check if lin switch is on
if n_elements(gh1_vec) eq 0 then begin
    gh1_vec=replicate(0d0, nparam_gh)
    lin_switch=0
endif else begin
    NA=where(gh1_vec ne 0, gh1count)
    if gh1count eq 0 then lin_switch=0 else lin_switch=1
endelse
;set sigma_lin
if lin_switch eq 0 then sigma_lin=0d0 else begin
    sigma_lin=gh1_vec[0]
    gh1_vec[0]=0d0
endelse

;Set oversamp/first_pix/npix_select
if n_elements(oversamp) eq 0 then oversamp=7L            
if n_elements(first_pix) eq 0 then first_pix=0L
if n_elements(npix_select) eq 0 then npix_select=npix-first_pix

;set npix_model
npix_model=npix_select*oversamp
;set npix_lsf
npix_lsf=(oversamp*10L)+1L

;start clock
lsf_tic=tic()

;populate lsf coordinates
x_lsf=rebin((dindgen(npix_lsf)-(npix_lsf/2))/oversamp, npix_lsf, npix_model)
pixel_2d=rebin(reform((((dindgen(npix_model)-(oversamp/2))/oversamp)+first_pix), 1, npix_model), npix_lsf, npix_model)
sig_arr=sigma+(sigma_lin*pixel_2d)
x_sig_2d=x_lsf/sig_arr
x_sig_3d=rebin(x_sig_2d, npix_lsf, npix_model, nparam_gh)

;Make non-normalized gh-polynomial cube (npix_lsf x npix_model x
;basis)

H_coeff=sg_hermite_coeff()
gh_cube_nonorm=dblarr(npix_lsf, npix_model, nparam_gh)
gh_cube_norm=dblarr(npix_lsf, npix_model, nparam_gh)
gh_cube=dblarr(npix_lsf, npix_model, nparam_gh)
basis_coeff=dblarr(npix_lsf, npix_model, nparam_gh)
for basis=0, nparam_gh-1 do begin
    ;make gh_polynomial for given basis
    gh_cube_nonorm[*,*,basis]=poly(x_sig_2d, H_coeff[*,basis])*exp(-0.5*(x_sig_2d^2))
    ;normalize each basis
    norm=rebin(reform(total(abs(gh_cube_nonorm[*,*,basis]), 1, /double), 1, npix_model), npix_lsf, npix_model)
    gh_cube_norm[*,*,basis]=gh_cube_nonorm[*,*,basis]/norm
    ;Adjust amplitudes according to gh parameters
    if lin_switch eq 1 then basis_coeff[*,*,basis]=gh0_vec[basis]+(gh1_vec[basis]*pixel_2d) $
      else if lin_switch eq 0 then basis_coeff[*,*,basis]=replicate(gh0_vec[basis], npix_lsf, npix_model)
    gh_cube[*,*,basis]=basis_coeff[*,*,basis]*gh_cube_norm[*,*,basis]
endfor

;Sum the bases
if nparam_gh gt 1 then lsf_nonorm=total(gh_cube, 3, /double) $
  else lsf_nonorm=gh_cube[*,*,0]

;account for unrealistically negative lsfs

min_vec=min(lsf_nonorm, dimension=1)
neg_ind=where(min_vec lt 0, neg_count)
if neg_count gt 0 then begin
    neg_vec=dblarr(npix_model)
    neg_vec[neg_ind]=min_vec[neg_ind]

    neg_arr=rebin(reform(neg_vec, 1, npix_model), npix_lsf, npix_model)

    lsf_nonorm_copy=lsf_nonorm
    lsf_nonorm=lsf_nonorm-neg_arr
endif

;normalize lsf
lsf_norm=rebin(reform(total(lsf_nonorm, 1, /double), 1, npix_model), npix_lsf, npix_model)
lsf=lsf_nonorm/lsf_norm
if lin_switch eq 0 then lsf=lsf[*,0]


;;;;
;Fix LSF center of mass
;;;;
lsf_com=rebin(reform(total(lsf*x_sig_2d, 1, /double), 1, npix_model), npix_lsf, npix_model)
x_sig_2d_com=x_sig_2d+lsf_com
x_sig_3d_com=rebin(x_sig_2d_com, npix_lsf, npix_model, nparam_gh)

;;;Remake LSF
;Make non-normalized gh-polynomial cube (npix_lsf x npix_model x
;basis)

gh_cube_nonorm_com=dblarr(npix_lsf, npix_model, nparam_gh)
gh_cube_norm_com=dblarr(npix_lsf, npix_model, nparam_gh)
gh_cube_com=dblarr(npix_lsf, npix_model, nparam_gh)
for basis=0, nparam_gh-1 do begin
    ;make gh_polynomial for given basis
    gh_cube_nonorm_com[*,*,basis]=poly(x_sig_2d_com, H_coeff[*,basis])*exp(-0.5*(x_sig_2d_com^2))
    ;normalize each basis
    norm_com=rebin(reform(total(abs(gh_cube_nonorm_com[*,*,basis]), 1, /double), 1, npix_model), npix_lsf, npix_model)
    gh_cube_norm_com[*,*,basis]=gh_cube_nonorm_com[*,*,basis]/norm_com
    
    gh_cube_com[*,*,basis]=basis_coeff[*,*,basis]*gh_cube_norm_com[*,*,basis]
endfor

;Sum the bases
if nparam_gh gt 1 then lsf_nonorm_com=total(gh_cube_com, 3, /double) $
  else lsf_nonorm_com=gh_cube_com[*,*,0]

;account for unrealistically negative lsfs

min_vec_com=min(lsf_nonorm_com, dimension=1)
neg_ind_com=where(min_vec_com lt 0, neg_count_com)
if neg_count_com gt 0 then begin
    neg_vec_com=dblarr(npix_model)
    neg_vec_com[neg_ind_com]=min_vec_com[neg_ind_com]

    neg_arr_com=rebin(reform(neg_vec_com, 1, npix_model), npix_lsf, npix_model)

    lsf_nonorm_com_copy=lsf_nonorm_com
    lsf_nonorm_com=lsf_nonorm_com-neg_arr_com
endif

;normalize lsf
lsf_norm_com=rebin(reform(total(lsf_nonorm_com, 1, /double), 1, npix_model), npix_lsf, npix_model)
lsf_com=lsf_nonorm_com/lsf_norm_com
if lin_switch eq 0 then lsf_com=lsf_com[*,0]

;make copy of old lsf 
;lsf_copy=lsf
lsf=lsf_com

;Decompose the LSF, if necessary
if decomp then begin
    if neg_count_com gt 0 then begin
        gh_decomp_pos=gh_cube_com-rebin(neg_arr_com, npix_lsf, npix_model, nparam_gh)
    endif else gh_decomp_pos=gh_cube_com
    ;normalize
    decomp_norm=rebin(lsf_norm_com, npix_lsf, npix_model, nparam_gh)
    gh_decomp=gh_decomp_pos/decomp_norm

    lsf_t=toc(lsf_tic)
    if verbose then begin
        print, "LSF construction took : ", lsf_t, " seconds."
        print, "Returning decomposed LSF"
    endif
        ;lsf=gh_decomp
endif else begin

;Stop clock and report if verbose
    lsf_t=toc(lsf_tic)
    if verbose then begin
        print, "LSF construction took : ", lsf_t, " seconds."
        print, "Returning LSF"
    endif

endelse
if decomp then return, gh_decomp else return, lsf
end
