pro IRCS_make_figs, chi2_dof, trial_soln, trial_lsf, oversamp, npix_lsf, int_obs, int_model, npix_trim_start, npix_trim_end, lin_switch

npix=1024L
;get rms
residuals=int_obs-int_model

rmse=sqrt(mean(residuals[npix_trim_start:-1*npix_trim_end+1]^2))


title_str="NH3 Model Fitting Process | chi2/degree of freedom : "+strtrim(chi2_dof,2)
title_str_tex="$\textrm{NH}_3 $ Model Fitting Process $\vert$ $\chi^2/\textrm{degree of freedom}$ : "+strtrim(chi2_dof,2)
phi=findgen(32)*(!PI*2/32.)
phi = [ phi, phi(0) ]
usersym, cos(phi), sin(phi), /fill



set_plot, 'ps'
device, /encapsulated, filename = 'fit_temp.eps'
!p.font=0
device, /color, bits=8
loadct, 12


!p.multi=0


plot, trial_soln, int_obs, /xs, title=title_str,xtitle="Wavelength (microns)", thick=5.0, xr=[2.31,2.33], yr=[0.4,1.1];yr=[0.35,1.15], /ys, yticks=16, ytickname=[' ',' ','-0.1',' ','0.0',' ','0.1',' ',' ', '0.8',' ','0.9',' ', '1.0',' ', '1.1', ' ']

;oplot, trial_soln[dindgen(npix/3)*3], int_model[dindgen(npix/3)*3], ps=8, symsize=0.6, color=200 
oplot, trial_soln, int_model, ps=8, symsize=0.4, color=200 
oplot, trial_soln, residuals+0.6, ps=8, symsize=0.2
al_legend, ['NH3 observed', 'NH3 model', 'Residuals'], psym=[-0,88,88], color=['black','red','black']
xyouts,2.32, 0.56, 'RMS : '+strtrim(rmse,2)
;axgap, [0,0.7]
;plot, [trial_soln[0],trial_soln[-1]],[-0.1d0,0.1d0], position=[0, 0, 1, 0.3], xstyle=5, ystyle=17, /nodata, /noerase

device, /close
latexify, 'fit_temp.eps', [title_str, 'NH3 observed', 'NH3 model','Residuals', 'Wavelength (microns)'], [title_str_tex, '$\textrm{NH}_3 $ obs', '$\textrm{NH}_3 $ model', 'Residuals', '$\lambda (\mu m)$'], /full


x1=100L
x2=500L
x3=900L
;x2=250L
;x3=450L

xx1=x1*oversamp+(oversamp/2)
xx2=x2*oversamp+(oversamp/2)
xx3=x3*oversamp+(oversamp/2)
if lin_switch eq 1 then begin
    lsf1=trial_lsf[*,xx1]
    lsf2=trial_lsf[*,xx2]
    lsf3=trial_lsf[*,xx3]
    
endif else begin
    lsf1=trial_lsf
    lsf2=trial_lsf
    lsf3=trial_lsf
endelse

max_lsf=max(lsf1)>max(lsf2)>max(lsf3)




!p.multi=0
device, /encapsulated, filename =  'lsf_temp.eps'
!p.font=0
device, /color, bits=8
loadct, 12


if lin_switch eq 1 then begin
    plot, dindgen(npix_lsf)/oversamp-npix_lsf/(oversamp*2), lsf1, xtitle='Pixel (relative to x_o)', title="Line Spread Function", /xs
    oplot, dindgen(npix_lsf)/oversamp-npix_lsf/(oversamp*2), lsf2, linestyle=2
    oplot, dindgen(npix_lsf)/oversamp-npix_lsf/(oversamp*2), lsf3, linestyle=4
    al_legend, ['LSF xxxxx 100', 'LSF xxxxx 500', 'LSF xxxxx 900'], linestyle=[0,2,4]
endif else begin
    plot, dindgen(npix_lsf), trial_lsf
endelse
device, /close
latexify, 'lsf_temp.eps',['Pixel (relative to x_o)','LSF xxxxx 100', 'LSF xxxxx 500', 'LSF xxxxx 900'], ['Pixel (relative to $x_o$)','$\textrm{LSF}_{x_o = 100}$','$\textrm{LSF}_{x_o = 500}$','$\textrm{LSF}_{x_o = 900}$'], /full

end
