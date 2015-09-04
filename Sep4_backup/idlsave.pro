; IDL Version 8.2.2 (linux x86_64 m64)
; Journal File for stgilhool@iroquois.physics.upenn.edu
; Working directory: /RAID/home/stgilhool/RV_projects/IRCS_rv/IDL_pros
; Date: Thu Jun 18 12:19:13 2015
 
.reset
template_comparison.pro
ot='../data/GJ273/GJ237_template.dat'
nt='../data/GJ273/GJ237_template_tapas.dat'
readcol, ot, wl_old, template_old, format='D,D'
; % FILE_LINES: Error opening file. File: ../data/GJ273/GJ237_template.dat
;   No such file or directory
ot='../data/GJ273/GJ273_template.dat'
nt='../data/GJ273/GJ273_template_tapas.dat'
readcol, ot, wl_old, template_old, format='D,D'
readcol, nt, wl_new, template_new, format='D,D'
plot, wl_old, template_old, /xs
plot, wl_old, template_old, /xs, yr=[0.5,1.2]
oplot, wl_new, template_new, color=200
otnorm=continuum_fit(wl_old, template_old, low_rej=0.53, high_rej=3.0)
template_oldn=template_old/otnorm
plot, wl_old, template_oldn, /xs, yr=[0.5,1.2]
oplot, wl_new, template_new, color=200
plot, wl_old, template_oldn, /xs, yr=[0.5,1.2], xr=[2.99,2.30]
plot, wl_old, template_oldn, /xs, yr=[0.5,1.2], xr=[2.29,2.30]
oplot, wl_new, template_new, color=200
plot, wl_old, template_oldn, /xs, yr=[0.5,1.2], xr=[2.30,2.31]
oplot, wl_new, template_new, color=200
plot, wl_old, template_oldn, /xs, yr=[0.5,1.2], xr=[2.31,2.32]
oplot, wl_new, template_new, color=200
