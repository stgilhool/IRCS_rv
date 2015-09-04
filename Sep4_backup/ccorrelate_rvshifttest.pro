; IDL Version 8.2.2 (linux x86_64 m64)
; Journal File for stgilhool@iroquois.physics.upenn.edu
; Working directory: /RAID/home/stgilhool/RV_projects/IRCS_rv/IDL_pros
; Date: Mon Jul 13 10:01:43 2015
 
;s1 is 11 times oversampled first 128 pixels of spectrum from exp 0
;s2 is the same thing for exp 1
;w is rough wl scale
;testing to see what best rv shift is
ntrials=100L
totres=dblarr(ntrials)
ntrials=101L
totres=dblarr(ntrials)
rvshift=dindgen(ntrials)-(ntrials/2)
rvshift=rvshift*20d0
bigc=3d8
for try=0,ntrials-1 do begin & $
rv=rvshift[try] & $
wlshift=w*(1d0 + (rv/bigc)) & $
s1shift=interpol(s1, wlshift, w) & $
res=s1shift-s2 & $
totres[try]=total(res^2, /double) & $
endfor
print, totres
;      0.39347971      0.39335275      0.39337847      0.39355687      0.39388794      0.39437169      0.39500812      0.39579722      0.39673900      0.39783346      0.39908060      0.40048041      0.40203290      0.40373806      0.40559590      0.407606
;42      0.40976962      0.41208549
;      0.41455404      0.41717527      0.41994917      0.42300727      0.42626628      0.42967785      0.43324196      0.43695863      0.44082785      0.44484963      0.44902396      0.45335084      0.45783027      0.46246226      0.46724680      0.472183
;89      0.47727353      0.48251573
;      0.48791048      0.49345778      0.49915763      0.50501003      0.51101499      0.51717250      0.52348256      0.52994517      0.53656034      0.54332806      0.55024833      0.55732115      0.56454652      0.57192444      0.57945492      0.587547
;84      0.59579319      0.60419097
;      0.61274118      0.62144382      0.63029889      0.63930638      0.64846630      0.65777865      0.66724343      0.67686063      0.68663027      0.69655233      0.70662681      0.71685373      0.72723307      0.73776484      0.74844904      0.759285
;66      0.77027472      0.78141620
;      0.79271010      0.80415644      0.81575520      0.82750638      0.83941000      0.85146604      0.86367451      0.87603540      0.88854594      0.90120479      0.91401595      0.92697940      0.94009517      0.95336323      0.96678360      0.980356
;27      0.99408124       1.0079585
;       1.0219881       1.0361700       1.0505042       1.0649907       1.0796295       1.0944205       1.1093639       1.1244596       1.1397076       1.1551080       1.1706606
plot, totres, ps=6
ntrials=101L
totres=dblarr(ntrials)
rvshift=dindgen(ntrials)-(ntrials/3)
rvshift=dindgen(ntrials)-(2*ntrials/3)
rvshift=rvshift*20d0
for try=0,ntrials-1 do begin & $
rv=rvshift[try] & $
wlshift=w*(1d0 + (rv/bigc)) & $
s1shift=interpol(s1, wlshift, w) & $
res=s1shift-s2 & $
totres[try]=total(res^2, /double) & $
endfor
plot, rvshift, totres, ps=6
print, rvshift[where(totres eq min(totres))]
;      -980.00000
