; IDL Version 8.2.2 (linux x86_64 m64)
; Journal File for stgilhool@iroquois.physics.upenn.edu
; Working directory: /RAID/home/stgilhool/RV_projects/IRCS_rv/IDL_pros
; Date: Sun Mar  8 22:02:32 2015
 
rv_t2=dblarr(14)
rv_r2=dblarr(14)
rv_r1=dblarr(14)
chi2_r1=dblarr(14)
chi2_r2=dblarr(14)
mjd=dblarr(14)

for i=0, 13 do begin & $
r1=mrdfits('/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_127_temp7_r1.fits', i+1) & $
r2=mrdfits('/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_127_temp7_r3.fits', i+1) & $
mjd[i]=r1.mjd & $
rv_r1[i]=r1.final_rv & $
rv_r2[i]=r2.final_rv & $
chi2_r1[i]=r1.chi2 & $
chi2_r2[i]=r2.chi2 & $
endfor


plot, mjd, rv_r2, ps=6, /xs
oplot, mjd, rv_r1, ps=6, color=200
good=[0,1,2,3,6,8,9,11,12,13]
;plot, mjd[good], rv_r2[good], ps=6
;plot, mjd[good], rv_r1[good], ps=6
print, max(rv_r2[good])-min(rv_r2[good])
;       216.63056
print, max(rv_r1[good])-min(rv_r1[good])
;       203.73069
print, chi2_r1[good]
;       2440.8364       2867.7920       1281.1978       1642.9250       921.47741       2467.1303       2157.9666       1408.6296       1229.1379       1346.4744
print, chi2_r2[good]
;       2346.2997       2656.0447       1168.9762       1478.0118       867.28117       2197.0645       1876.8996       1337.3171       1145.1890       1272.8544
stop
end
