pro temp

f1='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/rvfit_Jul10ksfit1_0.fits'
f2='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/rvfit_Jul10ksfit1_1.fits'

s1=mrdfits(f1,1)
o=replicate(s1, 5)
s2=mrdfits(f2,1)
t=replicate(s2, 5)

for i=0, 4 do begin
    o[i]=mrdfits(f1,i+1)
    t[i]=mrdfits(f2,i+1)
endfor

orv=o.final_rv
omjd=o.mjd

trv=t.final_rv
tmjd=t.mjd

stop

end
