pro test_hermites


for i=0l, 999 do begin

xx=dindgen(201)-100.
sigma=12.

t=cb_hermite(0, xx, sigma)
;start with 0, which is a GAussian, with sigma "sigma"

;how many additional ones to include? Maybe just

nbasis=3l

n=randomn(seed, nbasis)*0.1

;start with order 1?
ostart=1l

for j=0l, nbasis-1 do t=t+cb_hermite(j+ostart, xx, sigma)*n[j]


plot, t, yr=[-.05, 0.25]

wait, .2

endfor

end
