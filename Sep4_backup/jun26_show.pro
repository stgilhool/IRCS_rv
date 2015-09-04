function r_stat, r, ind
s=r[ind, *]

nwin2=n_elements(s[*,0])
nexp=14L

weight=1/s.chi2

avg=mean(s.delta_rv, dimension=1)

med=median(s.delta_rv, dimension=1)

wavg=total(weight*s.delta_rv, 1)/(total(weight, 1))

;var=total((sqrt(s.chi2)/rebin(reform(total(weight,1), 1, nexp), nwin2, nexp))^2, 1)
var=(1d0/(nwin2-1))*total((s.delta_rv-rebin(reform(wavg, 1, nexp), nwin2, nexp))^2, 1)

sig=sqrt(var)

rms=sqrt(mean((s.delta_rv-rebin(reform(avg, 1, nexp), nwin2, nexp))^2, dimension=1))

big_weight=1/var
big_wavg=total((big_weight*wavg))/total(big_weight)
big_var=(1d0/nexp)*total((big_wavg-wavg)^2);total(((1/sig)/total(big_weight))^2)
big_sig=sqrt(big_var)

out={avg:avg, $
     med:med, $
     wavg:wavg, $
     var:var, $
     sig:sig, $
     rms:rms, $
     allrv:s.delta_rv, $
     rv:big_wavg, $
     big_var:big_var, $
     big_sig:big_sig $
     }
ymax=max(wavg+sig*1.1)
ymin=min(wavg-sig*1.1)
plot, wavg, ps=6, yr=[ymin, ymax], title="Weighted mean of RV results by exposure using " + strtrim(nwin2,2) + " chunks.  RV is " + strtrim(big_wavg, 2) + " +- " + strtrim(big_sig,2)
errplot, wavg-sig, wavg+sig

return, out

end

pro jun26_show

a=mrdfits('../data/rvshift1_results/Jun23_ac_manychunk_allexp.fits', 1)


nwin=21L
nexp=14L

t=replicate(a[0], nwin, nexp)

for exp=0,nexp-1 do t[*,exp]= a[exp*nwin:(exp+1)*nwin-1]


;ALL WINDOWS
allwi=lindgen(22)
allw=r_stat(t,allwi)

;First half
fhi=lindgen(10)
fh=r_stat(t,fhi)

;First half even
fhei=lindgen(5)*2
fhe=r_stat(t,fhei)

;First half odd
fhoi=lindgen(5)*2 + 1
fho=r_stat(t,fhoi)

;Second half
shi=lindgen(10)+12
sh=r_stat(t,shi)

;Second half even
shei=lindgen(5)*2 +12
she=r_stat(t,shei)

;Second half odd
shoi=lindgen(5)*2 +13
sho=r_stat(t,shoi)

;All but mask
alli=[fhi,shi]
all=r_stat(t,alli)

;All but mask even
allei=[fhei, shei]
alle=r_stat(t,allei)

;All but mask odd
alloi=[fhoi, shoi]
allo=r_stat(t,alloi)





stop 

end
