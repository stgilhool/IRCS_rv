function r_stat, r, ind
s=r[ind, *]

nwin2=n_elements(s[*,0])
nexp=14L

weight=1/s.chi2

avg=mean(s.final_rv, dimension=1)

med=median(s.final_rv, dimension=1)

wavg=total(weight*s.final_rv, 1)/(total(weight, 1))

;var=total((sqrt(s.chi2)/rebin(reform(total(weight,1), 1, nexp), nwin2, nexp))^2, 1)
var=(1d0/(nwin2-1))*total((s.final_rv-rebin(reform(wavg, 1, nexp), nwin2, nexp))^2, 1)

sig=sqrt(var)

rms=sqrt(mean((s.final_rv-rebin(reform(avg, 1, nexp), nwin2, nexp))^2, dimension=1))

big_weight=1/var
big_wavg=total((big_weight*wavg))/total(big_weight)
big_var=(1d0/(nexp-1))*total((big_wavg-wavg)^2);total(((1/sig)/total(big_weight))^2)
big_sig=sqrt(big_var)

out={avg:avg, $
     med:med, $
     wavg:wavg, $
     var:var, $
     sig:sig, $
     rms:rms, $
     allrv:s.final_rv, $
     rv:big_wavg, $
     big_var:big_var, $
     big_sig:big_sig $
     }
ymax=max(wavg+sig*1.1)
ymin=min(wavg-sig*1.1)
plot, wavg, ps=6, yr=[ymin, ymax], xtitle='Exposure number', ytitle='RV (km/s)', title="Weighted mean of RV results by exposure using " + strtrim(nwin2,2) + " chunks.  RV is " + strtrim(big_wavg, 2) + " +- " + strtrim(big_sig,2)
errplot, wavg-sig, wavg+sig


return, out

end


pro plot_exp_rv, expnum, stopopt=stopopt, ovplot=ovplot
if n_elements(stopopt) eq 0 then stopopt=0
if n_elements(ovplot) eq 0 then ovplot=0

common data, t
tt=t
exprv=tt[*,expnum].final_rv
meanrv=mean(exprv)
rvrms=sqrt(mean((exprv-meanrv)^2))

winchi2=tt[*,expnum].chi2
chi2str=strtrim(winchi2,2)

nwind=n_elements(exprv)
x=lindgen(nwind)+1

if ovplot then begin
    oplot,x, exprv, ps=6, color=200
    ovplot=0
endif else begin
    plot,x, exprv, ps=6, xtitle='window number', ytitle='RV km/s', title='RV scatter of: ' + strtrim(rvrms,2) + ' km/s for exposure number: ' + strtrim(expnum,2), xr=[0,nwind+1]
endelse


xyouts, x, exprv*1.05d0, chi2str

if stopopt then begin
    stopopt=0
    stop
endif

end


pro jul2_show, model_tag

if n_elements(model_tag) eq 0 then model_tag='Jul1ptemp1'

;a=mrdfits('../data/rvshift1_results/Jun29_ac_manychunk_allexp.fits',1)
;a=mrdfits('../data/rvshift1_results/Jun29_ac_manychunk_smooth.fits',1)

;a=mrdfits('../data/rvshift1_results/Jun30_ac_manychunk_phoenix.fits',1)
;a=mrdfits('../data/rvshift1_results/Jun30_ac_manychunk_phoenixsmooth.fits',1)
a=mrdfits('../data/rvshift1_results/rvfit_'+model_tag+'_allexp.fits',1)

;fid='../data/rvshift1_results/Jun29_ac_manychunk_allexp_temp_1.fits'
;fid='../data/rvshift1_results/Jun29_ac_manychunk_smooth_temp_0.fits'
;fid='../data/rvshift1_results/Jun30_ac_manychunk_phoenix_temp_0.fits'
fid='../data/rvshift1_results/rvfit_'+model_tag+'_0.fits'


fits_info, fid, n_ext=nwin, /silent

nexp=14L

common data, t

t=replicate(a[0], nwin, nexp)

for exp=0,nexp-1 do t[*,exp]= a[exp*nwin:(exp+1)*nwin-1]


;ALL WINDOWS
allwi=lindgen(nwin)
allw=r_stat(t,allwi)

;First half
fhi=lindgen(nwin<4)
fh=r_stat(t,fhi)

; First half even
; fhei=lindgen(5)*2
; fhe=r_stat(t,fhei)

; First half odd
; fhoi=lindgen(5)*2 + 1
; fho=r_stat(t,fhoi)

; Second half
; shi=lindgen(10)+12
; sh=r_stat(t,shi)

; Second half even
; shei=lindgen(5)*2 +12
; she=r_stat(t,shei)

; Second half odd
; shoi=lindgen(5)*2 +13
; sho=r_stat(t,shoi)

; All but mask
; alli=[fhi,shi]
; all=r_stat(t,alli)

; All but mask even
; allei=[fhei, shei]
; alle=r_stat(t,allei)

; All but mask odd
; alloi=[fhoi, shoi]
; allo=r_stat(t,alloi)





stop 

end
