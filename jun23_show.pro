pro jun23_show
;Just a quick batch file to show recent results
f='../data/rvshift1_results/Jun22_ac_manychunk_gaussall.fits'

s=mrdfits(f, 1)



window, 0
;p=plot(s.delta_rv, ps=6)
plot, s.delta_rv, ps=6, title='RV for different windows', xtitle='starting pixel', ytitle='RV (km/s)'

;xaxis=axis(
window, 1
plot, s.first_pix, s.result[1], ps=6, title='H2O depth for each window', xtitle='starting pixel', ytitle='Tau H2O'



stop

window, 2
!p.multi=[0,1,2]
for i=18, 35 do begin
    plot, s[i].obs_select, title='Obs for '+strtrim(s[i].first_pix,2)+ ' to ' + strtrim(s[i].first_pix+s[i].npix_select-1,2)+ ' | TauH2O : ' + strtrim(s[i].result[1], 2) + ' | RV : ' + strtrim(s[i].delta_rv,2)
    oplot, s[i].model_select, ps=6, color=200
    plot, s[i].obs_select-s[i].model_select, ps=3, yr=[-0.03, 0.03]
    stop
    ;wait, 1
endfor
!p.multi=0

stop
end
