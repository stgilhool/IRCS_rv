function parscale_rescale, parinfo_str, llim=llim, ulim=ulim

npar = n_elements(parinfo_str)

;SCALE PARAMETERS TO BE BETWEEN 0 and 1, if limits are defined
if n_elements(ulim) eq 0 then ulim=replicate(1d0, npar) else $
  if n_elements(ulim) eq 1 then ulim=replicate(ulim, npar)
if n_elements(llim) eq 0 then llim=replicate(0d0, npar) else $
  if n_elements(llim) eq 1 then llim=replicate(llim, npar)


;read in true parameter value and limits
t_llim=parinfo_str.tlimits[0]
t_ulim=parinfo_str.tlimits[1]
t_val=parinfo_str.tvalue
t_step=parinfo_str.tstep
t_mpmaxstep=parinfo_str.tmpmaxstep

;check limits
highindex = where(t_val gt t_ulim, highcount)
lowindex = where(t_val lt t_llim, lowcount)

if highcount gt 0 then begin
    
    print, "value exceeds upper limit for parameter(s): " + parinfo_str[highindex].parname
    print, "changing upper limit..."

    print, t_ulim[highindex]

    diff = t_val[highindex] - t_ulim[highindex]
    t_ulim[highindex] = t_ulim[highindex] + (diff * 1.1d0)
    parinfo_str.tlimits[1] = t_ulim


    print, 'to'
    print, t_ulim[highindex]
    print, '======'
    print, 'Upper: ', t_ulim[highindex]
    print, 'Value: ', t_val[highindex]
    print, 'Lower: ', t_llim[highindex]
    
    wait, 3
endif

if lowcount gt 0 then begin
    
    print, "value below lower limit for parameter(s): " + parinfo_str[lowindex].parname
    print, "changing lower limit..."
    
    print, t_llim[lowindex]
    
    diff = t_llim[lowindex] - t_val[lowindex]
    t_llim[lowindex] = t_llim[lowindex] - (diff * 1.1d0)
    parinfo_str.tlimits[0] = t_llim
    
    print, 'to'
    print, t_llim[lowindex]
    print, '======'
    print, 'Upper: ', t_ulim[lowindex]
    print, 'Value: ', t_val[lowindex]
    print, 'Lower: ', t_llim[lowindex]

    wait, 3
endif




;define scale factors
pscale=1d0/(t_ulim-t_llim)
pshift=llim-(pscale*t_llim)

parinfo_str.parscale=transpose([[pscale],[pshift]])

;Scale the value and limits
parinfo_str.limits=transpose([[llim], [ulim]])

parinfo_str.value = (pscale * t_val) + pshift
parinfo_str.step = (pscale * t_step)
parinfo_str.mpmaxstep = (pscale * t_mpmaxstep)


return, parinfo_str

end


function param_scale_down, t_p, pscale, pshift

p = (pscale*t_p) + pshift

return, p

end


function param_scale_up, p, pscale, pshift

t_p = (p-pshift)/pscale

return, t_p

end


function parinfo_getindex, parinfo_str, all=all, used=used, free=free
;all setting finds indices for entire parinfo_str
;free finds indices for params with fixed = 0
;used finds indices for params with used = 1

if n_elements(all) eq 0 then all = 0

if n_elements(used) eq 0 then used = 0
if n_elements(free) eq 0 then free = 0

;ensure that only 1 option is set, and assign default if none is
if total([all,free,used]) eq 0 then begin
    print, "Setting free as default setting"
    free = 1
endif else if total([all,free,used]) gt 1 then message, 'choose all or used or free'
;


if all eq 0 then begin
    if used eq 1 then ind = where(parinfo_str.used eq 1, parcount)$
      else if free eq 1 then ind = where(parinfo_str.fixed eq 0, parcount)
    if parcount eq 0 then message, "error: must have used params"
endif else ind = lindgen(n_elements(parinfo_str))

parnames= parinfo_str[ind].parname

                                ;;Get rv info
rvi	= where(strmatch(parnames, 'd_rv *') eq 1, n_rv)

                                ;;Get wl vectors
wli	= where(strmatch(parnames, 'wl[0-4]*') eq 1, n_wl)

d_wli	= where(strmatch(parnames, 'd_wl[0-2]*') eq 1, n_d_wl)

                                ;;Get constant lsf params
ghi	= where(strmatch(parnames, 'gh[0-9]*') eq 1, n_gh)

                                ;;Get lsf shift params
d_ghi	= where(strmatch(parnames, 'd_gh[0-9]*') eq 1, n_d_gh)

                                ;;Get tau params
ti	= where(strmatch(parnames, 't\_*') eq 1, n_t)

                                ;;Get norm params
ki	= where(strmatch(parnames, 'k[0-9]*') eq 1, n_k)

outstr 	= {rvi:rvi, wli:wli, d_wli:d_wli, ghi:ghi, d_ghi:d_ghi, ti:ti, ki:ki, $
           n_rv:n_rv, n_wl:n_wl, n_d_wl:n_d_wl, n_gh:n_gh, n_d_gh:n_d_gh, $
           n_t:n_t, n_k:n_k}

return, outstr

end

function parinfo_changeval, parinfo_str, vals, parname=parname, index=index

parnames=strtrim(parinfo_str.parname,2)

if n_elements(index) eq 0 then begin
    
    npar  = n_elements(parname)
    index = lonarr(npar)
    
    for param_ind = 0, npar - 1 do begin
        index[param_ind] = where(parnames eq parname[param_ind], icount)
    endfor
endif else npar = n_elements(index)


;Change values


parinfo_str[index].tvalue = vals

parinfo_str = parscale_rescale(parinfo_str)

return, parinfo_str

end


function parinfo_changestep, parinfo_str, steps, parname=parname, index=index

parnames=strtrim(parinfo_str.parname,2)

if n_elements(index) eq 0 then begin
    
    npar  = n_elements(parname)
    index = lonarr(npar)
    
    for param_ind = 0, npar - 1 do begin
        index[param_ind] = where(parnames eq parname[param_ind], icount)
    endfor
endif else npar = n_elements(index)


;Change values


parinfo_str[index].tstep = steps

parinfo_str = parscale_rescale(parinfo_str)

return, parinfo_str

end




function parinfo_changelim, parinfo_str, lims, parname=parname, index=index, lower=lower, upper=upper, auto=auto

parnames=strtrim(parinfo_str.parname,2)

if n_elements(lower) ne 0 and n_elements(upper) ne 0 then begin
    lower = 0
    upper = 0
endif

if n_elements(lower) eq 0 then lower = 0
if n_elements(upper) eq 0 then upper = 0

if n_elements(index) eq 0 then begin
    
    npar  = n_elements(parname)
    index = lonarr(npar)
    
    for param_ind = 0, npar - 1 do begin
        index[param_ind] = where(parnames eq parname[param_ind], icount)
    endfor
endif else npar = n_elements(index)



;Change values
;check if using auto limits
if n_elements(auto) ne 0 then begin
    t_val = parinfo_str[index].tvalue
    t_delta = abs(t_val) * auto
    t_llim = t_val - t_delta
    t_ulim = t_val + t_delta
    if lower eq 1 then lims = t_llim else $
      if upper eq 1 then lims = t_ulim else $
      lims = transpose([[t_llim], [t_ulim]])
endif 
    
if lower eq 1 then parinfo_str[index].tlimits[*,0] = lims else $
  if upper eq 1 then parinfo_str[index].tlimits[0,*] = lims else $
  parinfo_str[index].tlimits = lims

parinfo_str = parscale_rescale(parinfo_str)

return, parinfo_str

end


function parinfo_randomval, parinfo_str, parname=parname, index=index, ranfactor=ranfactor

parnames=strtrim(parinfo_str.parname,2)

if n_elements(index) eq 0 then begin
    
    npar  = n_elements(parname)
    index = lonarr(npar)
    
    for param_ind = 0, npar - 1 do begin
        index[param_ind] = where(parnames eq parname[param_ind], icount)
    endfor
endif else npar = n_elements(index)

if n_elements(ranfactor) eq 0 then ranfactor=0.1d0

;Change values
start_tval = parinfo_str[index].tvalue

new_start_tval = start_tval * (1d0 + ((2d0*randomu(seed, npar, /double)-1d0) * ranfactor))

parinfo_str[index].tvalue = new_start_tval

parinfo_str = parscale_rescale(parinfo_str)

return, parinfo_str

end



function parinfo_readin, model_id, rawform=rawform

if n_elements(rawform) eq 0 then rawform = 0

path 		= ircsrv_paths()
model_path	= path.calib
parfile		= model_path+model_id

pinfo		= mrdfits(parfile, 2)

if rawform eq 0 then begin
    pinfo.tvalue = pinfo.tbestfit
    pinfo = parscale_rescale(pinfo)
endif

return, pinfo

end


pro parinfo_all_freepar, parinfo_all, n_rv=n_rv, n_wl=n_wl, n_d_wl=n_d_wl, n_gh=n_gh, n_d_gh=n_d_gh, n_t=n_t, n_k=n_k

if n_elements(n_rv) eq 0 then n_rv = 0
if n_elements(n_wl) eq 0 then n_wl = 0
if n_elements(n_d_wl) eq 0 then n_d_wl = 0
if n_elements(n_gh) eq 0 then n_gh = 0
if n_elements(n_d_gh) eq 0 then n_d_gh = 0
if n_elements(n_t) eq 0 then n_t = 0
if n_elements(n_k) eq 0 then n_k = 0

rvi=0
wli=lindgen(5)+1
d_wli=lindgen(3)+6
ghi=lindgen(11)+9
d_ghi=lindgen(11)+20
ti=lindgen(5)+31
ki=lindgen(12)+36

nfree = n_rv + n_wl + n_d_wl + n_gh + n_d_gh + n_t + n_k


;CREATE MASTER PARINFO STRUCTURE
;parinfo_all = parinfo_all_init()

;First, turn off all params
parinfo_all.fixed	= 1
parinfo_all.used	= 0

;TURN ON FREE PARAMETERS
if n_rv gt 0 then begin
    if n_rv gt 1 then message,"Only 1 RV parameter allowed"
    parinfo_all[rvi].fixed = 0
    parinfo_all[rvi].used = 1
endif

if n_wl gt 0 then begin
    if n_wl gt 5 then message, "Only 5 WL parameters allowed"
    parinfo_all[wli[0:n_wl-1]].fixed 	= 0
    parinfo_all[wli[0:n_wl-1]].used 	= 1
endif

if n_d_wl gt 0 then begin
    if n_d_wl gt 3 then message, "Only 3 D_WL parameters allowed"
    parinfo_all[d_wli[0:n_d_wl-1]].fixed 	= 0
    parinfo_all[d_wli[0:n_d_wl-1]].used 	= 1
endif

if n_gh gt 0 then begin
    if n_gh gt 11 then message, "Only 11 GH parameters allowed"
    parinfo_all[ghi[0:n_gh-1]].fixed 	= 0
    parinfo_all[ghi[0:n_gh-1]].used 	= 1
endif

if n_d_gh gt 0 then begin
    if n_d_gh gt 11 then message, "Only 11 D_GH parameters allowed"
    parinfo_all[d_ghi[0:n_d_gh-1]].fixed 	= 0
    parinfo_all[d_ghi[0:n_d_gh-1]].used 	= 1
endif

if n_t gt 0 then begin
    if n_t gt 5 then message, "Only 5 TAU parameters allowed"
    parinfo_all[ti[0:n_t-1]].fixed 	= 0
    parinfo_all[ti[0:n_t-1]].used 	= 1
endif

if n_k gt 0 then begin
    if n_k gt 12 then message, "Only 12 K parameters allowed"
    parinfo_all[ki[0:n_k-1]].fixed 	= 0
    parinfo_all[ki[0:n_k-1]].used 	= 1
endif



end


pro parinfo_all_init, parinfo_output

nparam_max=48L

parstr={parname:'null', $
        value:0d0, $
        step:0d0, $
        relstep:0d0, $
        fixed:1, $
        limited:[0,0], $
        limits:[0d0, 0d0], $
        mpside:0, $
        mpmaxstep:0, $
        mpprint:1, $
        bestfit:!values.d_nan, $
        tbestfit:!values.d_nan, $
        tvalue:0d0, $
        tlimits:[0d0, 0d0], $
        tmpmaxstep:0d0, $
        tstep:0d0, $
        parscale:[1d0,0d0], $
        used:0, $
        convflag:0, $
        fmode:'null' $
       }


parinfo_all = replicate(parstr, nparam_max)

rvi=0
wli=lindgen(5)+1
d_wli=lindgen(3)+6
ghi=lindgen(11)+9
d_ghi=lindgen(11)+20
ti=lindgen(5)+31
ki=lindgen(12)+36


parnames=['d_rv', $
          'wl0','wl1','wl2','wl3','wl4', $
          'd_wl0','d_wl1','d_wl2', $
          'gh0','gh1','gh2','gh3','gh4','gh5','gh6','gh7','gh8','gh9','gh10', $
          'd_gh0','d_gh1','d_gh2','d_gh3','d_gh4','d_gh5','d_gh6','d_gh7','d_gh8','d_gh9','d_gh10', $
          't_nh3','t_h2o','t_co2ch4','t_co2','t_ch4', $
          'k0','k1','k2','k3','k4','k5','k6','k7','k8','k9','k10','k11' $
          ]

parinfo_all.parname=parnames

;VALUE
values=dblarr(nparam_max)

values[rvi]   	=	[0d0]
values[wli]	=	[2.29107d0, 6.255d-5, -2.046d-9, -1.55d-12, 0d0]
values[d_wli]	=	[0d0, 0d0, 0d0]
values[ghi]	=	[0.7d0,0.05d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0]
values[d_ghi]	=	[0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0]
values[ti]	=	[0.95d0,0.3d0,0.3d0,0d0,0d0]
values[ki]	=	[1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0]

parinfo_all.tvalue=values

;STEP
steps=dblarr(nparam_max)          

;steps[rvi]	=	[1d-1]
;steps[wli]	=	[1d-1, 1d-6, 0d0, 0d0, 0d0]
;steps[d_wli]	=	[1d-6, 0d0, 0d0]
;steps[ghi]	=	[0.1d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0]
;steps[d_ghi]	=	[0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0]
;steps[ti]       =	[1d-2, 1d-1, 1d-1, 0d0, 0d0]
;steps[ki]	=	[1d-1,1d-1,1d-1,1d-1,1d-1,1d-1,1d-1,1d-1,1d-1,1d-1,1d-1,1d-1]

;parinfo_all.step=steps

;LIMITED
parinfo_all[rvi].limited	=	[0,0]
parinfo_all[wli].limited	=	[[1,1],[1,1],[0,0],[0,0],[0,0]]
parinfo_all[d_wli].limited	=	[[1,1],[0,0],[0,0]]
parinfo_all[ghi].limited	=	[[1,1],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
parinfo_all[d_ghi].limited	=	[[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
parinfo_all[ti].limited		=	[[1,1],[1,1],[1,1],[0,0],[0,0]]
parinfo_all[ki].limited		=	[[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1]]

;LIMITS
parinfo_all[rvi].tlimits	=	[-2d0,2d0]
parinfo_all[wli].tlimits	=	[[2.27d0,2.31d0],[4d-5,7d-5],[-1d-8,1d-8],[-1d-11,1d-11],[-1d-15,1d-15]]
parinfo_all[d_wli].tlimits	=	[[-1d5,1d-5],[-1d-6,1d-6],[-1d-7,1d-7]]
parinfo_all[ghi].tlimits	=	[[0.3d0,1.5d0],[-0.2d0,0.2d0],[-0.2d0,0.2d0],[-0.2d0,0.2d0],[-0.2d0,0.2d0],$
                                         [-0.2d0,0.2d0],[-0.2d0,0.2d0],[-0.2d0,0.2d0],[-0.2d0,0.2d0],[-0.2d0,0.2d0],[-0.2d0,0.2d0]]
parinfo_all[d_ghi].tlimits	=	[[-2d-1,2d-1],[-2d-1,2d-1],[-2d-1,2d-1],[-2d-1,2d-1],[-2d-1,2d-1],[-2d-1,2d-1],[-2d-1,2d-1],$
                                        [-2d-1,2d-1],[-2d-1,2d-1],[-2d-1,2d-1],[-2d-1,2d-1]]
parinfo_all[ti].tlimits		=	[[0d0,1d0],[0d0,1d0],[0d0,1d0],[0d0,1d0],[0d0,1d0]]
parinfo_all[ki].tlimits		=	[[0.8d0,1.2d0],[0.8d0,1.2d0],[0.8d0,1.2d0],[0.8d0,1.2d0],$
                                        [0.9d0,1.1d0],[0.9d0,1.1d0],[0.9d0,1.1d0],[0.9d0,1.1d0], $
					[0.9d0,1.1d0],[0.9d0,1.1d0],[0.9d0,1.1d0],[0.9d0,1.1d0]]
;MPSIDE
parinfo_all[rvi].mpside		=	2
parinfo_all[wli].mpside		=	[0, 0, 0, 0, 0]
parinfo_all[d_wli].mpside	=	[2, 0, 0]
parinfo_all[ghi].mpside		=	[2,0,0,0,0,0,0,0,0,0,0]



;parinfo_all_copy = parinfo_all

parinfo_all = parscale_rescale(parinfo_all)



parinfo_output = temporary(parinfo_all)

end
