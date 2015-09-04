pro ircsrv_calibrate_caller, parinfo_all, _extra=ex

!except=2

ircsrv_calibrate_jul15, parinfo_all, _extra=ex
;object=object, epoch=epoch, trace=trace, first_pix=first_pix, npix_select=npix_select, fmode=fmode, fit_measure=fit_measure, visualize=visualize, outfile=outfile


end



pro ircsrv_calibrate_wrapper

!except=2

;TUNE THESE
outfile='jul30ingabestoutmpfit.fits'
infile='5_2_0_6_best.fits'

object='GJ273'
epoch='18Jan2011'
trace='AB1'

first_pix=0L
npix_select=512L
fmode='mpfit'
fit_measure='chi_vec'
visualize=0


;THESE TOO
n_wl	= 5L
n_gh	= 2L
n_d_gh	= 0L
n_t	= 1L
n_k	= 6L

ex={object:object, $
    epoch:epoch, $
    trace:trace, $
    first_pix:first_pix, $
    npix_select:npix_select, $
    fmode:fmode, $
    fit_measure:fit_measure, $
    visualize:visualize, $
    outfile:outfile $
    }

iter=0


repeat begin

    if iter eq 0 then begin
        parinfo_all = parinfo_readin(infile, /rawform)
        i = where(finite(parinfo_all.tbestfit), icount)
        if icount gt 0 then begin
            parinfo_all[i].tvalue = parinfo_all[i].tbestfit
            parinfo_all = parscale_rescale(parinfo_all)
        endif
    endif else parinfo_all = parinfo_readin(infile)
    
    parinfo_all_freepar, parinfo_all, n_wl=n_wl, n_gh=n_gh, n_d_gh=n_d_gh, n_t = n_t, n_k=n_k
    
parname1=['wl0','wl1','wl2','wl3','wl4','t_nh3', 'k0','k1','k2','k3','k4','k5']
npar1 = n_elements(parname1)
auto1= replicate(1d-3, npar1)

ranfactor1 = replicate(5d-2, npar1)
ranfactor1[[0,1]] = 5d-5

;set relative limits
    parinfo_all = parinfo_changelim(parinfo_all, parname=parname1, auto=auto1)
    
    parinfo_all = parinfo_changelim(parinfo_all,  parname=['gh0'], auto=[1d-1])

    parinfo_all = parinfo_changelim(parinfo_all, [-0.15d0,0.15d0], parname=['gh1'])
   
;randomize starting guesses
;    parinfo_all = parinfo_randomval(parinfo_all, parname=parname1, ranfactor=ranfactor1)

;    parinfo_all = parinfo_randomval(parinfo_all, parname=['gh0','gh1'], ranfactor=[3d-1, 4d0])
    
;set step size to be 5% of ulim-llim
    parinfo_all.tstep = (parinfo_all.tlimits[1] - parinfo_all.tlimits[0]) * 5d-2
    parinfo_all = parscale_rescale(parinfo_all)

;fit
    ircsrv_calibrate_caller, parinfo_all, _extra=ex
    
    infile = '5_2_0_6_best.fits'

    iter++

    model_path = '/home/stgilhool/RV_projects/IRCS_rv/data/epoch/18Jan2011/calib_results/'
    model_infof = model_path+'5_2_0_6_best.fits'
    last_infof = model_path + outfile
    model_info = mrdfits(model_infof, 1)
    last_info = mrdfits(last_infof, 1)
    chi2 = model_info.fitness
    
    term = (chi2 lt 4000d0) or (iter ge 1)

    print, 'iter: ', iter
    print, 'best: ', chi2
    print, 'latest: ', last_info.fitness
    
    wait, 60
    
endrep until term eq 1


end
