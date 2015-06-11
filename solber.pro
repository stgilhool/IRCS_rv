;  ///////////////
; SOLution BreedER
; ///////////////

;EDITED BY STEVE GILHOOL
;Corrected sorting mistake at end of generation
;Changed save function to save only unique solutions
;Made all parameter values double precision

; USE
; solution = SOLBER('myfitfun', ndim, funa = funa)
; 
; INPUTS
; ----------
;   fitfun - name of the function to compute fitness, here defined as
;          chi2, so minimum of fitness is searched for.
;          FUNCTION myfitfun, param, nparam, funa = funa
;          ...
;          return, fitness
;          end
;
;          nparam - number of parameter vectors to try
;          funa   - structure with additional functional arguments, 
;                   e.g. funa = {y:y, e:e}, where
;                   x  - vector of independent variables or vector of measurements
;                   y  - vector of measurements at each X
;                   e  - vector of errors for each Y
; ----------
;
;
;   
;   ndim   - dimensionality of the parameter vector
;
; INPUT KEYWORDS 
;   ...are explained further below.
;   Level 1 keywords are most important and are worth investigating
;   Level 2 keywords are not so important and can be left alone
;
; OUTPUTS
;   gfit_best - fitness of the solution found
;   status    - termination status
;   ngen_tot  - total number of generations
;   save_gen  - all accepted individuals
;   save_gfit - fitness of all accepted individuals

; Random permutation of index
FUNCTION _Rand_Perm, numberOfElements, numberInPermutation
  x = Lindgen(numberOfElements)
  y = x[Sort(Randomu(seed, numberOfElements))]
  RETURN, y(0:numberInPermutation-1)
END

; MAIN
FUNCTION SolBer, fitfun, ndim, funa = funa, lim = lim, npop = npop, crossrate = crossrate, mutrate = mutrate, ngen_max = ngen_max, term_fit = term_fit, term_flag = term_flag, nrep_frac = nrep_frac, ngen_inbred = ngen_inbred, dfit_inbred = dfit_inbred, mrate_min = mrate_min, mrate_max = mrate_max, difgfit_min = difgfit_min, delta = delta, plot_flag = plot_flag, print_flag = print_flag, gfit_best = gfit_best, status = status, ngen_tot = ngen_tot, save_gen = save_gen, save_gfit = save_gfit, new_save_gen = new_save_gen, new_save_gfit = new_save_gfit, new_save_igen = new_save_igen

; LEVEL 1 keywords
; ----------------
; FUNA = structure with additional functional arguments for fitness calculation

; Parameter limits
IF NOT(keyword_set(lim)) THEN BEGIN 
  lim = dblarr(2, ndim)
  lim[0, *]= replicate(0.0d0, ndim)
  lim[1, *] = replicate(1.0d0, ndim)
ENDIF 
; Population size
IF NOT(keyword_set(npop)) THEN npop = 300
; Crossover rate
IF NOT(keyword_set(crossrate)) THEN crossrate = 0.6
; Mutation rate
IF NOT(keyword_set(mutrate)) THEN mutrate = 0.03
; Maximum number of generations
IF NOT(keyword_set(ngen_max)) THEN ngen_max = 500
; Level of fitness required at termination
IF NOT(keyword_set(term_fit)) THEN term_fit = -1

; LEVEL 2 keywords
; ----------------
; Termination type
IF NOT(keyword_set(term_flag)) THEN term_flag = 0
; Maximum fraction of generation to replace
IF NOT(keyword_set(nrep_frac)) THEN nrep_frac = 0.4
; Number of generations with inbreeding allowed
IF NOT(keyword_set(ngen_inbred)) THEN ngen_inbred = 30
; Minimum relative fitness difference to define inbreeding
IF NOT(keyword_set(dfit_inbred)) THEN dfit_inbred = 0.1
; Minimum mutation allowed
IF NOT(keyword_set(mrate_min)) THEN mrate_min = mutrate
; Maximum mutation allowed
IF NOT(keyword_set(mrate_max)) THEN mrate_max = 0.5
; Minimum relative fitness to start increasing mutation rate
IF NOT(keyword_set(difgfit)) THEN difgfit_min = 0.1
; Rate of increase of mutation rate
IF NOT(keyword_set(delta)) THEN delta = 1.5
; Plotting flag, 0 = no plotting, 1 = plotting while solving, 2 = plot
; a report in the end
IF NOT(keyword_set(plot_flag)) THEN plot_flag = 0
; Printing flag
IF NOT(keyword_set(print_flag)) THEN print_flag = 0


; SETUP ends
; -------------------------------------

; OTHER settings
; Scaling
p0 = rebin(reform(lim[0, *]), ndim, npop)
dp = rebin(reform(lim[1, *]-lim[0, *]), ndim, npop)
; Boost
boost = 4
nrepmax = fix(nrep_frac*npop)
; Plotting
title1 = 'Fitness: best __ and 50th ....'
title2 = 'Fraction replaced __ and mutation rate ....'
chs = 0.9
col = 0
thick = 2
winx = 300
winy = 400

; Graphic window?
IF plot_flag EQ 1 THEN BEGIN  
  window, 0, xs = winx, ys = winy, retain = 2
  device, decomposed = 0
ENDIF 

; First generation
gen = randomu(seed, ndim, npop, /double)

; Evaluate generation fitness
dum = p0+dp*gen
gfit = call_function(fitfun, dum, npop, funa = funa)

; Save
IF keyword_set(save_gen) THEN save_gen = reform(p0+dp*gen, ndim, npop)
IF keyword_set(save_gfit) THEN save_gfit = gfit



; Initialize histories
hist_fit0 = [min(gfit)]
hist_fit50 = [median(gfit)]
hist_nrep = [nrep_frac]
hist_mut = [mutrate]
hist_dfit = [-1]

; For each generation
igen = 0D
REPEAT BEGIN 

; Sort fitness
  rank = sort(gfit)
; Rank generation
  gfit = gfit[rank]
  gen = gen[*, rank]

; Update mutation rate?
  difgfit = (gfit[0.5*npop-1]-gfit[0])/(gfit[0]+gfit[0.5*npop-1])
  IF difgfit LE difgfit_min THEN mutrate = min([mrate_max, mutrate*delta]) ELSE IF mutrate GT mrate_min THEN mutrate = max([mrate_min, mutrate/delta])

; Best in generation
  gbest = gen[*, 0]
  bestfit = gfit[0]

; Breeding Probability
  gbprob = (npop-findgen(npop))/npop
;  gbprob = gfit/max(gfit)
; probability folded 
  gbp2 = reform(rebin(gbprob, npop, boost), boost*npop)
; index folded
  ind = reform(rebin(indgen(npop), npop, boost), boost*npop)
; Select 2 parents
  nreptot = 0

; REPRODUCE
  uran = randomu(seed, boost*npop)
  sel = where(gbp2 GT uran, nsel)
  par1 = (ind[sel])[_rand_perm(nsel, npop)]
  par2 = (ind[sel])[_rand_perm(nsel, npop)]
; Crossover
  cross = randomu(seed, ndim, npop, /double)
  wcross = where(cross LT crossrate, complement = wncross, ncomplement = nwncross)
  IF nwncross GT 0 THEN cross[wncross] = 1
  ch1 = cross*gen[*, par1]+(1-cross)*gen[*, par2]
  ch2 = (1-cross)*gen[*, par1]+cross*gen[*, par2] 
  
; Mutate
  mutest1 = randomu(seed, ndim, npop)
  mut1 = randomu(seed, ndim, npop, /double)
  wmut1 = where(mutest1 LT mutrate, nmut1)
  mutest2 = randomu(seed, ndim, npop)
  mut2 = randomu(seed, ndim, npop, /double)
  wmut2 = where(mutest2 LT mutrate, nmut2)
  IF nmut1 GT 0 THEN ch1[wmut1] = mut1[wmut1]
  IF nmut2 GT 0 THEN ch2[wmut2] = mut2[wmut2]

; Evaluate children's fitness
  dum = p0+dp*ch1
  ch1fit = call_function(fitfun, dum, npop, funa = funa)
  dum = p0+dp*ch2
  ch2fit = call_function(fitfun, dum, npop, funa = funa)

; Keep the best children of each couple
  chfit = ch2fit
  ch = ch2
  wkeep1 = where(ch1fit LT ch2fit, n12)
  IF n12 GT 0 THEN BEGIN 
    chfit[wkeep1] = ch1fit[wkeep1]
    ch[*, wkeep1] = ch1[*, wkeep1]
  ENDIF 

; Rank children
  chrank = sort(chfit)

; Check if generation best needs replacing
    chbest = chfit[chrank[0]]
    IF chbest LT gfit[0] THEN BEGIN 
      bestfit = chbest
      gbest = ch[*, chrank[0]]
    ENDIF 

; Replace
    wrep1 = where(chfit LT gfit[par1], nrep1)
    IF nrep1 GT 0 THEN BEGIN 
      nrep1 = min([0.5*nrepmax, nrep1])
      gen[*, par1[wrep1[0:nrep1-1]]] = ch[*, wrep1[0:nrep1-1]]
      gfit[par1[wrep1[0:nrep1-1]]] = chfit[wrep1[0:nrep1-1]]
    ENDIF   
    wrep2 = where(chfit LT gfit[par2], nrep2)
    IF nrep2 GT 0 THEN BEGIN 
      nrep2 = min([0.5*nrepmax, nrep2])
      gen[*, par2[wrep2[0:nrep2-1]]] = ch[*, wrep2[0:nrep2-1]]
      gfit[par2[wrep2[0:nrep2-1]]] = chfit[wrep2[0:nrep2-1]]
    ENDIF 

    nreptot = nreptot+nrep1+nrep2

; Mutate all?
    IF igen GT ngen_inbred+1 THEN BEGIN 
      inbred = total(hist_dfit[igen-ngen_inbred-1:igen-1] GT dfit_inbred)
      IF inbred EQ 0 THEN BEGIN 
        mutest = randomu(seed, ndim, npop)
        mut = randomu(seed, ndim, npop, /double)
        wmut = where(mutest LT mutrate, nmut)
        IF nmut GT 0 THEN BEGIN 
          gen[wmut] = mut[wmut]
; Re-evaluate fitness
          dum = p0+dp*gen
          gfit = call_function(fitfun, dum, npop, funa = funa)
        ENDIF 
      enDIF 
    ENDIF 

    ; Re-sort the solutions, since children might outrank their
    ; parents' betters, and also because Mutate-all could change the
    ; ranking
    
    newrank=sort(gfit)
    gfit=gfit(newrank)
    gen=gen[*,newrank]

; Keep best
    gen[*, 0] = gbest
    gfit[0] = bestfit

; History
  hist_fit0 = [hist_fit0, bestfit]
  hist_fit50 = [hist_fit50, gfit[0.5*npop]]
  hist_mut = [hist_mut, mutrate]
  hist_dfit = [hist_dfit, difgfit]
  hist_nrep = [hist_nrep, float(nreptot)/npop]

; Save
  IF keyword_set(save_gen) THEN save_gen = [[save_gen], [reform(p0+dp*gen, ndim, npop)]]

  IF keyword_set(save_gfit) THEN save_gfit = [save_gfit, gfit]

; New save function (save unique solutions within 10 percent of best solution)
  IF (keyword_set(new_save_gen) OR keyword_set(new_save_gfit) $
      OR keyword_set(new_save_igen)) THEN BEGIN
      ;nsaves=n_elements(new_save_gen)
      gdiff=(gfit-bestfit)/bestfit
      potentials=where(gdiff LT 0.1, npotentials)
      IF igen EQ 0 THEN BEGIN
          new_save_gen=gbest
          new_save_gfit=bestfit
          new_save_igen=0
          FOR ipot = 1, npotentials-1 DO BEGIN
              ;Check to see if unique by checking first parameter
              ;If any parameter matches, then check entire solution
              sol_match=0
              par0_match_row=where((p0[0,0]+dp[0,0]*gen[0,potentials[ipot]]) EQ new_save_gen[0,*], n_match)
              IF n_match GT 0 THEN BEGIN
                  
                  FOREACH row_num, par0_match_row DO BEGIN
                      ;Check equality of entire solution
                      IF total((p0[*,0]+dp[*,0]*gen[*,potentials[ipot]]) NE new_save_gen[*,row_num]) EQ 0 THEN BEGIN
                          sol_match=1
                          break
                      ENDIF
                  ENDFOREACH
              ENDIF
              
              IF sol_match EQ 0 THEN BEGIN
                 new_save_gen=[[new_save_gen], [(p0[*,0]+dp[*,0]*gen[*,potentials[ipot]])]]    
                 new_save_gfit=[new_save_gfit, gfit[potentials[ipot]]]
                 new_save_igen=[new_save_igen, 0]
              ENDIF
          ENDFOR
          new_save_gen = dblarr(ndim)
          new_save_gfit= 0d0
          new_save_igen= 0
       ENDIF ELSE IF npotentials GT 0 THEN BEGIN
          keep_flag=replicate(1, npotentials)
          FOR ipot = 0, npotentials-1 DO BEGIN
              ;Check to see if unique by checking first parameter
              ;If any parameter matches, then check entire solution
              sol_match=0
              par0_match_row=where((p0[0,0]+dp[0,0]*gen[0,potentials[ipot]]) EQ new_save_gen[0,*], n_match)
              IF n_match GT 0 THEN BEGIN
                  
                  FOREACH row_num, par0_match_row DO BEGIN
                      ;Check equality of entire solution
                      IF total((p0[*,0]+dp[*,0]*gen[*,potentials[ipot]]) NE new_save_gen[*,row_num]) eq 0 THEN BEGIN
                          sol_match = 1
                          keep_flag[ipot]=0
                          break
                      ENDIF
                  ENDFOREACH
              ENDIF
          ENDFOR
          keep_ind=where(keep_flag EQ 1, nkeep)
          IF nkeep GT 0 THEN BEGIN
              keepers=potentials[keep_ind]
              this_save_gen=  rebin(p0[*,0], ndim, nkeep) + $
                (rebin(dp[*,0], ndim, nkeep) * gen[*,keepers])
              new_save_gen = [[new_save_gen], [this_save_gen]]
              
              
              new_save_gfit = [new_save_gfit, gfit[keepers]]
              
              new_save_igen = [new_save_igen, replicate(igen, nkeep)]
              
              ;SORT BY FITNESS
              save_rank=sort(new_save_gfit)
              new_save_gen=new_save_gen[*,save_rank]
              new_save_gfit=new_save_gfit[save_rank]
              new_save_igen=new_save_igen[save_rank]
              
          ENDIF
      ENDIF
  ENDIF
  



; PLOT
  IF plot_flag EQ 1 THEN BEGIN 
    !p.multi = [0, 1, 2]
    !x.margin = [6, 2]
    !y.margin = [3.5, 2]

; fitness
    yr = [0.8*min(hist_fit0), 1.2*max(hist_fit0)]
    plot, hist_fit0, yr = yr, ys = 1, /ylog, title = title1, background = 255, color = col, charsize = chs, yticks = 3
    oplot, hist_fit50, linestyle = 1, color = col, thick = thick
; number of genes replaced
    plot, hist_nrep, yr = [0, 1.], ys = 1, title = title2, color = col, charsize = chs
; mutation rate
    oplot, hist_mut, linestyle = 1, color = col, thick = thick
    !p.multi = 0
  ENDIF 

; PRINT
  IF print_flag EQ 1 THEN BEGIN 
    print, 'Generation', igen
    print, 'Best fitness', gfit[0]
    print, 'Best solution', reform(lim[0, *])+reform(lim[1, *]-lim[0, *])*gen[*, 0]
  ENDIF 

; Terminate
   CASE term_flag OF 
     0: BEGIN 
       ter1 = (igen EQ ngen_max) 
       ter2 = (bestfit LE term_fit)
       terminate = ter1 OR ter2
       status = ter1+2*ter2
     END 
     1: terminate = igen EQ ngen_max
     2: terminate = bestfit LE term_fit
   ENDCASE 
   igen = igen+1
  ENDREP UNTIL terminate

; Return results
gfit_best = gfit[0]
ngen_tot = fix(igen-1)

IF plot_flag EQ 2 THEN BEGIN  
  window, 0, xs = winx, ys = winy, retain = 2
  device, decomposed = 0
  !x.margin = [6, 2]
  !y.margin = [3.5, 2]
  !p.multi = [0, 1, 2]
; fitness
  yr = [0.8*min(hist_fit0), 1.2*max(hist_fit0)]
  plot, indgen(ngen_tot)+1, hist_fit0[1:ngen_tot], yr = yr, ys = 1, title = title1, background = 255, color = col, charsize = chs, /ylog, yticks = 3
  oplot, indgen(ngen_tot)+1, hist_fit50[1:ngen_tot], linestyle = 1, color = col, thick = thick
; number of genes replaced
  plot, indgen(ngen_tot)+1, hist_nrep[1:ngen_tot], yr = [0, 1.], ys = 1, title = title2, color = col, charsize = chs
; mutation rate
  oplot, indgen(ngen_tot)+1, hist_mut[1:ngen_tot], linestyle = 1, color = col, thick = thick
  !p.multi = 0
ENDIF 

return, reform(lim[0, *])+reform(lim[1, *]-lim[0, *])*gen[*, 0]
END 
