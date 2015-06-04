; IDL Version 8.2.2 (linux x86_64 m64)
; Journal File for stgilhool@iroquois.physics.upenn.edu
; Working directory: /RAID/home/stgilhool/RV_projects/IRCS_rv/IDL_pros
; Date: Wed Mar  4 15:18:25 2015
 
rfile='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_127_smooth30_r2.fits'
fits_info, rfile, /silent, n_ext=n_ext
tstr=ptrarr(n_ext, /allocate_heap)
for ext=0,n_ext-1 do begin & $
*tstr[ext]=mrdfits(rfile, ext+1)
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
help, tstr
help, *tstr[12]
rfile='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/GJ273_18Jan2011_AB_0_127_ttest_r2.fits'
fits_info, rfile, /silent, n_ext=n_ext
astr=ptrarr(n_ext, /allocate_heap)
for ext=0,n_ext-1 do begin & $
*astr[ext]=mrdfits(rfile, ext+1)
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
tdwl=(*tstr[*]).guess[3,4]
; % Expression must be a scalar in this context: <POINTER   Array[14]>.
tdwl=*tstr[*].guess[3,4]
; % Expression must be a structure in this context: TSTR.
tdwl=*(tstr[*]).guess[3,4]
; % Expression must be a structure in this context: <POINTER   Array[14]>.
help, *tstr
; % Expression must be a scalar in this context: TSTR.
help, *tstr[0]
help, *tstr[0].obj
; % Expression must be a structure in this context: TSTR.
help, (*tstr[0]).obj
print, (*tstr[*]).obj
; % Expression must be a scalar in this context: <POINTER   Array[14]>.
print, ((*tstr[*]).obj)
; % Expression must be a scalar in this context: <POINTER   Array[14]>.
for i=0, n_ext-1 do begin & $
endfor
tdwl=dblarr(n_ext)
adwl=dblarr(n_ext)
for i=0, n_ext-1 do begin & $
tdwl[i]=(*tstr[i]).guess[3,4] & $
adwl[i]=(*astr[i]).guess[3,4] & $
endfor
; % Illegal subscript range: <No name>.
for i=0, n_ext-1 do begin & $
adwl[i]=(*astr[i]).(guess[3,4]) & $
endfor
; % Variable is undefined: GUESS.
print, (*astr[0]).guess
;      -308.96101     0.050667181      0.58432496   4.5803957e-06  -6.5111987e-08      0.65557457     -0.11203322   0.00049931185
;  -0.00081985216      0.90375530       1.2046356       1.1691917       1.1624649       1.1686384       1.1684193       1.1680520
;       1.1634166       1.1638898       1.1710953       1.1664465       1.1664963       1.1632989
print, ((*astr[0]).guess)[3,4]
; % Attempt to subscript <DOUBLE    Array[22]> with <INT      (       4)> is out of range.
print, ((*astr[0]).guess[3,4])
; % Illegal subscript range: <No name>.
print, (((*astr[0]).guess)[3,4])
; % Attempt to subscript <DOUBLE    Array[22]> with <INT      (       4)> is out of range.
print, (*astr[0]).guess
;      -308.96101     0.050667181      0.58432496   4.5803957e-06  -6.5111987e-08      0.65557457     -0.11203322   0.00049931185
;  -0.00081985216      0.90375530       1.2046356       1.1691917       1.1624649       1.1686384       1.1684193       1.1680520
;       1.1634166       1.1638898       1.1710953       1.1664465       1.1664963       1.1632989
help, (*astr[0]).guess
for i=0, n_ext-1 do begin & $
aguess=(*astr.[i]).guess & $
; % Syntax error.
aguess=(*astr[i]).guess & $
adwl[i]=aguess[3,4] & $
endfor
print, adwl
;       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
;       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
for i=0, n_ext-1 do begin & $
aguess=(*astr[i]).guess & $
adwl[i]=aguess[3,4] & $
endfor
; % Attempt to subscript AGUESS with <INT      (       4)> is out of range.
print, aguess
;      -308.96101     0.050667181      0.58432496   4.5803957e-06  -6.5111987e-08      0.65557457     -0.11203322   0.00049931185
;  -0.00081985216      0.90375530       1.2046356       1.1691917       1.1624649       1.1686384       1.1684193       1.1680520
;       1.1634166       1.1638898       1.1710953       1.1664465       1.1664963       1.1632989
print, ((*astr[0]).guess[3:4])
;   4.5803957e-06  -6.5111987e-08
for i=0, n_ext-1 do begin & $
adwl[i]=(*astr[i]).guess[3:4] & $
tdwl[i]=(*tstr[i]).guess[3:4] & $
endfor
; % Out of range subscript encountered: ADWL.
help, adwl
help, tdwl
adwl=dblarr(2,n_ext)
tdwl=dblarr(2,n_ext)
for i=0, n_ext-1 do begin & $
adwl[*, i]=(*astr[i]).guess[3:4] & $
tdwl[*, i]=(*tstr[i]).guess[3:4] & $
endfor
print, adwl[0,*]-tdwl[0,*]
;   2.4164756e-07
;  -4.0355939e-08
;   3.2949784e-07
;   2.5155602e-07
;   1.6109404e-07
;  -1.2370664e-06
;   3.0508005e-07
;  -8.0367961e-07
;   2.2293895e-08
;  -2.2146470e-07
;   3.4663323e-07
;   7.8874558e-07
;   1.0775984e-06
;   9.6413531e-07
print, adwl[1,*]-tdwl[1,*]
;  -2.0972718e-08
;  -2.3633493e-08
;  -2.1776113e-08
;  -2.2733626e-08
;  -1.9843293e-08
;   6.9544562e-11
;  -2.1524669e-08
;   2.3134764e-09
;  -1.9643477e-08
;  -2.0748534e-08
;  -2.4443888e-08
;  -3.1454524e-08
;  -3.3190138e-08
;  -2.3337288e-08
.run testing
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Type of end does not match statement (ENDFOREACH expected).
; % Syntax error.
; % Type of end does not match statement (END expected).
; % 11 Compilation error(s) in module TESTING.
; % Syntax error.
; % Type of end does not match statement (END expected).
; % 2 Compilation error(s) in module $MAIN$.
; % Syntax error.
; % 1 Compilation error(s) in module $MAIN$.
.run testing
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Type of end does not match statement (ENDFOREACH expected).
; % Syntax error.
; % Type of end does not match statement (END expected).
; % 11 Compilation error(s) in module TESTING.
; % Syntax error.
; % Type of end does not match statement (END expected).
; % 2 Compilation error(s) in module $MAIN$.
; % Syntax error.
; % 1 Compilation error(s) in module $MAIN$.
.run testing
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % 7 Compilation error(s) in module TESTING.
.run testing
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % 7 Compilation error(s) in module TESTING.
.run testing
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % Syntax error.
; % 6 Compilation error(s) in module TESTING.
.run testing
testing
;0 Result of first template, 2 iterations RV fit
;1 Result of fitting for a better stellar template
;2 Result of 1 iteration with new stellar template (RV fixed)
;3 Result of 2nd iteration with new stellar template (all parameters free)
;4 Result of fitting for a better telluric template
;5 Result of 1 iteration with new stellar and telluric (RV fixed)
;6 Result of 2nd iteration with new stellar and telluric (all parameters free)
; % MRDFITS: File access error
; % MRDFITS: OPENR: Null filename not allowed.
; % Expression must be a structure in this context: STR.
help, list
print, list
;/home/stgilhool/RV_projects/IRCS_rv/rvshift1_results/GJ273_18Jan2011_AB_0_127_run2sign.fits
;/home/stgilhool/RV_projects/IRCS_rv/data/smooth_penalty_test/test16/penaltytest.fits
;/home/stgilhool/RV_projects/IRCS_rv/rvshift1_results/GJ273_18Jan2011_AB_0_127_smooth30.fits
;/home/stgilhool/RV_projects/IRCS_rv/rvshift1_results/GJ273_18Jan2011_AB_0_127_smooth30_r2.fits
;/home/stgilhool/RV_projects/IRCS_rv/data/telluric_test/test5/tellurictest.fits
;/home/stgilhool/RV_projects/IRCS_rv/rvshift1_results/GJ273_18Jan2011_AB_0_127_ttest2_r1.fits
;/home/stgilhool/RV_projects/IRCS_rv/rvshift1_results/GJ273_18Jan2011_AB_0_127_ttest2_r2.fits
print, entry
;/home/stgilhool/RV_projects/IRCS_rv/rvshift1_results/GJ273_18Jan2011_AB_0_127_run2sign.fits
help, entry
.run testing
testing
;0 Result of first template, 2 iterations RV fit
;1 Result of fitting for a better stellar template
;2 Result of 1 iteration with new stellar template (RV fixed)
;3 Result of 2nd iteration with new stellar template (all parameters free)
;4 Result of fitting for a better telluric template
;5 Result of 1 iteration with new stellar and telluric (RV fixed)
;6 Result of 2nd iteration with new stellar and telluric (all parameters free)
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  20 columns by  1 rows.
; % Expression must be a scalar or 1 element array in this context: <BYTE      Array[22]>.
help, str
help, strspec
help, strspec.params
help, str.result
print, strspec.params eq str.result
;   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
.run testing
testing
;0 Result of first template, 2 iterations RV fit
;1 Result of fitting for a better stellar template
;2 Result of 1 iteration with new stellar template (RV fixed)
;3 Result of 2nd iteration with new stellar template (all parameters free)
;4 Result of fitting for a better telluric template
;5 Result of 1 iteration with new stellar and telluric (RV fixed)
;6 Result of 2nd iteration with new stellar and telluric (all parameters free)
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  20 columns by  1 rows.
;   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
;       2231.5962
;       1030.1764
.run testing
; % Type of end does not match statement (ENDFOREACH expected).
; % Syntax error.
; % Type of end does not match statement (END expected).
; % 3 Compilation error(s) in module TESTING.
; % Syntax error.
; % Type of end does not match statement (END expected).
; % 2 Compilation error(s) in module $MAIN$.
; % Syntax error.
; % 1 Compilation error(s) in module $MAIN$.
.run testing
testing
;0 Result of first template, 2 iterations RV fit
;1 Result of fitting for a better stellar template
;2 Result of 1 iteration with new stellar template (RV fixed)
;3 Result of 2nd iteration with new stellar template (all parameters free)
;4 Result of fitting for a better telluric template
;5 Result of 1 iteration with new stellar and telluric (RV fixed)
;6 Result of 2nd iteration with new stellar and telluric (all parameters free)
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  20 columns by  1 rows.
;   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
;MRDFITS: Binary table.  16 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  20 columns by  1 rows.
;   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  20 columns by  1 rows.
;   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
;MRDFITS: Binary table.  24 columns by  1 rows.
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  20 columns by  1 rows.
;   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
;MRDFITS: Binary table.  31 columns by  1 rows.
;MRDFITS: Binary table.  20 columns by  1 rows.
;   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
; % Program caused arithmetic error: Floating divide by 0
print, chi2
;       1030.1764        Infinity       3955.3442       3678.6849        Infinity       5663.4401       5607.8464
