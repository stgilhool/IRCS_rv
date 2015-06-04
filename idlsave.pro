; IDL Version 8.2.2 (linux x86_64 m64)
; Journal File for stgilhool@iroquois.physics.upenn.edu
; Working directory: /RAID/home/stgilhool/RV_projects/IRCS_rv/IDL_pros
; Date: Mon Feb 23 06:41:13 2015
 
myjournal.pro
a=lindgen(5, 2)
print, a
;           0           1           2           3           4
;           5           6           7           8           9
b=reform(a, 10, 1)
print, b
;           0           1           2           3           4           5
;           6           7           8           9
c=reform(transpose(a), 10, 1)
print, c
;           0           5           1           6           2           7
;           3           8           4           9
print, transpose(a)
;           0           5
;           1           6
;           2           7
;           3           8
;           4           9
help, c
d=reform(transpose(a), 10)
print, d
;           0           5           1           6           2           7
;           3           8           4           9
help, d
print, n_elements(a)
;          10
