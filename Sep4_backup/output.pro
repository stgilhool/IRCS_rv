; IDL Version 8.2.2 (linux x86_64 m64)
; Journal File for stgilhool@iroquois.physics.upenn.edu
; Working directory: /RAID/home/stgilhool/RV_projects/IRCS_rv/IDL_pros
; Date: Thu Jul  9 16:42:03 2015
 
af='../data/GJ273/18Jan2011/obs/IRCA00249843.0001.fits'
bf='../data/GJ273/18Jan2011/obs/IRCA00249845.0001.fits'
sf='../data/GJ273/18Jan2011/obs/AB1_sig.0001.fits'
a=mrdfits(af, 1)
; % MRDFITS: ERROR - Extension past EOF
a=mrdfits(af, 0)
;MRDFITS: Image array (1024,11)  Type=Real*4
help, a
b=mrdfits(bf, 0)
;MRDFITS: Image array (1024,11)  Type=Real*4
s=mrdfits(sf, 0)
;MRDFITS: Image array (1024,11)  Type=Real*4
ab=a-b
ab1d=total(ab, 2)
ab1d=total(ab, 2, /double)
s1d=sqrt(total(s^2, 2, /double))
plot, ab1d
bf='../data/GJ273/18Jan2011/obs/IRCA00249845A.0001.fits'
b=mrdfits(bf, 0)
;MRDFITS: Image array (1024,11)  Type=Real*4
ab=a-b
ab1d=total(ab, 2, /double)
plot, ab1d
plot, ab1d, /xs, yr=[2d5,5d5]
errplot, ab1d-s1d, ab1d+s1d
print, ab[5,*]
;      15504.2
;      16981.5
;      19455.3
;      20953.3
;      22100.8
;      22486.2
;      22357.7
;      20906.7
;      19790.1
;      17263.9
;      15894.5
print, s[5,*]
;      130.914
;      137.180
;      145.297
;      151.642
;      154.916
;      157.243
;      156.466
;      152.316
;      148.314
;      138.372
;      131.981
print, sqrt(ab[5,*])
;      124.516
;      130.313
;      139.482
;      144.753
;      148.663
;      149.954
;      149.525
;      144.592
;      140.677
;      131.392
;      126.073
print, ab[5,*]/s[5,*]
;      118.430
;      123.790
;      133.901
;      138.176
;      142.663
;      143.003
;      142.892
;      137.259
;      133.433
;      124.764
;      120.431
print, s[5,*]/ab[5,*]
;   0.00844377
;   0.00807821
;   0.00746820
;   0.00723713
;   0.00700952
;   0.00699286
;   0.00699831
;   0.00728550
;   0.00749437
;   0.00801512
;   0.00830354
test=mrdfits('../data/epoch/18Jan2011/GJ273_18Jan2011_AB1_1D.fits', 1)
; % MRDFITS: File access error
; % MRDFITS: OPENR: Null filename not allowed.
test=mrdfits('../data/epoch/18Jan2011/final_spectra/GJ273_18Jan2011_AB1_1D.fits', 1)
;MRDFITS: Binary table.  4 columns by  1 rows.
help, test
test=mrdfits('../data/epoch/18Jan2011/GJ273_18Jan2011_AB1.fits', 1)
; % MRDFITS: ERROR - Extension past EOF
test=mrdfits('../data/epoch/18Jan2011/GJ273_18Jan2011_AB1.fits', 0)
;MRDFITS: Image array (1024,1024)  Type=Real*4
help, test
test=mrdfits('../data/epoch/18Jan2011/GJ273_18Jan2011_AB1.fits', 1)
; % MRDFITS: ERROR - Extension past EOF
test=mrdfits('../data/epoch/18Jan2011/final_spectra/GJ273_18Jan2011_AB1_1D.fits', 1)
;MRDFITS: Binary table.  4 columns by  1 rows.
print, test.sigma[5]
;       1297.4070
print, test.spectrum[5]
;       869275.92
print, test.sigma[5]/test.spectrum[5]
;    0.0014925146
print, s1d[5]/ab1d[5]
;    0.0022686667
print, test.sigma[-6]/test.spectrum[-6]
;    0.0012637334
print, s1d[-6]/ab1d[-6]
;    0.0018792921
