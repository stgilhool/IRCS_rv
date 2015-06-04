; IDL Version 8.2.2 (linux x86_64 m64)
; Journal File for stgilhool@iroquois.physics.upenn.edu
; Working directory: /RAID/home/stgilhool/RV_projects/IRCS_rv/IDL_pros
; Date: Fri Feb 27 13:59:20 2015
pro rvreadin_journal
path='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/'
base='GJ273_18Jan2011_AB_0_127_'
ending='sign.fits'

file=path+base+ending

fits_info, file, /silent, n_ext=n_ext

rv=dblarr(n_ext)
;sample=mrdfits(file, 1)
structs=ptrarr(n_ext, /allocate_heap)

for i=0, n_ext-1 do begin
    *structs[i]=mrdfits(file, i+1)
    rv[i]=(*structs[i]).delta_rv
endfor

stop

end

