pro show_results

file_base='/home/stgilhool/RV_projects/IRCS_rv/data/rvshift1_results/'

large='1000'
med2='100'
med1='50'
small='10'

np=large
ng=small
file_name=file_base+'GJ273_18Jan2011_AB_310_398_Jun12_ga2_'+np+'_'+ng+'_ga.fits'
ls=mrdfits(file_name,1)

np=small
ng=large
file_name=file_base+'GJ273_18Jan2011_AB_310_398_Jun12_ga2_'+np+'_'+ng+'_ga.fits'
sl=mrdfits(file_name,1)

np=med2
ng=small
file_name=file_base+'GJ273_18Jan2011_AB_310_398_Jun12_ga2_'+np+'_'+ng+'_ga.fits'
m2s=mrdfits(file_name,1)

np=small
ng=small
file_name=file_base+'GJ273_18Jan2011_AB_310_398_Jun12_ga2_'+np+'_'+ng+'_ga.fits'
ss=mrdfits(file_name,1)

np=med1
ng=med2
file_name=file_base+'GJ273_18Jan2011_AB_310_398_Jun12_ga2_'+np+'_'+ng+'_ga.fits'
m1m2=mrdfits(file_name,1)

np='500'
ng='150'
file_name=file_base+'GJ273_18Jan2011_AB_310_398_Jun15_iter_'+np+'_'+ng+'_mpfit.fits'
mpbig=mrdfits(file_name, 1)
file_name=file_base+'GJ273_18Jan2011_AB_310_398_Jun14_ga2_'+np+'_'+ng+'_mpfit.fits'
mpbig0=mrdfits(file_name, 1)

ng='20'
file_name=file_base+'GJ273_18Jan2011_AB_310_398_Jun15_iter_'+np+'_'+ng+'_mpfit.fits'
mpsmall=mrdfits(file_name, 1)
file_name=file_base+'GJ273_18Jan2011_AB_310_398_Jun14_ga2_'+np+'_'+ng+'_mpfit.fits'
mpsmall0=mrdfits(file_name, 1)

file_name=file_base+'i50020.fits'
gasmall=mrdfits(file_name, 1)

file_name=file_base+'i500150.fits'
gabig=mrdfits(file_name, 1)

print, "Many gens (150)"
print, 'MP 0: ', mpbig0.chi2
print, 'GA : ', gabig.chi2
print, 'MP 1: ', mpbig.chi2
print, ''
print, '' 
print, "Few gens (20)"
print, 'MP 0: ', mpsmall0.chi2
print, 'GA : ', gasmall.chi2
print, 'MP 1: ', mpsmall.chi2

stop

end

