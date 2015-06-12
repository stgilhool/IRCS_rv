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


help, ls

stop

end

