pro ircsrv_templategenerate, model_tag

if n_elements(model_tag) eq 0 then model_tag='default'

 phoenix_path='/home/stgilhool/RV_projects/IRCS_rv/data/supplemental/stellar_template_models/'

wl_file='WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'

temp_file1='lte03300-5.00-1.5.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
temp_file2='lte03300-5.00-1.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
temp_file3='lte03300-5.00-0.5.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
temp_file4='lte03300-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
temp_file5='lte03300-5.00-2.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
temp_file6='lte03300-5.00-3.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
temp_file7='lte03300-5.00+0.5.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
temp_file8='lte03300-5.00+1.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'

tf1=phoenix_path+temp_file1
tf2=phoenix_path+temp_file2
tf3=phoenix_path+temp_file3
tf4=phoenix_path+temp_file4
tf5=phoenix_path+temp_file5
tf6=phoenix_path+temp_file6
tf7=phoenix_path+temp_file7
tf8=phoenix_path+temp_file8

wlf=phoenix_path+wl_file


wl=mrdfits(wlf, 0)

t1=mrdfits(tf1, 0)
t2=mrdfits(tf2, 0)
t3=mrdfits(tf3, 0)
t4=mrdfits(tf4, 0)
t5=mrdfits(tf5, 0)
t6=mrdfits(tf6, 0)
t7=mrdfits(tf7, 0)
t8=mrdfits(tf8, 0)

wrange=where(wl le 23500d0 and wl ge 22900d0)

stop

end


