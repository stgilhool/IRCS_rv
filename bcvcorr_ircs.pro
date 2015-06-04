pro bcvcorr_ircs, head, params

;bcvcorr modified for IRCS observations
;RA and DEC in degrees
rasex=sxpar(head, 'RA2000')
decsex=sxpar(head, 'DEC2000')

ra=double(strsplit(rasex, ':', /extract))
;ra[0]=ra[0]*(360d0/24d0) ;this converts HH to DD
;stop

dec0=ten(decsex)

dateobs=SXPAR(head,'DATE-OBS')
date=strsplit(dateobs, '-', /extract)

ut=sxpar(head, 'UT')
utday2=double(strsplit(ut, ':', /extract))

openw, file_unit, 'bcv_temp.input', /GET_LUN

printf,file_unit,'2'


printf, file_unit, date[1], date[2], date[0], format='(I2, 2x, I2, 2x, I4)'


printf, file_unit, utday2[0]+utday2[1]/60d0+utday2[2]/3600d0+sxpar(head, 'EXPTIME')*.5/3600d0, format='(D)'


printf, file_unit, ra[0], ra[1], ra[2], format='(I2, 2x, I2, 2x, I2)'

dec=dec0

printf, file_unit, dec0, format='(D)'

printf,file_unit,STRTRIM(STRING(2000),2)

;these are the Subaru coords

longitude=STRTRIM(STRING(155.48056),2)
latitude=STRTRIM(STRING(19.828611),2)
elevation=STRTRIM(STRING(4139),2)


printf,file_unit,latitude
printf,file_unit,longitude
printf,file_unit,elevation

printf,file_unit,STRTRIM(STRING(0),2)

close,file_unit
free_lun,file_unit

spawn, './runbcv < bcv_temp.input', bcv_temp


v_bits=STRSPLIT(bcv_temp(12),':',/extract)
v=FLOAT(STRTRIM(v_bits(1),2))

params=DBLARR(3)


params(0)=v
;params[1]=helio_jd(date_conv(ut, 'JULIAN')-2400000d, ra0, dec0)+2400000d
;params[2]=date_conv(ut, 'JULIAN')+sxpar(head, 'EXPTIME')*.5/(3600d0*24d0)

end
