function paths, epoch=epoch, object=object, trace=trace

if n_elements(epoch) eq 0 then epoch='18Jan2011'
if n_elements(object) eq 0 then object='GJ273'
if n_elements(trace) eq 0 then trace='AB'

root='/home/stgilhool/RV_projects/IRCS_rv/'

dataroot=root+'data/'

objectroot=dataroot+object+'/'
flatroot=dataroot+'cal_files/'
epochroot=dataroot+'epoch/'

objectfold=objectroot+epoch+'/'
flatfold=flatroot+epoch+'/'
epochfold=epochroot+epoch+'/'

obsfold=objectfold+'obs/'
flatONINfold=flatfold+'calONcellIN/'
flatOFFINfold=flatfold+'calOFFcellIN/'
flatOFFOUTfold=flatfold+'calOFFcellOUT/'
flatONOUTfold=flatfold+'calONcellOUT/'

calfold=epochfold+'calib_results/'
tempfold=epochfold+'temp_results/'
rvfold=epochfold+'rv_results/'
specfold=epochfold+'final_spectra/'


pathstr={root:root, $
	dataroot:dataroot, $
	objectroot:objectroot, $
	flatroot:flatroot, $
	epochroot:epochroot, $
	objectfold:objectfold, $
	flatfold:flatfold, $
	epochfold:epochfold, $
	obsfold:obsfold, $
	ONINfold:flatONINfold, $
	OFFINfold:flatOFFINfold, $
	OFFOUTfold:flatOFFOUTfold, $
	ONOUTfold:flatONOUTfold, $
	calfold:calfold, $
	tempfold:tempfold, $
	rvfold:rvfold, $
	specfold:specfold $
	}

return, pathstr

end