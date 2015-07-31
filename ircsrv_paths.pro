function ircsrv_paths, epoch=epoch, object=object, trace=trace

if n_elements(epoch) eq 0 then epoch='18Jan2011'
if n_elements(object) eq 0 then object='GJ273'
if n_elements(trace) eq 0 then trace='AB'

root	=	'/home/stgilhool/RV_projects/IRCS_rv/'

data	=	root+'data/'

sup	=	data+'supplemental/'
objroot	=	data+object+'/'
flatroot=	data+'cal_files/'
eproot	=	data+'epoch/'

obj	=	objroot+epoch+'/'
flat	=	flatroot+epoch+'/'
ep	=	eproot+epoch+'/'

obs	=	obj+'obs/'
onin	=	flat+'calONcellIN/'
offin	=	flat+'calOFFcellIN/'
offout	=	flat+'calOFFcellOUT/'
onout	=	flat+'calONcellOUT/'

calib	=	ep+'calib_results/'
temp	=	ep+'temp_results/'
rv	=	ep+'rv_results/'
spec	=	ep+'final_spectra/'


pathstr={root:root, $
         data:data, $
         objroot:objroot, $
         flatroot:flatroot, $
         eproot:eproot, $
         obj:obj, $
         flat:flat, $
         ep:ep, $
         obs:obs, $
         onin:onin, $
         offin:offin, $
         offout:offout, $
         onout:onout, $
         calib:calib, $
         sup:sup, $
         temp:temp, $
         rv:rv, $
         spec:spec $
        }

return, pathstr

end
