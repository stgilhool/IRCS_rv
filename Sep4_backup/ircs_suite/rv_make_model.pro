function rv_make_model, ideal_model, lsf, norm_vector
;+
; NAME: rv_make_model
;
; PURPOSE: Make a model observation from nh3, telluric, and stellar models
; convolved with an LSF and corrected for normalization errors
;
;
; CATEGORY:
;
; CALLING SEQUENCE: model=rv_make_model(ideal_model, lsf, norm_vector)
;
; INPUTS: 
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMON BLOCKS:
;
;
; SIDE EFFECTS:
;
;
; RESTRICTIONS:
;
;
; PROCEDURE:
;
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;
;-

;;;Add section to check inputs

;;;Convolve
model_over=ircsrv_convolve(ideal_model, lsf) 
