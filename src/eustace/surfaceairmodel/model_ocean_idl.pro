; IDL script ported from inside of loop in original model
; This is to be run from a command line like:
;
; idl -e ".run model_ocean_idl.pro" -args {jsonfile}
;
; where {jsonfile} is path to text file in JSON format containing:
;
; { "model": { "filename_parameters": "filename",
;              "filename_uncertainty": "filename", 
;              "filename_stdev": "filename" },
;   "data" : [ { "date": { "year": yyyy, "month": mm, "day": dd },
;                "filename_input": "filename",
;                "filename_output": "filename" },
;              ... ,
;            ]
; }

@model_ocean_idl_setup_satellite.pro

;---------- read config file which tells us what to process

; Get command line arguments
args = command_line_args()

; Check single argument (jsonfile)
if n_elements(args) ne 1 then message, 'ERROR: expected 1 argument (jsonfile)'

; Read text file to input string
inputstring = ''
openr, file_handle, string(args[0]), /GET_LUN
while not eof(file_handle) do begin
   line = ''
   readf, file_handle, line
   inputstring = inputstring + line
endwhile
free_lun, file_handle

; Parse JSON
json = json_parse(inputstring)
json_model = json['model']
json_data_list = json['data']

; Extract model file details from JSON - example contents should be like:
; filename_model_clim_parameters='/gws/nopw/j04/eustace/data/internal/surfaceair_model_parameters/ocean/hires_AST_clim_parameters_v.5.0.dat'
; filename_model_clim_parameters_uncertainty='/gws/nopw/j04/eustace/data/internal/surfaceair_model_parameters/ocean/hires_AST_clim_parameters_uncertainty_v.5.0.dat'
; filename_model_STDEV_clim_parameters='/gws/nopw/j04/eustace/data/internal/surfaceair_model_parameters/ocean/hires_AST_STDEV_clim_parameters_v.5.0.dat'
filename_model_clim_parameters = json_model['filename_parameters']
filename_model_clim_parameters_uncertainty = json_model['filename_uncertainty']
filename_model_STDEV_clim_parameters = json_model['filename_stdev']

;---------- Preparation for daily processing from top of original script

; Open IDL variable information
restore, filename_model_clim_parameters
restore, filename_model_clim_parameters_uncertainty
restore, filename_model_STDEV_clim_parameters

; Constants from top of original script
tstandard_name =  'surface_temperature'
tunits         =  'K'
tfill          = -32768L
tadd_offset    =  0.0
tscale_factor  =  0.001
tvalid_min     = -100000L
tvalid_max     =  100000L
uadd_offset    =  0.
uscale_factor  =  0.001
uvalid_min     =  0L
uvalid_max     =  10000L
coordinates    =  'lat lon' 
tcomment       =  'gridded marine surface air temperatures estimated from satellite data'
geo_fill       = -32768.
npts = 365
neofs = 4

; Fourier initialisation from top of original script
get_fourier_dailies, neofs, C, C_daily
mn = get_hires_offset(C, hi_res_params,       neofs, 0.0, -179.5, 0)
sd = get_hires_offset(C, hi_res_stdev_params, neofs, 0.0, -179.5, 0)

;----------  Loop over input files

foreach json_data, json_data_list do begin

    ; Extract daily data information from JSON - example contents should be like:
    ; yyy=2007
    ; mmm=1
    ; ddd=1
    ; filename_input='/gws/nopw/j04/eustace/data/incoming/UoR_sst/l3c/AATSR/2007/01/200701011200-ESACCI-L3C_GHRSST-AATSR-EXP1.2-v02.0-fv1.0.nc'
    ; filename_output='ocean.nc'
    json_date = json_data['date']
    yyy = json_date['year']
    mmm = json_date['month']
    ddd = json_date['day']
    filename_input = json_data['filename_input']
    filename_output = json_data['filename_output']

;---------- Inner part of daily loop from original script

;open and read in satellite data for today
            Cdfid = NCDF_OPEN(filename_input)
            NCDF_VARGET, Cdfid, 'sea_surface_temperature_depth', sst
            NCDF_VARGET, Cdfid, 'uncorrelated_uncertainty', random_unc
            NCDF_VARGET, Cdfid, 'synoptically_correlated_uncertainty', local_systematic
            NCDF_VARGET, Cdfid, 'large_scale_correlated_uncertainty', systematic
            NCDF_VARGET, Cdfid, 'lat', latitudes
            NCDF_VARGET, Cdfid, 'lon', longitudes
            NCDF_VARGET, Cdfid, 'time', time
            NCDF_CLOSE, Cdfid

            sst = sst * 0.01
            random_unc = random_unc * uscale_factor
            local_systematic = local_systematic * uscale_factor
            systematic = systematic * uscale_factor
            
            day_in_year = ymd2dn(yyy, mmm, ddd)
            if leap(yyy) and day_in_year gt 59 then begin
               day_in_year = day_in_year - 1
            endif
            
            mat = sst
            mat[*] = -32768
            mat_local_systematic = mat
            mat_systematic = mat
            mat_systematic[*,*] = 0.1

            mat_param0_systematic = mat
            mat_param1_systematic = mat
            mat_param2_systematic = mat
            mat_param3_systematic = mat
            mat_param4_systematic = mat

            param0_systematic = mat
            param1_systematic = mat
            param2_systematic = mat
            param3_systematic = mat
            param4_systematic = mat
            
;loop over all grid boxes
            min_sst = min(sst)

            for xx = 0,n_elements(longitudes)-1 do begin
            for yy = 0,n_elements(latitudes)-1 do begin

               if sst[xx,yy] ne min_sst then begin

                  lat = latitudes[yy]
                  lon = longitudes[xx]
                  
                  mn_conv = get_hires_offset(C, hi_res_params,       neofs, lat, lon, day_in_year-1)
                  sd_conv = get_hires_offset(C, hi_res_stdev_params, neofs, lat, lon, day_in_year-1)

;Uncertainties in the fourier coefficients
                  C_flat = fltarr(365,5)+1.0
                  param_0_unc = get_single_hires_offset(C_flat, hi_res_params_unc,0, neofs, lat, lon, day_in_year-1)
                  param_1_unc = get_single_hires_offset(C_flat, hi_res_params_unc,1, neofs, lat, lon, day_in_year-1)
                  param_2_unc = get_single_hires_offset(C_flat, hi_res_params_unc,2, neofs, lat, lon, day_in_year-1)
                  param_3_unc = get_single_hires_offset(C_flat, hi_res_params_unc,3, neofs, lat, lon, day_in_year-1)
                  param_4_unc = get_single_hires_offset(C_flat, hi_res_params_unc,4, neofs, lat, lon, day_in_year-1)

;uncertainty in temperatures associated with each coefficeint
                  param_0_temp_unc = get_single_hires_offset(C, hi_res_params_unc,0, neofs, lat, lon, day_in_year-1)
                  param_1_temp_unc = get_single_hires_offset(C, hi_res_params_unc,1, neofs, lat, lon, day_in_year-1)
                  param_2_temp_unc = get_single_hires_offset(C, hi_res_params_unc,2, neofs, lat, lon, day_in_year-1)
                  param_3_temp_unc = get_single_hires_offset(C, hi_res_params_unc,3, neofs, lat, lon, day_in_year-1)
                  param_4_temp_unc = get_single_hires_offset(C, hi_res_params_unc,4, neofs, lat, lon, day_in_year-1)

;The fourier series themselves
                  param_0 = C[day_in_year-1,0] 
                  param_1 = C[day_in_year-1,1] 
                  param_2 = C[day_in_year-1,2] 
                  param_3 = C[day_in_year-1,3] 
                  param_4 = C[day_in_year-1,4] 
               
                  if abs(mn_conv) lt 10 and abs(sd_conv) le 10 then begin

                     mat[xx,yy] = sst[xx,yy] + mn_conv
                     if sd_conv lt 0.3 then sd_conv = 0.3
                     mat_local_systematic[xx,yy] = sd_conv

                     param0_systematic[xx,yy] = param_0_unc 
                     param1_systematic[xx,yy] = param_1_unc 
                     param2_systematic[xx,yy] = param_2_unc 
                     param3_systematic[xx,yy] = param_3_unc 
                     param4_systematic[xx,yy] = param_4_unc 

                     mat_param0_systematic[xx,yy] = abs(param_0_temp_unc) 
                     mat_param1_systematic[xx,yy] = abs(param_1_temp_unc) 
                     mat_param2_systematic[xx,yy] = abs(param_2_temp_unc) 
                     mat_param3_systematic[xx,yy] = abs(param_3_temp_unc)
                     mat_param4_systematic[xx,yy] = abs(param_4_temp_unc) 

                  endif

               endif

            endfor
            endfor

            index = where(mat eq -32768)

            param0_systematic[index] = -32768
            param1_systematic[index] = -32768
            param2_systematic[index] = -32768
            param3_systematic[index] = -32768
            param4_systematic[index] = -32768

            mat_param0_systematic[index] = -32768
            mat_param1_systematic[index] = -32768
            mat_param2_systematic[index] = -32768
            mat_param3_systematic[index] = -32768
            mat_param4_systematic[index] = -32768

            random_unc[index] = -32768
            local_systematic[index] = -32768
            systematic[index] = -32768
            mat_systematic[index] = -32768


;make output file
; Create a new NetCDF file: 
            id = NCDF_CREATE(filename_output, /CLOBBER)
            NCDF_CONTROL, id, /FILL
            
; Make dimensions for NetCDF file
            xid=NCDF_DIMDEF(id, 'x', 360*4) ; dims in longitude direction
            yid=NCDF_DIMDEF(id, 'y', 180*4) ; dimes in latitude direction
            aid=NCDF_DIMDEF(id, 'a', 1)   ; 1-element data set.
            
; Define variables and assign attributes
;----------------------------------------
            NCDF_ATTPUT, id, /GLOBAL, 'Title', 'EUSTACE Surface Air Temperature Estimates' 
            NCDF_ATTPUT, id, /GLOBAL, 'Institution', 'Met Office Hadley Centre, UK' 
            NCDF_ATTPUT, id, /GLOBAL, 'History', 'Produced at Met Office Hadley Centre on '+SYSTIME() 
            NCDF_ATTPUT, id, /GLOBAL, 'References', 'None' 
            NCDF_ATTPUT, id, /GLOBAL, 'Comment', 'Marine surface air temperatures produced using '+$
                         'University of Reading AATSR SST retrievals'
            
; Time
            timeid=NCDF_VARDEF(id, 'time', [aid], /FLOAT)
            NCDF_ATTPUT, id, timeid, 'long_name', 'Time (days)' 
            NCDF_ATTPUT, id, timeid, 'standard_name', 'time' 
            NCDF_ATTPUT, id, timeid, 'units', 'days since 1850-1-1 0:0:0'      
            
; Latitudes
            latid=NCDF_VARDEF(id, 'latitude', [yid], /FLOAT)
            NCDF_ATTPUT, id, latid, 'long_name', 'Centre latitude (deg)' 
            NCDF_ATTPUT, id, latid, 'standard_name', 'latitude' 
            NCDF_ATTPUT, id, latid, 'units', 'degrees_north' 
            NCDF_ATTPUT, id, latid, '_FillValue', geo_fill 
            NCDF_ATTPUT, id, latid, 'valid_min', -90.
            NCDF_ATTPUT, id, latid, 'valid_max', 90.
            
; Longitudes
            lonid=NCDF_VARDEF(id, 'longitude', [xid], /FLOAT)
            NCDF_ATTPUT, id, lonid, 'long_name', 'Centre longitude (deg)' 
            NCDF_ATTPUT, id, lonid, 'standard_name', 'longitude' 
            NCDF_ATTPUT, id, lonid, 'units', 'degrees_east' 
            NCDF_ATTPUT, id, lonid, '_FillValue', geo_fill 
            NCDF_ATTPUT, id, lonid, 'valid_min', -180.
            NCDF_ATTPUT, id, lonid, 'valid_max', 180.
            
; Tmean
            tmean_id=NCDF_VARDEF(id, 'tas', [xid,yid], /LONG)
            NCDF_ATTPUT, id, tmean_id, 'long_name', 'Mean daily surface air temperature (K)'
            NCDF_ATTPUT, id, tmean_id, 'standard_name', 'air_temperature' 
            NCDF_ATTPUT, id, tmean_id, 'units', tunits 
            NCDF_ATTPUT, id, tmean_id, '_FillValue', tfill 
            NCDF_ATTPUT, id, tmean_id, 'valid_min', tvalid_min
            NCDF_ATTPUT, id, tmean_id, 'valid_max', tvalid_max
            NCDF_ATTPUT, id, tmean_id, 'scale_factor', tscale_factor
            NCDF_ATTPUT, id, tmean_id, 'add_offset', tadd_offset
            NCDF_ATTPUT, id, tmean_id, 'coordinates', coordinates
            NCDF_ATTPUT, id, tmean_id, 'cell_methods', 'time: mean'
            
; Tmean random uncertainty
            unc_rand_tmean_id=NCDF_VARDEF(id, 'unc_rand_tas', [xid,yid], /LONG)
            NCDF_ATTPUT, id, unc_rand_tmean_id, 'long_name', 'Random uncertainty on mean daily surface air temperature (K)'
            NCDF_ATTPUT, id, unc_rand_tmean_id, 'standard_name', 'air_temperature' 
            NCDF_ATTPUT, id, unc_rand_tmean_id, 'units', tunits 
            NCDF_ATTPUT, id, unc_rand_tmean_id, '_FillValue', tfill
            NCDF_ATTPUT, id, unc_rand_tmean_id, 'valid_min', uvalid_min
            NCDF_ATTPUT, id, unc_rand_tmean_id, 'valid_max', uvalid_max
            NCDF_ATTPUT, id, unc_rand_tmean_id, 'scale_factor', uscale_factor
            NCDF_ATTPUT, id, unc_rand_tmean_id, 'add_offset', uadd_offset
            NCDF_ATTPUT, id, unc_rand_tmean_id, 'coordinates', coordinates
            NCDF_ATTPUT, id, unc_rand_tmean_id, 'cell_methods', 'time: mean'
            NCDF_ATTPUT, id, unc_rand_tmean_id, 'uncertainty_type', 'uncorrelated'

            
; Tmean locally correlated uncertainty - component 1 from ATSR SST
            unc_corr_tmean_id=NCDF_VARDEF(id, 'unc_corr_tas', [xid,yid], /LONG)
            NCDF_ATTPUT, id, unc_corr_tmean_id, 'long_name', 'Locally correlated uncertainty on mean daily surface air temperature (K) from satellite retrieval'
            NCDF_ATTPUT, id, unc_corr_tmean_id, 'standard_name', 'air_temperature' 
            NCDF_ATTPUT, id, unc_corr_tmean_id, 'units', tunits 
            NCDF_ATTPUT, id, unc_corr_tmean_id, '_FillValue', tfill
            NCDF_ATTPUT, id, unc_corr_tmean_id, 'valid_min', uvalid_min
            NCDF_ATTPUT, id, unc_corr_tmean_id, 'valid_max', uvalid_max
            NCDF_ATTPUT, id, unc_corr_tmean_id, 'scale_factor', uscale_factor
            NCDF_ATTPUT, id, unc_corr_tmean_id, 'add_offset', uadd_offset
            NCDF_ATTPUT, id, unc_corr_tmean_id, 'coordinates', coordinates
            NCDF_ATTPUT, id, unc_corr_tmean_id, 'cell_methods', 'time: mean'
            NCDF_ATTPUT, id, unc_corr_tmean_id, 'uncertainty_type', 'locally correlated'
            NCDF_ATTPUT, id, unc_corr_tmean_id, 'length_scale', 100.0
            NCDF_ATTPUT, id, unc_corr_tmean_id, 'time_scale', 1.0

            
; Tmean locally correlated uncertainty - component 2 from SST-MAT conversion
            unc_corr2_tmean_id=NCDF_VARDEF(id, 'unc_corr2_tas', [xid,yid], /LONG)
            NCDF_ATTPUT, id, unc_corr2_tmean_id, 'long_name', 'Locally correlated uncertainty on mean daily surface air temperature (K)'
            NCDF_ATTPUT, id, unc_corr2_tmean_id, 'standard_name', 'air_temperature'
            NCDF_ATTPUT, id, unc_corr2_tmean_id, 'units', tunits
            NCDF_ATTPUT, id, unc_corr2_tmean_id, '_FillValue', tfill
            NCDF_ATTPUT, id, unc_corr2_tmean_id, 'valid_min', uvalid_min
            NCDF_ATTPUT, id, unc_corr2_tmean_id, 'valid_max', uvalid_max
            NCDF_ATTPUT, id, unc_corr2_tmean_id, 'scale_factor', uscale_factor
            NCDF_ATTPUT, id, unc_corr2_tmean_id, 'add_offset', uadd_offset
            NCDF_ATTPUT, id, unc_corr2_tmean_id, 'coordinates', coordinates
            NCDF_ATTPUT, id, unc_corr2_tmean_id, 'cell_methods', 'time: mean'
            NCDF_ATTPUT, id, unc_corr2_tmean_id, 'uncertainty_type', 'locally correlated'
            NCDF_ATTPUT, id, unc_corr2_tmean_id, 'length_scale', 1200.0
            NCDF_ATTPUT, id, unc_corr2_tmean_id, 'time_scale', 1.0

; Tmean systematic uncertainty - component 1 from ATSR SST
            unc_syst_tmean_id=NCDF_VARDEF(id, 'unc_syst_tas', [xid,yid], /LONG)
            NCDF_ATTPUT, id, unc_syst_tmean_id, 'long_name', 'Systematic uncertainty on mean daily surface air temperature (K)'
            NCDF_ATTPUT, id, unc_syst_tmean_id, 'standard_name', 'air_temperature'
            NCDF_ATTPUT, id, unc_syst_tmean_id, 'units', tunits
            NCDF_ATTPUT, id, unc_syst_tmean_id, '_FillValue', tfill
            NCDF_ATTPUT, id, unc_syst_tmean_id, 'valid_min', uvalid_min
            NCDF_ATTPUT, id, unc_syst_tmean_id, 'valid_max', uvalid_max
            NCDF_ATTPUT, id, unc_syst_tmean_id, 'scale_factor', uscale_factor
            NCDF_ATTPUT, id, unc_syst_tmean_id, 'add_offset', uadd_offset
            NCDF_ATTPUT, id, unc_syst_tmean_id, 'coordinates', coordinates
            NCDF_ATTPUT, id, unc_syst_tmean_id, 'cell_methods', 'time: mean'
            NCDF_ATTPUT, id, unc_syst_tmean_id, 'uncertainty_type', 'systematic'            

; Tmean systematic uncertainty - component 2 from in situ bias uncertainty 
            unc_syst2_tmean_id=NCDF_VARDEF(id, 'unc_syst2_tas', [xid,yid], /LONG)
            NCDF_ATTPUT, id, unc_syst2_tmean_id, 'long_name', 'Systematic uncertainty on mean daily surface air temperature (K)'
            NCDF_ATTPUT, id, unc_syst2_tmean_id, 'standard_name', 'air_temperature' 
            NCDF_ATTPUT, id, unc_syst2_tmean_id, 'units', tunits 
            NCDF_ATTPUT, id, unc_syst2_tmean_id, '_FillValue', tfill
            NCDF_ATTPUT, id, unc_syst2_tmean_id, 'valid_min', uvalid_min
            NCDF_ATTPUT, id, unc_syst2_tmean_id, 'valid_max', uvalid_max
            NCDF_ATTPUT, id, unc_syst2_tmean_id, 'scale_factor', uscale_factor
            NCDF_ATTPUT, id, unc_syst2_tmean_id, 'add_offset', uadd_offset
            NCDF_ATTPUT, id, unc_syst2_tmean_id, 'coordinates', coordinates
            NCDF_ATTPUT, id, unc_syst2_tmean_id, 'cell_methods', 'time: mean'
            NCDF_ATTPUT, id, unc_syst2_tmean_id, 'uncertainty_type', 'systematic'            

; Tmean systematic uncertainty - parameters in climatology 0
            unc_param0_tmean_id=NCDF_VARDEF(id, 'unc_parameter_0_tas', [xid,yid], /LONG)
            NCDF_ATTPUT, id, unc_param0_tmean_id, 'long_name', 'Systematic uncertainty mean offset on mean daily surface air temperature (K)'
            NCDF_ATTPUT, id, unc_param0_tmean_id, 'standard_name', 'air_temperature'
            NCDF_ATTPUT, id, unc_param0_tmean_id, 'units', tunits
            NCDF_ATTPUT, id, unc_param0_tmean_id, '_FillValue', tfill
            NCDF_ATTPUT, id, unc_param0_tmean_id, 'valid_min', uvalid_min
            NCDF_ATTPUT, id, unc_param0_tmean_id, 'valid_max', uvalid_max
            NCDF_ATTPUT, id, unc_param0_tmean_id, 'scale_factor', uscale_factor
            NCDF_ATTPUT, id, unc_param0_tmean_id, 'add_offset', uadd_offset
            NCDF_ATTPUT, id, unc_param0_tmean_id, 'coordinates', coordinates
            NCDF_ATTPUT, id, unc_param0_tmean_id, 'cell_methods', 'time: mean'
            NCDF_ATTPUT, id, unc_param0_tmean_id, 'uncertainty_type', 'systematic'

; Tmean systematic uncertainty - parameters in climatology 1
            unc_param1_tmean_id=NCDF_VARDEF(id, 'unc_parameter_1_tas', [xid,yid], /LONG)
            NCDF_ATTPUT, id, unc_param1_tmean_id, 'long_name', 'Systematic uncertainty first fourier component on mean daily surface air temperature (K)'
            NCDF_ATTPUT, id, unc_param1_tmean_id, 'standard_name', 'air_temperature'
            NCDF_ATTPUT, id, unc_param1_tmean_id, 'units', tunits
            NCDF_ATTPUT, id, unc_param1_tmean_id, '_FillValue', tfill
            NCDF_ATTPUT, id, unc_param1_tmean_id, 'valid_min', uvalid_min
            NCDF_ATTPUT, id, unc_param1_tmean_id, 'valid_max', uvalid_max
            NCDF_ATTPUT, id, unc_param1_tmean_id, 'scale_factor', uscale_factor
            NCDF_ATTPUT, id, unc_param1_tmean_id, 'add_offset', uadd_offset
            NCDF_ATTPUT, id, unc_param1_tmean_id, 'coordinates', coordinates
            NCDF_ATTPUT, id, unc_param1_tmean_id, 'cell_methods', 'time: mean'
            NCDF_ATTPUT, id, unc_param1_tmean_id, 'uncertainty_type', 'systematic'

            
; Tmean systematic uncertainty - parameters in climatology 2
            unc_param2_tmean_id=NCDF_VARDEF(id, 'unc_parameter_2_tas', [xid,yid], /LONG)
            NCDF_ATTPUT, id, unc_param2_tmean_id, 'long_name', 'Systematic uncertainty second fourier component on mean daily surface air temperature (K)'
            NCDF_ATTPUT, id, unc_param2_tmean_id, 'standard_name', 'air_temperature'
            NCDF_ATTPUT, id, unc_param2_tmean_id, 'units', tunits
            NCDF_ATTPUT, id, unc_param2_tmean_id, '_FillValue', tfill
            NCDF_ATTPUT, id, unc_param2_tmean_id, 'valid_min', uvalid_min
            NCDF_ATTPUT, id, unc_param2_tmean_id, 'valid_max', uvalid_max
            NCDF_ATTPUT, id, unc_param2_tmean_id, 'scale_factor', uscale_factor
            NCDF_ATTPUT, id, unc_param2_tmean_id, 'add_offset', uadd_offset
            NCDF_ATTPUT, id, unc_param2_tmean_id, 'coordinates', coordinates
            NCDF_ATTPUT, id, unc_param2_tmean_id, 'cell_methods', 'time: mean'
            NCDF_ATTPUT, id, unc_param2_tmean_id, 'uncertainty_type', 'systematic'

            
; Tmean systematic uncertainty - parameters in climatology 3
            unc_param3_tmean_id=NCDF_VARDEF(id, 'unc_parameter_3_tas', [xid,yid], /LONG)
            NCDF_ATTPUT, id, unc_param3_tmean_id, 'long_name', 'Systematic uncertainty third fourier component on mean daily surface air temperature (K)'
            NCDF_ATTPUT, id, unc_param3_tmean_id, 'standard_name', 'air_temperature'
            NCDF_ATTPUT, id, unc_param3_tmean_id, 'units', tunits
            NCDF_ATTPUT, id, unc_param3_tmean_id, '_FillValue', tfill
            NCDF_ATTPUT, id, unc_param3_tmean_id, 'valid_min', uvalid_min
            NCDF_ATTPUT, id, unc_param3_tmean_id, 'valid_max', uvalid_max
            NCDF_ATTPUT, id, unc_param3_tmean_id, 'scale_factor', uscale_factor
            NCDF_ATTPUT, id, unc_param3_tmean_id, 'add_offset', uadd_offset
            NCDF_ATTPUT, id, unc_param3_tmean_id, 'coordinates', coordinates
            NCDF_ATTPUT, id, unc_param3_tmean_id, 'cell_methods', 'time: mean'
            NCDF_ATTPUT, id, unc_param3_tmean_id, 'uncertainty_type', 'systematic'

            
; Tmean systematic uncertainty - parameters in climatology 4
            unc_param4_tmean_id=NCDF_VARDEF(id, 'unc_parameter_4_tas', [xid,yid], /LONG)
            NCDF_ATTPUT, id, unc_param4_tmean_id, 'long_name', 'Systematic uncertainty fourth fourier component on mean daily surface air temperature (K)'
            NCDF_ATTPUT, id, unc_param4_tmean_id, 'standard_name', 'air_temperature'
            NCDF_ATTPUT, id, unc_param4_tmean_id, 'units', tunits
            NCDF_ATTPUT, id, unc_param4_tmean_id, '_FillValue', tfill
            NCDF_ATTPUT, id, unc_param4_tmean_id, 'valid_min', uvalid_min
            NCDF_ATTPUT, id, unc_param4_tmean_id, 'valid_max', uvalid_max
            NCDF_ATTPUT, id, unc_param4_tmean_id, 'scale_factor', uscale_factor
            NCDF_ATTPUT, id, unc_param4_tmean_id, 'add_offset', uadd_offset
            NCDF_ATTPUT, id, unc_param4_tmean_id, 'coordinates', coordinates
            NCDF_ATTPUT, id, unc_param4_tmean_id, 'cell_methods', 'time: mean'
            NCDF_ATTPUT, id, unc_param4_tmean_id, 'uncertainty_type', 'systematic'
            
; uncertainty in parameters in climatology 0                                                                           
            unc_param0_id=NCDF_VARDEF(id, 'unc_parameter_0', [xid,yid], /LONG)
            NCDF_ATTPUT, id, unc_param0_id, 'long_name', 'Systematic uncertainty mean offset parameter'
            NCDF_ATTPUT, id, unc_param0_id, 'standard_name', 'air_temperature'
            NCDF_ATTPUT, id, unc_param0_id, 'units', tunits
            NCDF_ATTPUT, id, unc_param0_id, '_FillValue', tfill
            NCDF_ATTPUT, id, unc_param0_id, 'valid_min', uvalid_min
            NCDF_ATTPUT, id, unc_param0_id, 'valid_max', uvalid_max
            NCDF_ATTPUT, id, unc_param0_id, 'scale_factor', tscale_factor
            NCDF_ATTPUT, id, unc_param0_id, 'add_offset', uadd_offset
            NCDF_ATTPUT, id, unc_param0_id, 'coordinates', coordinates
            NCDF_ATTPUT, id, unc_param0_id, 'cell_methods', 'time: mean'
            NCDF_ATTPUT, id, unc_param0_id, 'uncertainty_type', 'systematic'
            
; uncertainty - parameters in climatology 1                                                                           
            unc_param1_id=NCDF_VARDEF(id, 'unc_parameter_1', [xid,yid], /LONG)
            NCDF_ATTPUT, id, unc_param1_id, 'long_name', 'Systematic uncertainty first fourier component parameter'
            NCDF_ATTPUT, id, unc_param1_id, 'standard_name', 'air_temperature'
            NCDF_ATTPUT, id, unc_param1_id, 'units', tunits
            NCDF_ATTPUT, id, unc_param1_id, '_FillValue', tfill
            NCDF_ATTPUT, id, unc_param1_id, 'valid_min', uvalid_min
            NCDF_ATTPUT, id, unc_param1_id, 'valid_max', uvalid_max
            NCDF_ATTPUT, id, unc_param1_id, 'scale_factor', tscale_factor
            NCDF_ATTPUT, id, unc_param1_id, 'add_offset', uadd_offset
            NCDF_ATTPUT, id, unc_param1_id, 'coordinates', coordinates
            NCDF_ATTPUT, id, unc_param1_id, 'cell_methods', 'time: mean'
            NCDF_ATTPUT, id, unc_param1_id, 'uncertainty_type', 'systematic'
            
; uncertainty - parameters in climatology 2                                                                           
            unc_param2_id=NCDF_VARDEF(id, 'unc_parameter_2', [xid,yid], /LONG)
            NCDF_ATTPUT, id, unc_param2_id, 'long_name', 'Systematic uncertainty second fourier component parameter'
            NCDF_ATTPUT, id, unc_param2_id, 'standard_name', 'air_temperature'
            NCDF_ATTPUT, id, unc_param2_id, 'units', tunits
            NCDF_ATTPUT, id, unc_param2_id, '_FillValue', tfill
            NCDF_ATTPUT, id, unc_param2_id, 'valid_min', uvalid_min
            NCDF_ATTPUT, id, unc_param2_id, 'valid_max', uvalid_max
            NCDF_ATTPUT, id, unc_param2_id, 'scale_factor', tscale_factor
            NCDF_ATTPUT, id, unc_param2_id, 'add_offset', uadd_offset
            NCDF_ATTPUT, id, unc_param2_id, 'coordinates', coordinates
            NCDF_ATTPUT, id, unc_param2_id, 'cell_methods', 'time: mean'
            NCDF_ATTPUT, id, unc_param2_id, 'uncertainty_type', 'systematic'
            
; uncertainty - parameters in climatology 3                                                                           
            unc_param3_id=NCDF_VARDEF(id, 'unc_parameter_3', [xid,yid], /LONG)
            NCDF_ATTPUT, id, unc_param3_id, 'long_name', 'Systematic uncertainty third fourier component parameter'
            NCDF_ATTPUT, id, unc_param3_id, 'standard_name', 'air_temperature'
            NCDF_ATTPUT, id, unc_param3_id, 'units', tunits
            NCDF_ATTPUT, id, unc_param3_id, '_FillValue', tfill
            NCDF_ATTPUT, id, unc_param3_id, 'valid_min', uvalid_min
            NCDF_ATTPUT, id, unc_param3_id, 'valid_max', uvalid_max
            NCDF_ATTPUT, id, unc_param3_id, 'scale_factor', tscale_factor
            NCDF_ATTPUT, id, unc_param3_id, 'add_offset', uadd_offset
            
; uncertainty - parameters in climatology 4                                                                           
            unc_param4_id=NCDF_VARDEF(id, 'unc_parameter_4', [xid,yid], /LONG)
            NCDF_ATTPUT, id, unc_param4_id, 'long_name', 'Systematic uncertainty fourth fourier component parameter'
            NCDF_ATTPUT, id, unc_param4_id, 'standard_name', 'air_temperature'
            NCDF_ATTPUT, id, unc_param4_id, 'units', tunits
            NCDF_ATTPUT, id, unc_param4_id, '_FillValue', tfill
            NCDF_ATTPUT, id, unc_param4_id, 'valid_min', uvalid_min
            NCDF_ATTPUT, id, unc_param4_id, 'valid_max', uvalid_max
            NCDF_ATTPUT, id, unc_param4_id, 'scale_factor', tscale_factor
            NCDF_ATTPUT, id, unc_param4_id, 'add_offset', uadd_offset
            NCDF_ATTPUT, id, unc_param4_id, 'coordinates', coordinates
            NCDF_ATTPUT, id, unc_param4_id, 'cell_methods', 'time: mean'
            NCDF_ATTPUT, id, unc_param4_id, 'uncertainty_type', 'systematic'


; Put file in data mode: 
            NCDF_CONTROL, id, /ENDEF 
            
; Scale data for file-writing
            invalid = WHERE(mat EQ -32768)

            mat = ROUND(mat/tscale_factor)
            mat[invalid] = tfill

            random_unc = ROUND(random_unc/uscale_factor)
            random_unc[invalid] = tfill

            local_systematic = ROUND(local_systematic/uscale_factor)
            local_systematic[invalid] = tfill

            systematic = ROUND(systematic/uscale_factor)
            systematic[invalid] = tfill

            mat_local_systematic = ROUND(mat_local_systematic/uscale_factor)
            mat_local_systematic[invalid] = tfill
            
            mat_systematic = ROUND(mat_systematic/uscale_factor)
            mat_systematic[invalid] = tfill

            mat_param0_systematic = ROUND(mat_param0_systematic/uscale_factor)
            mat_param0_systematic[invalid] = tfill

            mat_param1_systematic = ROUND(mat_param1_systematic/uscale_factor)
            mat_param1_systematic[invalid] = tfill

            mat_param2_systematic = ROUND(mat_param2_systematic/uscale_factor)
            mat_param2_systematic[invalid] = tfill

            mat_param3_systematic = ROUND(mat_param3_systematic/uscale_factor)
            mat_param3_systematic[invalid] = tfill

            mat_param4_systematic = ROUND(mat_param4_systematic/uscale_factor)
            mat_param4_systematic[invalid] = tfill

            param0_systematic = ROUND(param0_systematic/tscale_factor)
            param0_systematic[invalid] = tfill

            param1_systematic = ROUND(param1_systematic/tscale_factor)
            param1_systematic[invalid] = tfill

            param2_systematic = ROUND(param2_systematic/tscale_factor)
            param2_systematic[invalid] = tfill

            param3_systematic = ROUND(param3_systematic/tscale_factor)
            param3_systematic[invalid] = tfill

            param4_systematic = ROUND(param4_systematic/tscale_factor)
            param4_systematic[invalid] = tfill


; Input data: 
            NCDF_VARPUT, id, timeid, DOUBLE(JULDAY(mmm, ddd, yyy)-JULDAY(1, 1, 1850))
            NCDF_VARPUT, id, latid, REFORM(latitudes) 
            NCDF_VARPUT, id, lonid, REFORM(longitudes)

            NCDF_VARPUT, id, tmean_id, mat

            NCDF_VARPUT, id, unc_rand_tmean_id, random_unc

            NCDF_VARPUT, id, unc_corr_tmean_id, local_systematic
            NCDF_VARPUT, id, unc_syst_tmean_id, systematic             
            
            NCDF_VARPUT, id, unc_corr2_tmean_id, mat_local_systematic
            NCDF_VARPUT, id, unc_syst2_tmean_id, systematic
            
            NCDF_VARPUT, id, unc_param0_tmean_id, mat_param0_systematic
            NCDF_VARPUT, id, unc_param1_tmean_id, mat_param1_systematic
            NCDF_VARPUT, id, unc_param2_tmean_id, mat_param2_systematic
            NCDF_VARPUT, id, unc_param3_tmean_id, mat_param3_systematic
            NCDF_VARPUT, id, unc_param4_tmean_id, mat_param4_systematic
            
            NCDF_VARPUT, id, unc_param0_id, param0_systematic
            NCDF_VARPUT, id, unc_param1_id, param1_systematic
            NCDF_VARPUT, id, unc_param2_id, param2_systematic
            NCDF_VARPUT, id, unc_param3_id, param3_systematic
            NCDF_VARPUT, id, unc_param4_id, param4_systematic
            
; Close the NetCDF file & zip
            NCDF_CLOSE, id

endforeach

end
