; IDL script ported from inside of loop in original model
; This is to be run from a command line like:
;
; idl -e ".run model_ocean_idl.pro" -args {jsonfile}
;
; where {jsonfile} is path to text file in JSON format containing:
;
; { "model": { "filename_model_clim_parameters": ["filenames"],
;              "frac_sat_obs_threshold":float , 
;              "sampling_threshold":float  },
;   "data" :  [{ "date": { "year": yyyy, "month": mm, "day": dd },
;                "dir_input_fvc": "dir"
;                "filename_input_dem": "filename",
;                "filename_input_mask" :["filename north","filename south"],
;                "filename_input_snow": "filename"
;                "filename_input_lst": ["filename day","filename night"],
;                "filename_output_id": "filename" },
;                "dir_output":"dir"},
;		...
;		]
; }     

@model_land_idl_setup_satellite.pro

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
; filename_model_clim_parameters='/gws/nopw/j04/eustace/data/internal/surfaceair_model_parameters/land/MYG_model{1,2,3}_coefs.txt
; frac_sat_obs_threshold=0.2'
; sampling_threshold=3.0'
coef_files = json_model['filename_model_clim_parameters']
frac_sat_obs_thresh = json_model['frac_sat_obs_threshold']
sampling_thresh = json_model['sampling_threshold']

;PRINT, "coef_files =", coef_files
;PRINT, "fsatobs_thresh =", frac_sat_obs_thresh
;PRINT, "sampling_thresh =", sampling_thresh

    

;======================
; Hardcoded settings
;======================

; Systematic uncertainty term value for file.  Value chosen matches input LST.
tmax_sys_value=0.1
tmin_sys_value=0.1

; Threshold for sampling uncertainty.  If sampling uncertainty exceeds this then cell is set
; to invalid because of likely cloud contamination.
IF N_ELEMENTS(sampling_thresh) EQ 1 THEN $
  spl_thresh=sampling_thresh $
ELSE $
  spl_thresh=10.
  
; Threshold for fraction of cloud-free observations, below which observations
; are rejected.
IF N_ELEMENTS(frac_sat_obs_thresh) EQ 1 THEN $
  fsatobs_thresh=frac_sat_obs_thresh $
ELSE $
  fsatobs_thresh=0.1
;PRINT, "fsatobs_thresh =", fsatobs_thresh
;PRINT, "sampling_thresh =", spl_thresh

; Fill values, add offsets and scale factors for LST variables.
; Save time by not reading in for every file.
tsmean_fill=-32768
tsmean_sf=0.005
tsmean_ao=273.15
tsmax_fill=-32768
tsmax_sf=0.005
tsmax_ao=273.15
tsmin_fill=-32768
tsmin_sf=0.005
tsmin_ao=273.15
tsvar_fill=-1
tsvar_sf=5.e-06
tsvar_ao=0.
tsmean_unc_ran_fill=-32768
tsmean_unc_ran_sf=0.001
tsmean_unc_ran_ao=32.767
tsmean_unc_loc_atm_fill=-32768
tsmean_unc_loc_atm_sf=0.001
tsmean_unc_loc_atm_ao=32.767
tsmean_unc_loc_sfc_fill=-32768
tsmean_unc_loc_sfc_sf=0.001
tsmean_unc_loc_sfc_ao=32.767
tsmean_unc_sys_fill=-32768
tsmean_unc_sys_sf=0.001
tsmean_unc_sys_ao=32.767
tsmean_unc_spl_fill=-32768
tsmean_unc_spl_sf=0.001
tsmean_unc_spl_ao=32.767
fvc_sf=0.004
fvc_ao=0.0
fvc_unc_sf=0.005
fvc_unc_ao=0.0

foreach json_data, json_data_list do begin
  
    ; Extract model file details from JSON - example contents should be like:
    ; dir_input_fvc=/gws/nopw/j04/eustace/users/ejn2/Output/L3/FVC/with_unc/global/res_0.25_0.25/
    ; filename_input_dem=/gws/nopw/j04/eustace/users/ejn2/Output/L3/DEM/with_unc/global/DEM_global_0.25_0.25.nc
    ; filename_input_snow=/gws/nopw/j04/eustace/users/ejn2/Output/L3/snow/global/res_0.25_0.25/snow_global_0.25_0.25_YYYYMMDD.nc
    ; filename_input_lst=/gws/nopw/j04/eustace/data/internal/satgrid_lst/res_0.25/satgrid.{day,night}.YYYYMMDD.nc
    ; filename_output='eustace_satellite_4.100'

    json_date = json_data['date']
    yyyy = json_date['year']
    mm = json_date['month']
    dd = json_date['day']

    fvc_dir = json_data['dir_input_fvc']
    dem_file = json_data['filename_input_dem']
    mask_files = json_data['filename_input_mask']
    filename_input_snow = json_data['filename_input_snow']
    filename_input_lst = json_data['filename_input_lst']
    outfile_id = json_data['filename_output_id']
    outdir = json_data['dir_output']

    ;date=yyyy+mm+dd 
    
    
    ;PRINT, "date = ", yyyy,"/",mm,"/",dd
    ;PRINT, "dir_input_fvc = ", fvc_dir
    ;PRINT, "filename_input_dem = ", dem_file 
    ;PRINT, "filename_input_mask = ", mask_files 
    ;PRINT, "filename_input_snow = ", filename_input_snow
    ;PRINT, "filename_input_lst = ", filename_input_lst
    ;PRINT, "filename_output_id = ", outfile_id
    ;PRINT, "dir_output = ", outdir

    ; Idl code extracted from lst_to_lsat_global_coefs.pro, written by Elizabeth Good.
    ; The original code has been modified to process a single day, and with a modified input. 
    ; Now no hardcoded instruction are needed for input data files, since these are stored into 
    ; the input json file created by modela_land.py

    ;======================
    ; Housekeeping
    ;======================
    IF FILE_TEST(outdir, /DIRECTORY) NE 1 THEN $
      FILE_MKDIR, outdir
      
    IF NOT KEYWORD_SET(outfile_id) THEN $
      outfile_id='sat_lsat'

    ;-------------------------
    ; Derivatives from inputs
    ;-------------------------
    ;date_range=[date,date]
    ;yyyy_range=STRMID(date_range, 0, 4)
    yyyy_range=[yyyy,yyyy]
    ;mm_range=STRMID(date_range, 4, 2)
    mm_range=[mm,mm]
    ;dd_range=STRMID(date_range, 6, 2)
    dd_range=[dd,dd]
    ; !!!!!NEED TO WATCH HERE - range may cut off last day!!!!!!!!!!!
    jd_range=JULDAY(FIX(mm_range), FIX(dd_range), FIX(yyyy_range))
    njd=jd_range[1]-jd_range[0]+1
    list_jds=jd_range[0]+LINDGEN(njd)

    ;-------------------------
    ; Output NetCDF file info
    ;-------------------------

    tunits = 'K' 
    tfill=-32767s
    tadd_offset = 273.15
    tscale_factor = 0.005
    ufill=-32767s
    uadd_offset = 32.767
    uscale_factor = 0.001
    coordinates = "lon lat" 
    tcomment = "gridded land surface air temperatures estimated from satellite data" ;

    geo_fill=-32767.
    ;GOTO, test_netcdf

    ;=========================
    ; Get model data
    ;=========================

    ;----------------------------------
    ; Read in coefficient files
    ;----------------------------------

    read_asc, coef_files[0], data, data_names, header, head_str='#'
    ngt_coefs1=data[0:6, 0]
    day_coefs1=data[0:6, 1]
    fcal_tmin1=data[7, 0]
    fcal_tmax1=data[7, 1]

    read_asc, coef_files[1], data, data_names, header, head_str='#'
    ngt_coefs2=data[0:6, 0]
    day_coefs2=data[0:6, 1]
    fcal_tmin2=data[7, 0]
    fcal_tmax2=data[7, 1]

    read_asc, coef_files[2], data, data_names, header, head_str='#'
    ngt_coefs3=data[0:6, 0]
    day_coefs3=data[0:6, 1]
    fcal_tmin3=data[7, 0]
    fcal_tmax3=data[7, 1]

    ;=========================
    ; Get Static DEM data
    ;=========================

    ; Lats and longs (will be same for all data sets, although we will check)
    lat=ncdf_quick_Read(dem_file[0], 'lat')
    lon=ncdf_quick_Read(dem_file[0], 'lon')
    dims=SIZE(lat, /DIMENSIONS)
    nx=dims[0]
    ny=dims[1]

    ; Data
    dem=ncdf_quick_Read(dem_file[0], 'dem')
    dem_ran=MAKE_ARRAY(nx, ny, VALUE=0.)
    dem_loc=MAKE_ARRAY(nx, ny, VALUE=0.)
    dem_sys=MAKE_ARRAY(nx, ny, VALUE=0.)


    ;=========================
    ; List of FVC files
    ;=========================

    list_fvc_files=FILE_SEARCH(fvc_dir[0]+'*')
    fvc_jd=yyyymmdd_to_jd(STRMID(FILE_BASENAME(list_fvc_files), 21, 8))
    test=WHERE(fvc_jd GT MIN(jd_range)-10 AND fvc_jd LT MAX(jd_range)+10, ntest)
    IF ntest EQ 0 THEN $
      MESSAGE, 'Cannot find any FVC files within specified date range'

    ;==================================================
    ; Lat and Long Grid and for SATGRID global data set
    ;==================================================

    dy=ROUND((lat[0, 1]-lat[0, 0])*100)/100.
    dx=ROUND((lon[1, 0]-lon[0, 0])*100)/100.
    generate_lat_long_grid, -90, 90, -180, 180, $
	    dy, dx, lat_grid, lon_grid

    ;===========================
    ; Get cells to mask out LSAT
    ;===========================

    lat_mask1=ncdf_quick_Read(mask_files[0], 'lat')
    lon_mask1=ncdf_quick_Read(mask_files[0], 'lon')
    surface_type_mask1=ncdf_quick_Read(mask_files[0], 'surface_type')
    lat_mask2=ncdf_quick_Read(mask_files[1], 'lat')
    lon_mask2=ncdf_quick_Read(mask_files[1], 'lon')
    surface_type_mask2=ncdf_quick_Read(mask_files[1], 'surface_type')

    ; Masks should span full global longitude, but check.
    IF MAX(lon_mask1) NE MAX(lon_grid) OR $
      MIN(lon_mask1) NE MIN(lon_grid) THEN $
      MESSAGE, 'Mask long grids do not match satellite long grids'
      
    ; Find where mask lats occur in global grid starts.
    match, REFORM(lat_mask1[0, *]), REFORM(lat_grid[0, *]), subs_mask1, subs_lat_grid1
    match, REFORM(lat_mask2[0, *]), REFORM(lat_grid[0, *]), subs_mask2, subs_lat_grid2
    global_mask_temp=BYTARR(nx, ny)
    global_mask=BYTARR(nx, ny)
    global_mask_temp[*, subs_lat_grid1]=surface_type_mask1[*, subs_mask1]
    global_mask_temp[*, subs_lat_grid2]=surface_type_mask2[*, subs_mask2]
    icecap_cells=WHERE(global_mask_temp EQ 1)
    global_mask[icecap_cells]=1

    ;======================
    ; Main processing loop
    ;======================

    ; dummy array for variables we can't read in because file is incomplete.
    nodata=MAKE_ARRAY(nx, ny, /INTEGER, VALUE=tsmean_fill)

    FOR i=jd_range[0], jd_range[1] DO BEGIN
	  
      ; This day in YYYYMMDD format for file finding
      yyyymmdd=jd_to_yyyymmdd(i)
      yyyy=STRMID(yyyymmdd, 0, 4)

      ; Index counter for results arrays.
      loopn=i-jd_range[0]

      ;----------------------
      ; Get all data
      ;----------------------
	  
      ; Read in LST data 
      ;----------------------------------------
      lst_file=[FILE_SEARCH(filename_input_lst[1], COUNT=nfiles_ngt), $
		FILE_SEARCH(filename_input_lst[0], COUNT=nfiles_day)]
      IF nfiles_ngt+nfiles_day GE 1 THEN BEGIN
      ;PRINT, "lst_file", lst_file  
	  ; Read in data from satgrid files.
	  IF nfiles_ngt EQ 1 THEN BEGIN
	    lst_ngt=(ncdf_quick_read(lst_file[0], 'tsmean'))
	    num_lst_ngt_clear=(ncdf_quick_read(lst_file[0], 'ts_number_of_observations'))
	    num_lst_ngt_all=(ncdf_quick_read(lst_file[0], 'total_number_of_observations'))
	    lst_ngt_ran0=(ncdf_quick_read(lst_file[0], 'tsmean_unc_ran'))
	    lst_ngt_sys=(ncdf_quick_read(lst_file[0], 'tsmean_unc_sys'))
	    lst_ngt_loc_atm=(ncdf_quick_read(lst_file[0], 'tsmean_unc_loc_atm'))
	    lst_ngt_loc_sfc=(ncdf_quick_read(lst_file[0], 'tsmean_unc_loc_sfc'))
	    lst_ngt_spl=(ncdf_quick_read(lst_file[0], 'tsmean_unc_spl'))
	  ENDIF ELSE BEGIN
	    lst_ngt=nodata
	    num_lst_ngt_clear=nodata
	    num_lst_ngt_all=nodata
	    lst_ngt_ran0=nodata
	    lst_ngt_sys=nodata
	    lst_ngt_loc_atm=nodata
	    lst_ngt_loc_sfc=nodata
	    lst_ngt_spl=nodata
	  ENDELSE      
	  IF nfiles_day EQ 1 THEN BEGIN
	    lst_day=(ncdf_quick_read(lst_file[1], 'tsmean'))
	    num_lst_day_clear=(ncdf_quick_read(lst_file[1], 'ts_number_of_observations'))
	    num_lst_day_all=(ncdf_quick_read(lst_file[1], 'total_number_of_observations'))
	    lst_day_ran0=(ncdf_quick_read(lst_file[1], 'tsmean_unc_ran'))
	    lst_day_sys=(ncdf_quick_read(lst_file[1], 'tsmean_unc_sys'))
	    lst_day_loc_atm=(ncdf_quick_read(lst_file[1], 'tsmean_unc_loc_atm'))
	    lst_day_loc_sfc=(ncdf_quick_read(lst_file[1], 'tsmean_unc_loc_sfc'))
	    lst_day_spl=(ncdf_quick_read(lst_file[1], 'tsmean_unc_spl'))
	  ENDIF ELSE BEGIN
	    lst_day=nodata
	    num_lst_day_clear=nodata
	    num_lst_day_all=nodata
	    lst_day_ran0=nodata
	    lst_day_sys=nodata
	    lst_day_loc_atm=nodata
	    lst_day_loc_sfc=nodata
	    lst_day_spl=nodata
	  ENDELSE

	  ; Calculate observation fraction.
	  obs_frac_ngt=MAKE_ARRAY(nx, ny, VALUE=0.)
	  obs_frac_day=MAKE_ARRAY(nx, ny, VALUE=0.)
	  valid=WHERE(num_lst_ngt_all GT 0, nvalid)
	  IF nvalid GT 0 THEN $
	    obs_frac_ngt[valid]=num_lst_ngt_clear[valid]/FLOAT(num_lst_ngt_all[valid])
	  valid=WHERE(num_lst_day_all GT 0, nvalid)
	  IF nvalid GT 0 THEN $
	    obs_frac_day[valid]=num_lst_day_clear[valid]/FLOAT(num_lst_day_all[valid])

	  ; Scale and set invalids to Nans: Sampling uncertainty night
	  invalid=WHERE(lst_ngt_spl EQ tsmean_unc_spl_fill) 
	  lst_ngt_spl=lst_ngt_spl*tsmean_unc_spl_sf+tsmean_unc_spl_ao
	  lst_ngt_spl[invalid]=!VALUES.F_NAN
	  invalid=WHERE(lst_ngt_spl LT 0, ninvalid)
	  IF ninvalid GT 0 THEN $
	    lst_ngt_spl[invalid]=0
	  
	  ; Scale and set invalids to Nans: Sampling uncertainty day
	  invalid=WHERE(lst_day_spl EQ tsmean_unc_spl_fill) 
	  lst_day_spl=(lst_day_spl*tsmean_unc_spl_sf+tsmean_unc_spl_ao)
	  lst_day_spl[invalid]=!VALUES.F_NAN
	  invalid=WHERE(lst_day_spl LT 0, ninvalid)
	  IF ninvalid GT 0 THEN $
	    lst_day_spl[invalid]=0
	  
	  ; Scale and set invalids to Nans: LST night
	  invalid=WHERE(lst_ngt EQ tsmean_fill OR lst_ngt_spl GT spl_thresh OR $
                        obs_frac_ngt LT fsatobs_thresh OR global_mask EQ 1) 
	  lst_ngt=lst_ngt*tsmean_sf+tsmean_ao
	  lst_ngt[invalid]=!VALUES.F_NAN
	  
	  ; Scale and set invalids to Nans: LST day
	  invalid=WHERE(lst_day EQ tsmean_fill OR lst_day_spl GT spl_thresh OR $
                        obs_frac_day LT fsatobs_thresh OR global_mask EQ 1) 
	  lst_day=(lst_day*tsmean_sf+tsmean_ao)
	  lst_day[invalid]=!VALUES.F_NAN

	  ; Scale and set invalids to Nans: Random uncertainty night
	  invalid=WHERE(lst_ngt_ran0 EQ tsmean_unc_ran_fill OR lst_ngt_spl GT spl_thresh OR $
                        obs_frac_ngt LT fsatobs_thresh OR global_mask EQ 1) 
	  lst_ngt_ran0=lst_ngt_ran0*tsmean_unc_ran_sf+tsmean_unc_ran_ao
	  lst_ngt_ran0[invalid]=!VALUES.F_NAN
	  
	  ; Scale and set invalids to Nans: Random uncertainty day
	  invalid=WHERE(lst_day_ran0 EQ tsmean_unc_ran_fill OR lst_day_spl GT spl_thresh OR $
                        obs_frac_day LT fsatobs_thresh  OR global_mask EQ 1) 
	  lst_day_ran0=lst_day_ran0*tsmean_unc_ran_sf+tsmean_unc_ran_ao
	  lst_day_ran0[invalid]=!VALUES.F_NAN

	  ; Scale and set invalids to Nans: Locally correlated uncertainty night atm
	  invalid=WHERE(lst_ngt_loc_atm EQ tsmean_unc_loc_atm_fill OR lst_ngt_spl GT spl_thresh $
                        OR obs_frac_ngt LT fsatobs_thresh OR global_mask EQ 1) 
	  lst_ngt_loc_atm=lst_ngt_loc_atm*tsmean_unc_loc_atm_sf+tsmean_unc_loc_atm_ao
	  lst_ngt_loc_atm[invalid]=!VALUES.F_NAN
	  
	  ; Scale and set invalids to Nans: locally correlated uncertainty day atm
	  invalid=WHERE(lst_day_loc_atm EQ tsmean_unc_loc_atm_fill OR lst_day_spl GT spl_thresh $
                        OR obs_frac_day LT fsatobs_thresh OR global_mask EQ 1) 
	  lst_day_loc_atm=lst_day_loc_atm*tsmean_unc_loc_atm_sf+tsmean_unc_loc_atm_ao
	  lst_day_loc_atm[invalid]=!VALUES.F_NAN

	  ; Scale and set invalids to Nans: Locally correlated uncertainty night sfc
	  invalid=WHERE(lst_ngt_loc_sfc EQ tsmean_unc_loc_sfc_fill OR lst_ngt_spl GT spl_thresh OR $
                        obs_frac_ngt LT fsatobs_thresh OR global_mask EQ 1) 
	  lst_ngt_loc_sfc=lst_ngt_loc_sfc*tsmean_unc_loc_sfc_sf+tsmean_unc_loc_sfc_ao
	  lst_ngt_loc_sfc[invalid]=!VALUES.F_NAN
	  
	  ; Scale and set invalids to Nans: locally correlated uncertainty day sfc
	  invalid=WHERE(lst_day_loc_sfc EQ tsmean_unc_loc_sfc_fill OR lst_day_spl GT spl_thresh $
                        OR obs_frac_day LT fsatobs_thresh OR global_mask EQ 1) 
	  lst_day_loc_sfc=lst_day_loc_sfc*tsmean_unc_loc_sfc_sf+tsmean_unc_loc_sfc_ao
	  lst_day_loc_sfc[invalid]=!VALUES.F_NAN

	  ; Scale and set invalids to Nans: Systematic uncertainty night
	  invalid=WHERE(lst_ngt_sys EQ tsmean_unc_sys_fill OR lst_ngt_spl GT spl_thresh OR $
                        obs_frac_ngt LT fsatobs_thresh OR global_mask EQ 1) 
	  lst_ngt_sys=lst_ngt_sys*tsmean_unc_sys_sf+tsmean_unc_sys_ao
	  lst_ngt_sys[invalid]=!VALUES.F_NAN
	  
	  ; Scale and set invalids to Nans: Systematic uncertainty day
	  invalid=WHERE(lst_day_sys EQ tsmean_unc_sys_fill OR lst_day_spl GT spl_thresh OR $
                        obs_frac_day LT fsatobs_thresh OR global_mask EQ 1) 
	  lst_day_sys=lst_day_sys*tsmean_unc_sys_sf+tsmean_unc_sys_ao
	  lst_day_sys[invalid]=!VALUES.F_NAN

	  ; Add sampling and ran0 uncertainties in quadrature
	  lst_day_ran=SQRT(lst_day_ran0^2+lst_day_spl^2)
	  lst_ngt_ran=SQRT(lst_ngt_ran0^2+lst_ngt_spl^2)
	  
	  ; Flag variable to indicate we have valid LST data
	  lst_available=1

      ENDIF ELSE $
	  lst_available=0
      ;PRINT, 'lst_available', lst_available

      ; Read in FVC data 
      ;----------------------------------------
      IF MIN(ABS(i-fvc_jd)) LT 10 THEN BEGIN
	  ; Find FVC file with closest date to the data currently being processed.
	  best_fvc_file=WHERE(ABS(i-fvc_jd) EQ MIN(ABS(i-fvc_jd)), ntest)
	  fvc_file=list_fvc_files[best_fvc_file[0]]
	  file_ext=STRMID(fvc_file, 1, 2, /REVERSE_OFFSET)
	  IF file_ext EQ 'gz' THEN BEGIN
	    SPAWN, 'gunzip '+fvc_file
	    fvc_file=STRMID(fvc_file, 0, STRLEN(fvc_file)-3)
	  ENDIF
	  
	  ; Read in FVC and related parameters
	  fvc=FLOAT(ncdf_quick_read(fvc_file, 'fvc'))
	  fvc_unc_mean=FLOAT(ncdf_quick_read(fvc_file, 'fvc_unc'))
	  fvc_unc_max=FLOAT(ncdf_quick_read(fvc_file, 'fvc_unc_max'))
	  
	  ; Convert to floats and missings to NANs
	  ; Define uncertainties.
	  invalid=WHERE(fvc EQ 255 OR fvc_unc_mean EQ 255)
	  fvc[invalid]=!VALUES.F_NAN
	  fvc_unc_mean[invalid]=!VALUES.F_NAN
	  fvc_unc_max[invalid]=!VALUES.F_NAN
	  fvc=fvc*fvc_sf-fvc_ao
	  fvc_loc=fvc_unc_mean*fvc_unc_sf-fvc_unc_ao
	  fvc_ran=(fvc_unc_max*fvc_unc_sf-fvc_unc_ao)-fvc_loc
	  fvc_sys=MAKE_ARRAY(nx, ny, VALUE=0.004)       
	  IF file_ext EQ 'gz' THEN $
	    SPAWN, 'gzip '+fvc_file
	    
      ENDIF ELSE BEGIN
          fvc=MAKE_ARRAY(nx,ny,VALUE=!VALUES.F_NAN)
	  fvc_loc=MAKE_ARRAY(nx,ny,VALUE=!VALUES.F_NAN)
	  fvc_ran=MAKE_ARRAY(nx,ny,VALUE=!VALUES.F_NAN)
	  fvc_sys=MAKE_ARRAY(nx,ny,VALUE=!VALUES.F_NAN)
      ENDELSE
      ;PRINT, 'fvc_available', fvc_available

      ; Read in snow data
      ;------------------

      snow_file=FILE_SEARCH(filename_input_snow, COUNT=nfiles)
      ;PRINT, 'n_files', nfiles
      IF nfiles EQ 1 THEN BEGIN
	  file_ext=STRMID(snow_file, 1, 2, /REVERSE_OFFSET)
	  IF file_ext EQ 'gz' THEN BEGIN
	    SPAWN, 'gunzip '+snow_file
	    snow_file=STRMID(snow_file, 0, STRLEN(snow_file)-3)
	  ENDIF
	  
	  ; Read in data and associated parameters
	  snow=FLOAT(ncdf_quick_read(snow_file, 'snow'))
	  
	  ; Convert missing data to NANs.
	  ; Define uncertainties.
	  invalid=WHERE(snow EQ 255)
	  snow[invalid]=!VALUES.F_NAN
	  snow_ran=MAKE_ARRAY(nx, ny, VALUE=0.)
	  snow_loc=MAKE_ARRAY(nx, ny, VALUE=0.)
	  snow_sys=MAKE_ARRAY(nx, ny, VALUE=0.)       
	  IF file_ext EQ 'gz' THEN $
	    SPAWN, 'gzip '+snow_file
	  
      ENDIF ELSE BEGIN
          snow=MAKE_ARRAY(nx, ny, VALUE=!VALUES.F_NAN)
          snow_ran=MAKE_ARRAY(nx, ny, VALUE=!VALUES.F_NAN)
          snow_loc=MAKE_ARRAY(nx, ny, VALUE=!VALUES.F_NAN)
          snow_sys=MAKE_ARRAY(nx, ny, VALUE=!VALUES.F_NAN) 
      ENDELSE
      ;PRINT, 'snow_available', snow_available
      ;PRINT, 'snow', snow
      ; Set up arrays to write to file, but do not proceed further
      ; unless we have some valid LST data.
      tmax=MAKE_ARRAY(nx, ny, VALUE=!VALUES.F_NAN)
      tmin=MAKE_ARRAY(nx, ny, VALUE=!VALUES.F_NAN)
      tmax_ran=MAKE_ARRAY(nx, ny, VALUE=!VALUES.F_NAN)
      tmin_ran=MAKE_ARRAY(nx, ny, VALUE=!VALUES.F_NAN)
      tmax_loc_atm=MAKE_ARRAY(nx, ny, VALUE=!VALUES.F_NAN)
      tmin_loc_atm=MAKE_ARRAY(nx, ny, VALUE=!VALUES.F_NAN)
      tmax_loc_sfc=MAKE_ARRAY(nx, ny, VALUE=!VALUES.F_NAN)
      tmin_loc_sfc=MAKE_ARRAY(nx, ny, VALUE=!VALUES.F_NAN)
      tmax_sys=MAKE_ARRAY(nx, ny, VALUE=!VALUES.F_NAN)
      tmin_sys=MAKE_ARRAY(nx, ny, VALUE=!VALUES.F_NAN)
      tmax_model_num=MAKE_ARRAY(nx, ny, /BYTE, VALUE=0)
      tmin_model_num=MAKE_ARRAY(nx, ny, /BYTE, VALUE=0)
      IF lst_available EQ 1 THEN BEGIN
      
	  ; Get solar zenith angles
	  ;------------------------
	  ; Only calculate for where we have either LST day or LST night.
	  ; We will only need them for these locations.
	  doy2d=MAKE_ARRAY(nx, ny, VALUE=date_to_doy(yyyymmdd))
	  yyyy2d=MAKE_ARRAY(nx, ny, VALUE=FIX(yyyy))
	  noon2d=FLTARR(nx, ny)
	  valid=WHERE(FINITE(lst_day) EQ 1 OR FINITE(lst_ngt) EQ 1, nvalid)
	  FOR j=0, nvalid-1 DO $  
	      noon2d[valid[j]]=(sunrise(doy2d[valid[j]], yyyy2d[valid[j]], $
			    lat=lat_grid[valid[j]], lon=lon_grid[valid[j]]))[2]
	  test=WHERE(noon2d LT 0, ntest)
	  IF ntest GT 0 THEN $
	    noon2d[test]=noon2d[test]+24.
	  msec=noon2d*60*60*1000   
	  lonlat2solz, lon_grid[valid], lat_grid[valid], yyyy2d[valid], doy2d[valid], msec[valid], solz_valid, sola_valid
	  solz=FLTARR(nx, ny)
	  solz[valid]=solz_valid
      
	  ; Define uncertainties.
	  solz_ran=MAKE_ARRAY(nx, ny, VALUE=0.)
	  solz_loc=MAKE_ARRAY(nx, ny, VALUE=0.)
	  solz_sys=MAKE_ARRAY(nx, ny, VALUE=0.)       

	  ;----------------------------------------------
	  ; LSAT calculation with uncertainties
	  ;----------------------------------------------
	  
	  ; Apply model coefficients to observational data.  Do for model 3 first, then overwrite
	  ; with model 2, then model 1.  This means we have model 1 wherever possible, followed
	  ; by model 2 then model 3 where model 1 and 2 cannot be achieved.
	  
	  FOR modelv=1, 2 DO BEGIN
	    
	    CASE 1 OF 
		modelv EQ 0: BEGIN
		      ngt_coefs=ngt_coefs3
		      day_coefs=day_coefs3
		      fcal_tmin=fcal_tmin3
		      fcal_tmax=fcal_tmax3
		      this_model_num=3
		  END
		modelv EQ 1: BEGIN
		      ngt_coefs=ngt_coefs2
		      day_coefs=day_coefs2
		      fcal_tmin=fcal_tmin2
		      fcal_tmax=fcal_tmax2
		      this_model_num=2
		  END
		modelv EQ 2: BEGIN
		      ngt_coefs=ngt_coefs1
		      day_coefs=day_coefs1
		      fcal_tmin=fcal_tmin1
		      fcal_tmax=fcal_tmax1
		      this_model_num=1
		  END
	    ELSE:
	    ENDCASE
	  
	    ; LSAT calculation
	    ;------------------
	    ; Make a temporary copy of the variables used in the LSAT prediction so that we can
	    ; set all these values to zero if they are not used. Currently missing data are set
	    ; to NANs so if one parameter is missing, there will be a NAN in the LSAT calculation.
	    ; Note conversion of LSTs to deg C as this is what coefficients are trained on.
	  
	    ; Tmax
	    lst_day_temp=lst_day-273.15
	    lst_ngt_temp=lst_ngt-273.15
	    fvc_temp=fvc
	    dem_temp=dem
	    solz_temp=solz
	    snow_temp=snow      
	    IF day_coefs[1] EQ 0. THEN lst_day_temp[*]=0.
	    IF day_coefs[2] EQ 0. THEN lst_ngt_temp[*]=0.
	    IF day_coefs[3] EQ 0. THEN fvc_temp[*]=0.
	    IF day_coefs[4] EQ 0. THEN dem_temp[*]=0.
	    IF day_coefs[5] EQ 0. THEN solz_temp[*]=0.
	    IF day_coefs[6] EQ 0. THEN snow_temp[*]=0.
	    tmax_temp=day_coefs[0]+day_coefs[1]*lst_day_temp+day_coefs[2]*lst_ngt_temp+day_coefs[3]*fvc_temp+$
		  day_coefs[4]*dem_temp+day_coefs[5]*solz_temp+day_coefs[6]*snow_temp
	    tmax_temp=tmax_temp+273.15
	  
	    ; Tmin
	    lst_day_temp=lst_day-273.15
	    lst_ngt_temp=lst_ngt-273.15
	    fvc_temp=fvc
	    dem_temp=dem
	    solz_temp=solz
	    snow_temp=snow      
	    IF ngt_coefs[1] EQ 0. THEN lst_day_temp[*]=0.
	    IF ngt_coefs[2] EQ 0. THEN lst_ngt_temp[*]=0.
	    IF ngt_coefs[3] EQ 0. THEN fvc_temp[*]=0.
	    IF ngt_coefs[4] EQ 0. THEN dem_temp[*]=0.
	    IF ngt_coefs[5] EQ 0. THEN solz_temp[*]=0.
	    IF ngt_coefs[6] EQ 0. THEN snow_temp[*]=0.
	    tmin_temp=ngt_coefs[0]+ngt_coefs[1]*lst_day_temp+ngt_coefs[2]*lst_ngt_temp+ngt_coefs[3]*fvc_temp+$
		  ngt_coefs[4]*dem_temp+ngt_coefs[5]*solz_temp+ngt_coefs[6]*snow_temp
	    tmin_temp=tmin_temp+273.15

	    ; Uncertainties
	    ;--------------------   
	    ; Method: Finn's uncertainty document.  Again, set parameters we are not actually using
	    ; to zero so that they don't cause calculated results to be NANs where they are not actually used.
	  
	    ; Tmax Random
	    lst_day_ran_temp=lst_day_ran
	    lst_ngt_ran_temp=lst_ngt_ran
	    fvc_ran_temp=fvc_ran
	    dem_ran_temp=dem_ran
	    solz_ran_temp=solz_ran
	    snow_ran_temp=snow_ran    
	    IF day_coefs[1] EQ 0. THEN lst_day_ran_temp[*]=0.
	    IF day_coefs[2] EQ 0. THEN lst_ngt_ran_temp[*]=0.
	    IF day_coefs[3] EQ 0. THEN fvc_ran_temp[*]=0.
	    IF day_coefs[4] EQ 0. THEN dem_ran_temp[*]=0.
	    IF day_coefs[5] EQ 0. THEN solz_ran_temp[*]=0.
	    IF day_coefs[6] EQ 0. THEN snow_ran_temp[*]=0.
	    tmax_ran_temp=SQRT(day_coefs[1]^2*lst_day_ran_temp^2+day_coefs[2]^2*lst_ngt_ran_temp^2+day_coefs[3]^2*fvc_ran_temp^2+$
		  day_coefs[4]^2*dem_ran_temp^2+day_coefs[5]^2*solz_ran_temp^2+day_coefs[6]^2*snow_ran_temp^2)
	  
	    ; Tmax Locally correlated - ATM
	    lst_day_loc_temp=lst_day_loc_atm
	    lst_ngt_loc_temp=lst_ngt_loc_atm
	    fvc_loc_temp=fvc_loc
	    dem_loc_temp=dem_loc
	    solz_loc_temp=solz_loc
	    snow_loc_temp=snow_loc    
	    IF day_coefs[1] EQ 0. THEN lst_day_loc_temp[*]=0.
	    IF day_coefs[2] EQ 0. THEN lst_ngt_loc_temp[*]=0.
	    IF day_coefs[3] EQ 0. THEN fvc_loc_temp[*]=0.
	    IF day_coefs[4] EQ 0. THEN dem_loc_temp[*]=0.
	    IF day_coefs[5] EQ 0. THEN solz_loc_temp[*]=0.
	    IF day_coefs[6] EQ 0. THEN snow_loc_temp[*]=0.
            tmax_loc_atm_temp=SQRT(day_coefs[1]^2*lst_day_loc_temp^2+day_coefs[2]^2*lst_ngt_loc_temp^2+$
                                   day_coefs[4]^2*dem_loc_temp^2+day_coefs[5]^2*solz_loc_temp^2+fcal_tmax^2)

	    ; Tmax Locally correlated - SFC
	    lst_day_loc_temp=lst_day_loc_sfc
	    lst_ngt_loc_temp=lst_ngt_loc_sfc
	    fvc_loc_temp=fvc_loc
	    dem_loc_temp=dem_loc
	    solz_loc_temp=solz_loc
	    snow_loc_temp=snow_loc    
	    IF day_coefs[1] EQ 0. THEN lst_day_loc_temp[*]=0.
	    IF day_coefs[2] EQ 0. THEN lst_ngt_loc_temp[*]=0.
	    IF day_coefs[3] EQ 0. THEN fvc_loc_temp[*]=0.
	    IF day_coefs[4] EQ 0. THEN dem_loc_temp[*]=0.
	    IF day_coefs[5] EQ 0. THEN solz_loc_temp[*]=0.
	    IF day_coefs[6] EQ 0. THEN snow_loc_temp[*]=0.
            tmax_loc_sfc_temp=SQRT(day_coefs[1]^2*lst_day_loc_temp^2+day_coefs[2]^2*lst_ngt_loc_temp^2+$
                                   day_coefs[3]^2*fvc_loc_temp^2+day_coefs[6]^2*snow_loc_temp^2)

	    ; Tmin Random
	    lst_day_ran_temp=lst_day_ran
	    lst_ngt_ran_temp=lst_ngt_ran
	    fvc_ran_temp=fvc_ran
	    dem_ran_temp=dem_ran
	    solz_ran_temp=solz_ran
	    snow_ran_temp=snow_ran    
	    IF ngt_coefs[1] EQ 0. THEN lst_day_ran_temp[*]=0.
	    IF ngt_coefs[2] EQ 0. THEN lst_ngt_ran_temp[*]=0.
	    IF ngt_coefs[3] EQ 0. THEN fvc_ran_temp[*]=0.
	    IF ngt_coefs[4] EQ 0. THEN dem_ran_temp[*]=0.
	    IF ngt_coefs[5] EQ 0. THEN solz_ran_temp[*]=0.
	    IF ngt_coefs[6] EQ 0. THEN snow_ran_temp[*]=0.
	    tmin_ran_temp=SQRT(ngt_coefs[1]^2*lst_day_ran_temp^2+ngt_coefs[2]^2*lst_ngt_ran_temp^2+ngt_coefs[3]^2*fvc_ran_temp^2+$
		  ngt_coefs[4]^2*dem_ran_temp^2+ngt_coefs[5]^2*solz_ran_temp^2+ngt_coefs[6]^2*snow_ran_temp^2)
	  
	    ; Tmin Locally correlated - ATM
	    lst_day_loc_temp=lst_day_loc_atm
	    lst_ngt_loc_temp=lst_ngt_loc_atm
	    fvc_loc_temp=fvc_loc
	    dem_loc_temp=dem_loc
	    solz_loc_temp=solz_loc
	    snow_loc_temp=snow_loc    
	    IF ngt_coefs[1] EQ 0. THEN lst_day_loc_temp[*]=0.
	    IF ngt_coefs[2] EQ 0. THEN lst_ngt_loc_temp[*]=0.
	    IF ngt_coefs[3] EQ 0. THEN fvc_loc_temp[*]=0.
	    IF ngt_coefs[4] EQ 0. THEN dem_loc_temp[*]=0.
	    IF ngt_coefs[5] EQ 0. THEN solz_loc_temp[*]=0.
	    IF ngt_coefs[6] EQ 0. THEN snow_loc_temp[*]=0.
            tmin_loc_atm_temp=SQRT(ngt_coefs[1]^2*lst_day_loc_temp^2+ngt_coefs[2]^2*lst_ngt_loc_temp^2+$
                                   ngt_coefs[4]^2*dem_loc_temp^2+ngt_coefs[5]^2*solz_loc_temp^2+fcal_tmin^2)
 
	    ; Tmin Locally correlated - SFC
	    lst_day_loc_temp=lst_day_loc_sfc
	    lst_ngt_loc_temp=lst_ngt_loc_sfc
	    fvc_loc_temp=fvc_loc
	    dem_loc_temp=dem_loc
	    solz_loc_temp=solz_loc
	    snow_loc_temp=snow_loc    
	    IF ngt_coefs[1] EQ 0. THEN lst_day_loc_temp[*]=0.
	    IF ngt_coefs[2] EQ 0. THEN lst_ngt_loc_temp[*]=0.
	    IF ngt_coefs[3] EQ 0. THEN fvc_loc_temp[*]=0.
	    IF ngt_coefs[4] EQ 0. THEN dem_loc_temp[*]=0.
	    IF ngt_coefs[5] EQ 0. THEN solz_loc_temp[*]=0.
	    IF ngt_coefs[6] EQ 0. THEN snow_loc_temp[*]=0.
            tmin_loc_sfc_temp=SQRT(ngt_coefs[1]^2*lst_day_loc_temp^2+ngt_coefs[2]^2*lst_ngt_loc_temp^2+$
                                   ngt_coefs[3]^2*fvc_loc_temp^2+ngt_coefs[6]^2*snow_loc_temp^2)
	  
	    ; For the present, set systematic uncertainty to be the same everywhere.
	    tmax_sys_temp=MAKE_ARRAY(nx, ny, VALUE=tmax_sys_value)
	    tmin_sys_temp=MAKE_ARRAY(nx, ny, VALUE=tmin_sys_value)
	    
	    ; Assign data to final arrays
	    ;----------------------------
	    valid_tmax=WHERE(FINITE(tmax_temp) EQ 1 AND $
			      FINITE(tmax_ran_temp) EQ 1 AND $
			      FINITE(tmax_loc_atm_temp) EQ 1 AND $
			      FINITE(tmax_loc_sfc_temp) EQ 1, nvalid_tmax)
	    IF nvalid_tmax GT 0 THEN BEGIN
		tmax[valid_tmax]=tmax_temp[valid_tmax]
		tmax_ran[valid_tmax]=tmax_ran_temp[valid_tmax]
		tmax_loc_atm[valid_tmax]=tmax_loc_atm_temp[valid_tmax]
		tmax_loc_sfc[valid_tmax]=tmax_loc_sfc_temp[valid_tmax]
		tmax_sys[valid_tmax]=tmax_sys_temp[valid_tmax]
		tmax_model_num[valid_tmax]=this_model_num
	    ENDIF
	    valid_tmin=WHERE(FINITE(tmin_temp) EQ 1 AND $
			      FINITE(tmin_ran_temp) EQ 1 AND $
			      FINITE(tmin_loc_atm_temp) EQ 1 AND $
			      FINITE(tmin_loc_sfc_temp) EQ 1, nvalid_tmin)
	    IF nvalid_tmin GT 0 THEN BEGIN
		tmin[valid_tmin]=tmin_temp[valid_tmin]
		tmin_ran[valid_tmin]=tmin_ran_temp[valid_tmin]
		tmin_loc_atm[valid_tmin]=tmin_loc_atm_temp[valid_tmin]
		tmin_loc_sfc[valid_tmin]=tmin_loc_sfc_temp[valid_tmin]
		tmin_sys[valid_tmin]=tmin_sys_temp[valid_tmin]
		tmin_model_num[valid_tmin]=this_model_num
	    ENDIF                    
	  
	  ENDFOR

      ENDIF
	  
      ;-------------------------
      ; Write to NetCDF file
      ;--------------------------

      ; Create a new NetCDF file: 
      outfile=outdir+outfile_id;+'.nc';'_'+yyyymmdd+'.nc'
      id = NCDF_CREATE(outfile, /CLOBBER,  /NETCDF4_FORMAT)
      NCDF_CONTROL, id, /FILL

      ; Make dimensions for NetCDF file
      xid=NCDF_DIMDEF(id, 'lon', nx) ; dims in longitude direction
      yid=NCDF_DIMDEF(id, 'lat', ny) ; dims in latitude direction
      aid=NCDF_DIMDEF(id, 'time', /UNLIMITED)   ; Unlimited in time dimension.

      ; Define variables and assign attributes
      ;----------------------------------------

      NCDF_ATTPUT, id, /GLOBAL, 'Title', 'EUSTACE Surface Air Temperature Estimates', /CHAR 
      NCDF_ATTPUT, id, /GLOBAL, 'Institution', 'Met Office Hadley Centre, UK' , /CHAR 
      NCDF_ATTPUT, id, /GLOBAL, 'History', 'Produced at Met Office Hadley Centre on '+SYSTIME() , /CHAR 
      NCDF_ATTPUT, id, /GLOBAL, 'References', 'None' , /CHAR 
      NCDF_ATTPUT, id, /GLOBAL, 'Comment', 'Land surface air temperatures produced using '+$
				'GlobTemperature MODIS Aqua LST retrievals (MYG product). '+$
				'The tmin and tmax fields are produced using three model '+$
				'variants.  The model used to produce each T2m estimate is '+$
				'provided in file.  Model 1 is the best model, followed by '+$
				'model 2 or model 3(mutually exclusive).  Uncertainties should generally be '+$
				'larger for models 2 and 3 compared with model 1 predicted T2m. '+$
				'Sampling uncertainty threshold used to reject grid-cell LSTs is '+$
				num_to_string(spl_thresh)+'. Fraction of valid LSTs threshold used '+$
				'to reject grid-cell LSTs is '+num_to_string(fsatobs_thresh)+'.', /CHAR 
	  
      ; Time
      timeid=NCDF_VARDEF(id, 'time', [aid], /LONG, GZIP=9)
      NCDF_ATTPUT, id, timeid, 'long_name', 'Time (days)' , /CHAR 
      NCDF_ATTPUT, id, timeid, 'standard_name', 'time', /CHAR 
      NCDF_ATTPUT, id, timeid, 'units', 'days since 1850-1-1 0:0:0', /CHAR      

      ; Latitudes
      latid=NCDF_VARDEF(id, 'lat', [yid], /FLOAT, GZIP=9)
      NCDF_ATTPUT, id, latid, 'long_name', 'Centre latitude (deg)' , /CHAR
      NCDF_ATTPUT, id, latid, 'standard_name', 'latitude' , /CHAR
      NCDF_ATTPUT, id, latid, 'units', 'degrees_north', /CHAR
      NCDF_ATTPUT, id, latid, '_FillValue', geo_fill 
      NCDF_ATTPUT, id, latid, 'valid_min', -90.
      NCDF_ATTPUT, id, latid, 'valid_max', 90.

      ; Longitudes
      lonid=NCDF_VARDEF(id, 'lon', [xid], /FLOAT, GZIP=9)
      NCDF_ATTPUT, id, lonid, 'long_name', 'Centre longitude (deg)', /CHAR 
      NCDF_ATTPUT, id, lonid, 'standard_name', 'longitude', /CHAR
      NCDF_ATTPUT, id, lonid, 'units', 'degrees_east', /CHAR
      NCDF_ATTPUT, id, lonid, '_FillValue', geo_fill 
      NCDF_ATTPUT, id, lonid, 'valid_min', -180.
      NCDF_ATTPUT, id, lonid, 'valid_max', 180.

      ; Tmin
      tmin_id=NCDF_VARDEF(id, 'tasmin', [aid, xid,yid], /SHORT, GZIP=9)
      NCDF_ATTPUT, id, tmin_id, 'long_name', 'Minimum daily surface air temperature (K)', /CHAR
      NCDF_ATTPUT, id, tmin_id, 'standard_name', 'air_temperature', /CHAR 
      NCDF_ATTPUT, id, tmin_id, 'units', tunits, /CHAR 
      NCDF_ATTPUT, id, tmin_id, '_FillValue', tfill 
      NCDF_ATTPUT, id, tmin_id, 'scale_factor', tscale_factor
      NCDF_ATTPUT, id, tmin_id, 'add_offset', tadd_offset
      NCDF_ATTPUT, id, tmin_id, 'coordinates', coordinates, /CHAR

      ; Tmax
      tmax_id=NCDF_VARDEF(id, 'tasmax', [aid, xid,yid], /SHORT, GZIP=9)
      NCDF_ATTPUT, id, tmax_id, 'long_name', 'Maximum daily surface air temperature (K)', /CHAR
      NCDF_ATTPUT, id, tmax_id, 'standard_name', 'air_temperature', /CHAR 
      NCDF_ATTPUT, id, tmax_id, 'units', tunits, /CHAR 
      NCDF_ATTPUT, id, tmax_id, '_FillValue', tfill 
      NCDF_ATTPUT, id, tmax_id, 'scale_factor', tscale_factor
      NCDF_ATTPUT, id, tmax_id, 'add_offset', tadd_offset
      NCDF_ATTPUT, id, tmax_id, 'coordinates', coordinates, /CHAR

      ; Tmin random uncertainty
      unc_rand_tmin_id=NCDF_VARDEF(id, 'unc_rand_tasmin', [aid, xid,yid], /SHORT, GZIP=9)
      NCDF_ATTPUT, id, unc_rand_tmin_id, 'long_name', 'Random uncertainty on minimum daily surface air temperature (K)', /CHAR
      ;NCDF_ATTPUT, id, unc_rand_tmin_id, 'standard_name', 'uncertainty', /CHAR 
      NCDF_ATTPUT, id, unc_rand_tmin_id, 'units', tunits, /CHAR 
      NCDF_ATTPUT, id, unc_rand_tmin_id, '_FillValue', ufill
      NCDF_ATTPUT, id, unc_rand_tmin_id, 'scale_factor', uscale_factor
      NCDF_ATTPUT, id, unc_rand_tmin_id, 'add_offset', uadd_offset
      NCDF_ATTPUT, id, unc_rand_tmin_id, 'coordinates', coordinates, /CHAR

      ; Tmax random uncertainty
      unc_rand_tmax_id=NCDF_VARDEF(id, 'unc_rand_tasmax', [aid, xid,yid], /SHORT, GZIP=9)
      NCDF_ATTPUT, id, unc_rand_tmax_id, 'long_name', 'Random uncertainty on maximum daily surface air temperature (K)', /CHAR
      ;NCDF_ATTPUT, id, unc_rand_tmax_id, 'standard_name', 'uncertainty', /CHAR 
      NCDF_ATTPUT, id, unc_rand_tmax_id, 'units', tunits, /CHAR 
      NCDF_ATTPUT, id, unc_rand_tmax_id, '_FillValue', ufill 
      NCDF_ATTPUT, id, unc_rand_tmax_id, 'scale_factor', uscale_factor
      NCDF_ATTPUT, id, unc_rand_tmax_id, 'add_offset', uadd_offset
      NCDF_ATTPUT, id, unc_rand_tmax_id, 'coordinates', coordinates, /CHAR

      ; Tmin locally correlated uncertainty - ATM
      unc_corr_atm_tmin_id=NCDF_VARDEF(id, 'unc_corr_atm_tasmin', [aid, xid,yid], /SHORT, GZIP=9)
      NCDF_ATTPUT, id, unc_corr_atm_tmin_id, 'long_name', 'Locally correlated uncertainty (atmospheric scale lengths) on minimum daily surface air temperature (K)', /CHAR
      ;NCDF_ATTPUT, id, unc_corr_tmin_id, 'standard_name', 'uncertainty', /CHAR 
      NCDF_ATTPUT, id, unc_corr_atm_tmin_id, 'units', tunits, /CHAR
      NCDF_ATTPUT, id, unc_corr_atm_tmin_id, '_FillValue', ufill
      NCDF_ATTPUT, id, unc_corr_atm_tmin_id, 'scale_factor', uscale_factor
      NCDF_ATTPUT, id, unc_corr_atm_tmin_id, 'add_offset', uadd_offset
      NCDF_ATTPUT, id, unc_corr_atm_tmin_id, 'coordinates', coordinates, /CHAR

      ; Tmax locally correlated uncertainty - ATM
      unc_corr_atm_tmax_id=NCDF_VARDEF(id, 'unc_corr_atm_tasmax', [aid, xid,yid], /SHORT, GZIP=9)
      NCDF_ATTPUT, id, unc_corr_atm_tmax_id, 'long_name', 'Locally correlated uncertainty (atmospheric scale lengths) on maximum daily surface air temperature (K)', /CHAR
      ;NCDF_ATTPUT, id, unc_corr_tmax_id, 'standard_name', 'uncertainty', /CHAR 
      NCDF_ATTPUT, id, unc_corr_atm_tmax_id, 'units', tunits, /CHAR
      NCDF_ATTPUT, id, unc_corr_atm_tmax_id, '_FillValue', ufill 
      NCDF_ATTPUT, id, unc_corr_atm_tmax_id, 'scale_factor', uscale_factor
      NCDF_ATTPUT, id, unc_corr_atm_tmax_id, 'add_offset', uadd_offset
      NCDF_ATTPUT, id, unc_corr_atm_tmax_id, 'coordinates', coordinates, /CHAR

      ; Tmin locally correlated uncertainty - SFC
      unc_corr_sfc_tmin_id=NCDF_VARDEF(id, 'unc_corr_sfc_tasmin', [aid, xid,yid], /SHORT, GZIP=9)
      NCDF_ATTPUT, id, unc_corr_sfc_tmin_id, 'long_name', 'Locally correlated uncertainty (surface scale lengths) on minimum daily surface air temperature (K)', /CHAR
      ;NCDF_ATTPUT, id, unc_corr_tmin_id, 'standard_name', 'uncertainty', /CHAR 
      NCDF_ATTPUT, id, unc_corr_sfc_tmin_id, 'units', tunits, /CHAR
      NCDF_ATTPUT, id, unc_corr_sfc_tmin_id, '_FillValue', ufill
      NCDF_ATTPUT, id, unc_corr_sfc_tmin_id, 'scale_factor', uscale_factor
      NCDF_ATTPUT, id, unc_corr_sfc_tmin_id, 'add_offset', uadd_offset
      NCDF_ATTPUT, id, unc_corr_sfc_tmin_id, 'coordinates', coordinates, /CHAR

      ; Tmax locally correlated uncertainty - SFC
      unc_corr_sfc_tmax_id=NCDF_VARDEF(id, 'unc_corr_sfc_tasmax', [aid, xid,yid], /SHORT, GZIP=9)
      NCDF_ATTPUT, id, unc_corr_sfc_tmax_id, 'long_name', 'Locally correlated uncertainty (surface scale lengths) on maximum daily surface air temperature (K)', /CHAR
      ;NCDF_ATTPUT, id, unc_corr_tmax_id, 'standard_name', 'uncertainty', /CHAR 
      NCDF_ATTPUT, id, unc_corr_sfc_tmax_id, 'units', tunits, /CHAR
      NCDF_ATTPUT, id, unc_corr_sfc_tmax_id, '_FillValue', ufill 
      NCDF_ATTPUT, id, unc_corr_sfc_tmax_id, 'scale_factor', uscale_factor
      NCDF_ATTPUT, id, unc_corr_sfc_tmax_id, 'add_offset', uadd_offset
      NCDF_ATTPUT, id, unc_corr_sfc_tmax_id, 'coordinates', coordinates, /CHAR

      ; Tmin systematic uncertainty
      unc_sys_tmin_id=NCDF_VARDEF(id, 'unc_sys_tasmin', [aid, xid,yid], /SHORT, GZIP=9)
      NCDF_ATTPUT, id, unc_sys_tmin_id, 'long_name', 'Systematic uncertainty on minimum daily surface air temperature (K)', /CHAR
      ;NCDF_ATTPUT, id, unc_sys_tmin_id, 'standard_name', 'uncertainty', /CHAR
      NCDF_ATTPUT, id, unc_sys_tmin_id, 'units', tunits, /CHAR
      NCDF_ATTPUT, id, unc_sys_tmin_id, '_FillValue', ufill
      NCDF_ATTPUT, id, unc_sys_tmin_id, 'scale_factor', uscale_factor
      NCDF_ATTPUT, id, unc_sys_tmin_id, 'add_offset', uadd_offset
      NCDF_ATTPUT, id, unc_sys_tmin_id, 'coordinates', coordinates, /CHAR

      ; Tmax systematic uncertainty
      unc_sys_tmax_id=NCDF_VARDEF(id, 'unc_sys_tasmax', [aid, xid,yid], /SHORT, GZIP=9)
      NCDF_ATTPUT, id, unc_sys_tmax_id, 'long_name', 'Systematic uncertainty on maximum daily surface air temperature (K)', /CHAR
      ;NCDF_ATTPUT, id, unc_sys_tmax_id, 'standard_name', 'uncertainty', /CHAR 
      NCDF_ATTPUT, id, unc_sys_tmax_id, 'units', tunits, /CHAR 
      NCDF_ATTPUT, id, unc_sys_tmax_id, '_FillValue', ufill 
      NCDF_ATTPUT, id, unc_sys_tmax_id, 'scale_factor', uscale_factor
      NCDF_ATTPUT, id, unc_sys_tmax_id, 'add_offset', uadd_offset
      NCDF_ATTPUT, id, unc_sys_tmax_id, 'coordinates', coordinates, /CHAR

      ; Tmin model number
      mod_num_tmin_id=NCDF_VARDEF(id, 'tasmin_model_number', [aid, xid,yid], /SHORT, GZIP=9)
      NCDF_ATTPUT, id, mod_num_tmin_id, 'long_name', 'Model number used for estimating Tmin from satellite data', /CHAR
      NCDF_ATTPUT, id, mod_num_tmin_id, '_FillValue', 0
      NCDF_ATTPUT, id, mod_num_tmin_id, 'scale_factor', 1
      NCDF_ATTPUT, id, mod_num_tmin_id, 'add_offset', 0
      NCDF_ATTPUT, id, mod_num_tmin_id, 'coordinates', coordinates, /CHAR

      ; Tmax model number
      mod_num_tmax_id=NCDF_VARDEF(id, 'tasmax_model_number', [aid, xid,yid], /SHORT, GZIP=9)
      NCDF_ATTPUT, id, mod_num_tmax_id, 'long_name', 'Model number used for estimating Tmax from satellite data', /CHAR
      NCDF_ATTPUT, id, mod_num_tmax_id, '_FillValue', 0
      NCDF_ATTPUT, id, mod_num_tmax_id, 'scale_factor', 1
      NCDF_ATTPUT, id, mod_num_tmax_id, 'add_offset', 0
      NCDF_ATTPUT, id, mod_num_tmax_id, 'coordinates', coordinates, /CHAR

      ; Put file in data mode: 
      NCDF_CONTROL, id, /ENDEF 

      ; Scale data for file-writing
      invalid=WHERE(FINITE(tmin) EQ 0, ninvalid)
      tmin=ROUND((tmin-tadd_offset)/tscale_factor)
      IF ninvalid GT 0 THEN $
	  tmin[invalid]=tfill
      tmin_sys=ROUND((tmin_sys-uadd_offset)/uscale_factor)
      IF ninvalid GT 0 THEN $
	  tmin_sys[invalid]=tfill
      invalid=WHERE(FINITE(tmin_ran) EQ 0, ninvalid)
      tmin_ran=ROUND((tmin_ran-uadd_offset)/uscale_factor)
      IF ninvalid GT 0 THEN $
	  tmin_ran[invalid]=ufill
      invalid=WHERE(FINITE(tmin_loc_atm) EQ 0, ninvalid)
      tmin_loc_atm=ROUND((tmin_loc_atm-uadd_offset)/uscale_factor)
      IF ninvalid GT 0 THEN $
	  tmin_loc_atm[invalid]=tfill
      invalid=WHERE(FINITE(tmin_loc_sfc) EQ 0, ninvalid)
      tmin_loc_sfc=ROUND((tmin_loc_sfc-uadd_offset)/uscale_factor)
      IF ninvalid GT 0 THEN $
	  tmin_loc_sfc[invalid]=tfill
	  
      invalid=WHERE(FINITE(tmax) EQ 0, ninvalid)
      tmax=ROUND((tmax-tadd_offset)/tscale_factor)
      IF ninvalid GT 0 THEN $
	  tmax[invalid]=tfill
      tmax_sys=ROUND((tmax_sys-uadd_offset)/uscale_factor)
      IF ninvalid GT 0 THEN $
	  tmax_sys[invalid]=tfill
      invalid=WHERE(FINITE(tmax_ran) EQ 0, ninvalid)
      tmax_ran=ROUND((tmax_ran-uadd_offset)/uscale_factor)
      IF ninvalid GT 0 THEN $
	  tmax_ran[invalid]=ufill
      invalid=WHERE(FINITE(tmax_loc_atm) EQ 0, ninvalid)
      tmax_loc_atm=ROUND((tmax_loc_atm-uadd_offset)/uscale_factor)
      IF ninvalid GT 0 THEN $
	  tmax_loc_atm[invalid]=tfill
      invalid=WHERE(FINITE(tmax_loc_sfc) EQ 0, ninvalid)
      tmax_loc_sfc=ROUND((tmax_loc_sfc-uadd_offset)/uscale_factor)
      IF ninvalid GT 0 THEN $
	  tmax_loc_sfc[invalid]=tfill

      IF MIN(tmin) LT tfill OR $
	  MIN(tmax) LT tfill OR $
	  MIN(tmin_ran) LT tfill OR $
	  MIN(tmin_loc_atm) LT tfill OR $
	  MIN(tmin_loc_sfc) LT tfill OR $
	  MIN(tmax_ran) LT tfill OR $
	  MIN(tmax_loc_atm) LT tfill OR $
	  MIN(tmax_loc_sfc) LT tfill THEN $
	  MESSAGE, 'Warning: Integer wrap-around detected for date: '+yyyymmdd

      ; Input data: 
      NCDF_VARPUT, id, timeid, LONG(i-JULDAY(1, 1, 1850))
      NCDF_VARPUT, id, latid, REFORM(lat[0, *]) 
      NCDF_VARPUT, id, lonid, REFORM(lon[*, 0])
      NCDF_VARPUT, id, tmin_id, REFORM(FIX(tmin), 1, nx, ny)
      NCDF_VARPUT, id, tmax_id, REFORM(FIX(tmax), 1, nx, ny)
      NCDF_VARPUT, id, unc_rand_tmin_id, REFORM(FIX(tmin_ran), 1, nx, ny)
      NCDF_VARPUT, id, unc_rand_tmax_id, REFORM(FIX(tmax_ran), 1, nx, ny)
      NCDF_VARPUT, id, unc_corr_atm_tmin_id, REFORM(FIX(tmin_loc_atm), 1, nx, ny)
      NCDF_VARPUT, id, unc_corr_atm_tmax_id, REFORM(FIX(tmax_loc_atm), 1, nx, ny)
      NCDF_VARPUT, id, unc_corr_sfc_tmin_id, REFORM(FIX(tmin_loc_sfc), 1, nx, ny)
      NCDF_VARPUT, id, unc_corr_sfc_tmax_id, REFORM(FIX(tmax_loc_sfc), 1, nx, ny)
      NCDF_VARPUT, id, unc_sys_tmin_id, REFORM(FIX(tmin_sys), 1, nx, ny)
      NCDF_VARPUT, id, unc_sys_tmax_id, REFORM(FIX(tmax_sys), 1, nx, ny)
      NCDF_VARPUT, id, mod_num_tmin_id, REFORM(tmin_model_num, 1, nx, ny)
      NCDF_VARPUT, id, mod_num_tmax_id, REFORM(tmax_model_num, 1, nx, ny)
	    
      ; Close the NetCDF file & zip
      NCDF_CLOSE, id 
      
    ENDFOR
endforeach





END
