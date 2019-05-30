;-------------------------------------------------------------------------------------
;                            PROGRAM: model_land_idl_setup.pro
;
;
;-------------------------------------------------------------------------------------
; AUTHOR:  Francesco Capponi
; DATE:    16/10/2017 
;
; This is a collection of idl procedures and functions, to be used during the prediction of LSAT.
; Original authors of procedures and functions are reported

pro match, a, b, suba, subb, COUNT = count, SORT = sort, $
           COMPLEMENTA = complementa, COMPLEMENTB = complementb 
; 
; downloaded by Lizzie Good on 05/06/08 from: 
; http://astro.uni-tuebingen.de/software/idl/astrolib/misc/match.pro
;+
; NAME:
;       MATCH
; PURPOSE:
;       Routine to match values in two vectors.
;
; CALLING SEQUENCE:
;       match, a, b, suba, subb, [ COUNT =, /SORT ]
;
; INPUTS:
;       a,b - two vectors to match elements, numeric or string data types
;
; OUTPUTS:
;       suba - subscripts of elements in vector a with a match
;               in vector b
;       subb - subscripts of the positions of the elements in
;               vector b with matchs in vector a.
;
;       suba and subb are ordered such that a[suba] equals b[subb]
;
; OPTIONAL INPUT KEYWORD:
;       /SORT - By default, MATCH uses two different algorithm: (1) the 
;               /REVERSE_INDICES keyword to HISTOGRAM is used for integer data,
;               while a sorting algorithm is used for non-integer data.   The
;               histogram algorithm is usually faster, except when the input
;               vectors are sparse and contain very large numbers, possibly
;               causing memory problems.   Use the /SORT keyword to always use
;               the sort algorithm.
;               
; OPTIONAL KEYWORD OUTPUT:
;       COUNT - set to the number of matches, integer scalar
;       COMPLEMENTA - set to the complement of a where there are no matches
;       COMPLEMENTB - set to the complement of b where there are no matches
;
; SIDE EFFECTS:
;       The obsolete system variable !ERR is set to the number of matches;
;       however, the use !ERR is deprecated in favor of the COUNT keyword 
;
; RESTRICTIONS:
;       The vectors a and b should not have duplicate values within them.
;       You can use rem_dup function to remove duplicate values
;       in a vector
;
; EXAMPLE:
;       If a = [3,5,7,9,11]   & b = [5,6,7,8,9,10]
;       then 
;               IDL> match, a, b, suba, subb, COUNT = count
;
;       will give suba = [1,2,3], subb = [0,2,4],  COUNT = 3
;       and       suba[a] = subb[b] = [5,7,9]
;
; 
; METHOD:
;       For non-integer data types, the two input vectors are combined and
;       sorted and the consecutive equal elements are identified.   For integer
;       data types, the /REVERSE_INDICES keyword to HISTOGRAM of each array
;       is used to identify where the two arrays have elements in common.   
; HISTORY:
;       D. Lindler  Mar. 1986.
;       Fixed "indgen" call for very large arrays   W. Landsman  Sep 1991
;       Added COUNT keyword    W. Landsman   Sep. 1992
;       Fixed case where single element array supplied   W. Landsman Aug 95
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Use a HISTOGRAM algorithm for integer vector inputs for improved 
;             performance                W. Landsman         March 2000
;       Work again for strings           W. Landsman         April 2000
;       Use size(/type)                  W. Landsman         December 2002
;       Work for scalar integer input    W. Landsman         June 2003
;       Added COMPLEMENTA and COMPLEMENTB keywords - Elizabeth Good, Jun 2008
;-
;-------------------------------------------------------------------------
 On_error,2

 if N_params() LT 3 then begin
     print,'Syntax - match, a, b, suba, subb, [ COUNT = ]'
     print,'    a,b -- input vectors for which to match elements'
     print,'    suba,subb -- output subscript vectors of matched elements'
     return
 endif

 da = size(a,/type) & db =size(b,/type)
 if keyword_set(sort) then hist = 0b else $
 hist = (( da LE 3 ) or (da GE 12)) and  ((db LE 3) or (db GE 12 )) 

 na = N_elements(a)              ;number of elements in a
 nb = N_elements(b)             ;number of elements in b

 if not hist then begin           ;Non-integer calculation
 

; Check for a single element array

 if (na EQ 1) or (nb EQ 1) then begin
        if (nb GT 1) then begin
                subb = where(b EQ a[0], nw, COMPLEMENT=complementb)
                if (nw GT 0) then BEGIN
		   suba = replicate(0,nw)
		   complementa= REPLICATE(1,nw)
		ENDIF else BEGIN
		   suba = [-1]
		   complementa = [0]
		ENDELSE
        endif else begin
                suba = where(a EQ b[0], nw, COMPLEMENT=complementa)
                if (nw GT 0) then BEGIN
		   subb = replicate(0,nw) 
		   complementb = REPLICATE(1, nw)
		ENDIF else BEGIN
		   subb = [-1]
		   complementb = [0]
		ENDELSE
        endelse
        count = nw
        return
 endif
        
 c = [ a, b ]                   ;combined list of a and b
 ind = [ lindgen(na), lindgen(nb) ]       ;combined list of indices
 vec = [ bytarr(na), replicate(1b,nb) ]  ;flag of which vector in  combined 
                                         ;list   0 - a   1 - b

; sort combined list

 sub = sort(c)
 c = c[sub]
 ind = ind[sub]
 vec = vec[sub]

; find duplicates in sorted combined list

 n = na + nb                            ;total elements in c
 firstdup = where( (c EQ shift(c,-1)) and (vec NE shift(vec,-1)), Count)

 if Count EQ 0 then begin               ;any found?
        suba = lonarr(1)-1
        subb = lonarr(1)-1
	complementa = LINDGEN(na)
	complementb = LINDGEN(nb)
        return
 end
 
 dup = lonarr( Count*2 )                     ;both duplicate values
 even = lindgen( N_elements(firstdup))*2     ;Changed to LINDGEN 6-Sep-1991
 dup[even] = firstdup
 dup[even+1] = firstdup+1
 ind = ind[dup]                         ;indices of duplicates
 vec = vec[dup]                         ;vector id of duplicates
 suba = ind[ where( vec EQ 0)  ]       ;a subscripts
 subb = ind[ where( vec) ]             ;b subscripts

 endif else begin             ;Integer calculation using histogram.

 minab = min(a, MAX=maxa) > min(b, MAX=maxb) ;Only need intersection of ranges
 maxab = maxa < maxb

;If either set is empty, or their ranges don't intersect: 
;  result = NULL (which is denoted by integer = -1)
  !ERR = -1
  suba = -1
  subb = -1
  COUNT = 0L
  complementa = LINDGEN(na)
  complementb = LINDGEN(nb)
 if (maxab lt minab) or (maxab lt 0) then return
 
 ha = histogram([a], MIN=minab, MAX=maxab, reverse_indices=reva)
 hb = histogram([b], MIN=minab, MAX=maxab, reverse_indices=revb)
 
 r = where((ha ne 0) and (hb ne 0), count)
 if count gt 0 then begin
  suba = reva[reva[r]]
  subb = revb[revb[r]]
 endif 
 endelse 
 
 ; Establish complementary arrays of indices for a and b - i.e. where the 
 ; arrays do not match.
 arra=LINDGEN(na)
 IF suba[0] NE -1 THEN arra[suba]=-1
 complementa=WHERE(arra NE -1)
 arrb=LINDGEN(nb)
 IF subb[0] NE -1 THEN arrb[subb]=-1
 complementb=WHERE(arrb NE -1)
 
 return
 
 end


PRO read_asc, file, data, data_names, header, head_str=head_str, $
        strings=strings, double=double

;-------------------------------------------------------------------------------------
;                            PROGRAM: PRO read_asc.pro
;
;
;-------------------------------------------------------------------------------------
; AUTHOR:  Elizabeth Noyes
; DATE:    May 2004
;
; ABOUT THE PROGRAM:
;
; Procedure to read in an ascii file.  Ascii file may have n columns of data,
; but must have one line of column labels immediately before the data.  A header
; of any number of lines is permitted, but the first character of the line must
; be present, denoting the header lines.  E.g.:
;
; # This is an example header
; # It can contain any number of lines
; This   line    must   contain    data    labels
; data1  data2   data3  data4      data5   data6
;
; ....any number of rows of data/data records will work.
;
;
; LIST OF ARGUMENTS:
;
; INPUT:
; file - full pathname to data file
;
; OUTPUT:
; data - 2D floating point array containing data
; data_names - string array containing corresponding column labels
; header - string array containing header lines
;
; KEYWORDS
; head_str - can be set to indicate the first character of header lines.
;            Default is '!'
; strings - set to indicate that the columns of data should not be converted
; to floats.  Use this option if one or more columns of the data contain strings.
; double - set to indicate that the columns of data should be read in as doubles.
;
;
; LIMITATIONS:
;
; must contain n columns of data with one containing the column data labels.
;
; i.e.
;
; col1_name   col2_name   col3_name
; 1data_1      1data_2      1data_3
; 2data_1      2data_2      2data_3
; ...
; ndata_1      ndata_2      ndata_3
;
;
; MODIFICATIONS (version 2.0)
;
; BY:      Elizabeth Noyes
; DATE:    October 2004
;
; Elimation of variable ncols.
;
; MODIFICATIONS (version 3.0)
;
; BY:      Elizabeth Noyes
; DATE:    October 2004
;
; Inclusion of keyword to return data in strings and not numbers, as version
; 2 will only work where the data is numbers.  Also switch to return data as
; doubles.
;
; MODIFICATIONS (version 4.0) - E Noyes, June 2006
;
; Now checks every line of code read in to verify there are values in the data file. If not
; the program stops reading in the file.
;
; MODIFICATIONS (version 5.0) - E Noyes Feb 2008
; Modified to return -1 or '-1' if no data found in file in data variable.
;
;
;
;-------------------------------------------------------------------------------------

; Check keywords.  Default is 3 cols of data and header lines denoted by '!'

IF NOT KEYWORD_SET(head_str) THEN head_str='!'



OPENR, lun, file, /GET_LUN    ; Open file



; Define loop variables.

head_count=0      ; counts number of header lines
datanames_count=0 ; counts column header lines
data_count=0      ; counts number of data lines

dummy_str=' '     ; dummy string to house each line file
;dummy_data=FLTARR(ncols)    ; dummy floating point array to house data


;data_names=STRARR(ncols)      ; array to house column headers



; Get header information using the defined header string.  Do this by reading in
; one line of the file at a time while not end of file.

WHILE NOT EOF(lun) DO BEGIN

   READF, lun, dummy_str

   IF (STRCOMPRESS(dummy_str, /REMOVE_ALL) EQ '') THEN BREAK


   ; If first character of dummy_str is the character that marks a header line,
   ; then read this into the variable 'header'.  This set up also allows reading
   ; of files that do not contain any header information.

   mark=STRMID(dummy_str, 0, 1)

   IF (mark EQ head_str) THEN BEGIN

      IF (head_count EQ 0) THEN header=dummy_str ELSE header=[header, dummy_str]

      head_count=head_count+1

   ENDIF ELSE BEGIN


      ; If line of data file is not denoted by a header character then read in
      ; one line of data labels followed by n lines of data.  However, as some
      ; files have rows of blank lines at the end of the file, first check we
      ; actually have some values in this line of the file.  Otherwise break from
      ; the code.

      dummy_str2=STRCOMPRESS(dummy_str)    ; compress line of file (necessary)



         IF (datanames_count EQ 0) THEN BEGIN

            data_names=STRSPLIT(dummy_str2, /EXTRACT)  ; extract column labels
            datanames_count=1    ; increment counter to indicate col labels read


         ENDIF ELSE BEGIN


            IF KEYWORD_SET(strings) THEN dummy_data=STRSPLIT(dummy_str2, /EXTRACT) $
            ELSE IF KEYWORD_SET(double) THEN dummy_data=DOUBLE(STRSPLIT(dummy_str2, /EXTRACT)) $
            ELSE dummy_data=FLOAT(STRSPLIT(dummy_str2, /EXTRACT))

;STOP
         IF (data_count EQ 0) THEN data=dummy_data ELSE data=[data, dummy_data]

            data_count=data_count+1

         ENDELSE



   ENDELSE

ENDWHILE

CLOSE, lun

FREE_LUN, lun


; reform data into ncols.

IF data_count NE 0 THEN $
   data=REFORM(data, N_ELEMENTS(dummy_data), N_ELEMENTS(data)/N_ELEMENTS(dummy_data)) $
ELSE $
   IF KEYWORD_SET(strings) THEN $
      data='-1' $
   ELSE $
      data=-1


END


FUNCTION date_to_doy, yyyymmdd, string=string, start_zero=start_zero, julian=julian

;-----------------------------------------------------------------------------
; BY:  Elizabeth Good
; DATE: 22 March 2010
;
; ABOUT THE PROGRAM
; Quick function to convert a date in YYYYMMDD format (String) to a day of the
; year.
; 
; Example: 
; IDL> PRINT, date_to_doy('20100322')
;    81    
; IDL> PRINT, date_to_doy('20100322', /start_zero, /string)
; 080    
; 
; KEYWORDS:
; string - set this keyword (/string) to return day of year as three-digit
;          string with leading zeros, e.g. 001.
; start_zero - default is to start 01 Jan as day 1, but if you want to start
;              day of year with 0 (i.e. 01 Jan is day 0), then set this keyword.
; julian - set this keyword (/julian) to return Julian day instead of day of year (as 
;          defined by IDL JULDAY function.
;
; MODIFICATIONS:
; 11 July 2011, Elizabeth Good - added Julian keyword.
;-----------------------------------------------------------------------------

;---------------
; CHECK INPUTS
;---------------

IF MAX(STRLEN(yyyymmdd)) NE 8 AND MIN(STRLEN(yyyymmdd)) THEN $
   MESSAGE, 'Input date not recognised; must be in YYYYMMDD format as a string'

IF KEYWORD_SET(start_zero) THEN add_on=0 ELSE add_on=1

;---------------
; CALCULATION
;---------------

; Separate year, month and day
yyyy=STRMID(yyyymmdd, 0, 4)
mm=STRMID(yyyymmdd, 4, 2)
dd=STRMID(yyyymmdd, 6, 2)

IF KEYWORD_SET(julian) THEN $
   RETURN, JULDAY(FIX(mm), FIX(dd), FIX(yyyy)) 
doy=JULDAY(FIX(mm), FIX(dd), FIX(yyyy))-JULDAY(1, 1, FIX(yyyy))+add_on


;------------------------------
; CONVERT TO STRING IF REQUIRED
;------------------------------

IF KEYWORD_SET(string) THEN BEGIN
   doy='00'+STRCOMPRESS(doy, /REMOVE_ALL)
   doy=STRMID(doy, 2, 3, /REVERSE_OFFSET)
ENDIF


;---------------
; PROGRAM END
;---------------

RETURN, doy

END


; This code is sourced from: http://cheas.psu.edu/data/lib/general/adesai/sunrise.pro and was
; obtained by Elizabeth Good (Met Office, April 2008).
; Code verification: checked against NASA tool at 
; http://www.srrb.noaa.gov/highlights/sunrise/sunrise.html
; which gives consistent results.
; E Good, 2009
; Made small change to add on 0.0001 to longitudes if input longitude is 0
; as the code does not work for an exact 0 degrees.
; E Good, August 2009
; Added keyword 'instrumental' which allows definition of sunrise and sunset consistent with
;       first/last sunlight observed by an in situ instrument.
;
;---------------------------------------------------------------------------------------------

;Sunrise, Sunset, Solar noon calculator
;Source: http://www.srrb.noaa.gov/highlights/sunrise/sunrise.html
;Converted to IDL by Ankur Desai, 7 April 2003
;Main function is sunrise (at bottom)

FUNCTION DayofYear,mn,dy,yr,leapyear=leapyear
;Give month and day and year (can be arrays), return day of year
;year is optional, in which case keyword leapyear determines leapyearness

  IF NOT keyword_set(yr) THEN BEGIN
    IF keyword_set(leapyear) THEN lp = 1b ELSE lp = 0b
  ENDIF ELSE BEGIN 
    lp = ((((Yr mod 4) eq 0) and ((yr mod 100) ne 0)) or ((yr mod 400) eq 0))
  ENDELSE
  k = 2 - lp
  doy = floor((275*mn)/9) - k * floor((mn+9)/12) + dy - 30
  return,doy
END

FUNCTION CalcMDY,doy,yr,leapyear=leapyear
;Given the Day of Year and Year find the Month and Day of Month
;Year is optional

  IF NOT keyword_set(yr) THEN BEGIN
    IF keyword_set(leapyear) THEN lp = 1b ELSE lp = 0b
  ENDIF ELSE BEGIN 
    lp = ((((Yr mod 4) eq 0) and ((yr mod 100) ne 0)) or ((yr mod 400) eq 0))
  ENDELSE
  
  md = [31,28,31,30,31,30,31,31,30,31,30,31]
  IF lp THEN md[1] = 29
  modcum = fix(total(md,/cumulative))
  modcum = fix(doy)-(modcum - md)
  l = where(modcum le 0)
  if l[0] ne -1 then modcum[l] = 999
  dy = min(modcum)
  mo = (where(modcum eq dy))[0] + 1
  return,[mo,dy]
END 

FUNCTION calcJD,yr, mo, dy
;Calculate Julian Day (at 0 UTC) given Year, Month, Day
  year=long(yr)
  month = long(mo)
  IF month LE 2l THEN BEGIN
    year = year - 1l
    month = month + 12l
  ENDIF
  a = floor(year/100l)
  b = 2l - a + floor(a/4l)
  year = double(year)
  month = double(month)
  day = double(dy)
  JD = floor(365.25d*(year + 4716.0d)) + floor(30.6001d*(month+1.0d)) + day + B - 1524.5d
  return,jd
END

FUNCTION calcTimeJulianCent,jd
;convert Julian Day to centuries since J2000.0.	
  T = (double(jd) - 2451545.0d)/36525.0d
  return,t
END

FUNCTION calcJDFromJulianCent,t
;convert centuries since J2000.0 to Julian Day.	
  JD = double(t) * 36525.0d + 2451545.0d
  return,jd
END

FUNCTION calcGeomMeanLongSun,t
;calculate the Geometric Mean Longitude of the Sun (in degrees)
  L0 = (280.46646d + double(t) * (36000.76983d + 0.0003032d * double(t)))
  WHILE (l0 GT 360.0d) DO l0 = l0-360.0d
  WHILE (l0 LT 0.0d) DO l0 = l0+360.0d
  return,l0
END 

FUNCTION calcGeomMeanAnomalySun,t
;calculate the Geometric Mean Anomaly of the Sun (in degrees)
  M = 357.52911d + double(t) * (35999.05029d - 0.0001537d * double(t))
  return,m
END 

FUNCTION calcEccentricityEarthOrbit,t
;calculate the eccentricity of earth's orbit (unitless)
  e = 0.016708634d - double(t) * (0.000042037d + 0.0000001267d * double(t)) 
  return, e
END

FUNCTION calcSunEqOfCenter,t
;calculate the equation of center for the sun (in degrees)
  m = double(calcGeomMeanAnomalySun(t))
  mrad = (!dpi/180.0d)*m
  sinm = sin(mrad) 
  sin2m = sin(mrad+mrad) 
  sin3m = sin(mrad+mrad+mrad)
  C = sinm * (1.914602d - double(t) * (0.004817d + 0.000014d * double(t))) + sin2m * (0.019993d - 0.000101d * double(t)) + sin3m * 0.000289d 
  return,c
END

FUNCTION calcSunTrueLong,t
;calculate the true longitude of the sun (in degrees)
  l0 = calcGeomMeanLongSun(double(t))
  c = calcSunEqOfCenter(double(t))
  O = l0 + c               
  return,O
END

FUNCTION calcSunTrueAnomaly,t
;calculate the true anamoly of the sun (in degrees)
  m = calcGeomMeanAnomalySun(double(t))
  c = calcSunEqOfCenter(double(t)) 
  v = m + c
  return,v
END

FUNCTION calcSunRadVector,t
;calculate the distance to the sun in AU (in degrees)
  v = calcSunTrueAnomaly(double(t))
  e = calcEccentricityEarthOrbit(double(t)) 
  R = (1.000001018d * (1.0d - e * e)) / (1.0d + e * cos(((!dpi/180.0d)*v))) 
  return,r
END 

FUNCTION calcSunApparentLong,t
;calculate the apparent longitude of the sun (in degrees)
  o = calcSunTrueLong(double(t))
  omega = 125.04d - 1934.136d * double(t) 
  lambda = o - 0.00569d - 0.00478d * sin(((!dpi/180.0d)*omega)) 
  return,lambda
END 

FUNCTION calcMeanObliquityOfEcliptic,t
;calculate the mean obliquity of the ecliptic (in degrees)
  seconds = 21.448d - double(t)*(46.8150d + double(t)*(0.00059d - double(t)*(0.001813d))) 
  e0 = 23.0d + (26.0d + (seconds/60.0d))/60.0d 
  return,e0
END 

FUNCTION calcObliquityCorrection,t
;calculate the corrected obliquity of the ecliptic (in degrees)
  e0 = calcMeanObliquityOfEcliptic(double(t)) 
  omega = 125.04d - 1934.136d * double(t) 
  e = e0 + 0.00256d * cos(((!dpi/180.0d)*omega)) 
  return,e
END 

FUNCTION calcSunRtAscension,t
;calculate the right ascension of the sun (in degrees)
  e = calcObliquityCorrection(double(t)) 
  lambda = calcSunApparentLong(double(t)) 
  tananum = (cos(((!dpi/180.0d)*e)) * sin(((!dpi/180.0d)*lambda))) 
  tanadenom = (cos(((!dpi/180.0d)*lambda)))
  alpha = (180.0d/!dpi)*(atan(tananum, tanadenom))
  return,alpha
END 

FUNCTION calcSunDeclination,t
;calculate the declination of the sun (in degrees)
  e = calcObliquityCorrection(double(t)) 
  lambda = calcSunApparentLong(double(t)) 
  sint = sin(((!dpi/180.0d)*e)) * sin(((!dpi/180.0d)*lambda)) 
  theta = (180.0d/!dpi)*(asin(sint))
  return,theta
END 

FUNCTION calcEquationOfTime,t
;calculate the difference between true solar time and mean solar time
;(output: equation of time in minutes of time)	
  epsilon = calcObliquityCorrection(double(t)) 
  l0 = calcGeomMeanLongSun(double(t)) 
  e = calcEccentricityEarthOrbit(double(t)) 
  m = calcGeomMeanAnomalySun(double(t)) 
  y = tan(((!dpi/180.0d)*epsilon)/2.0d) 
  y = y * y
  sin2l0 = sin(2.0d * ((!dpi/180.0d)*l0)) 
  sinm   = sin(((!dpi/180.0d)*m)) 
  cos2l0 = cos(2.0d * ((!dpi/180.0d)*l0)) 
  sin4l0 = sin(4.0d * ((!dpi/180.0d)*l0))
  sin2m  = sin(2.0d * ((!dpi/180.0d)*m))
  Etime = y * sin2l0 - 2.0d * e * sinm + 4.0d * e * y * sinm * cos2l0 - 0.5d * y * y * sin4l0 - 1.25d * e * e * sin2m
  return,(180.0d/!dpi)*Etime*4.0d
END 

FUNCTION calcHourAngleSunrise,lat,solarDec,civil=civil,nautical=nautical,astronomical=astronomical,instrumental=instrumental 
;calculate the hour angle of the sun at sunrise for the latitude (in radians)
  angle = 90.0d + (5.0d/6.0d)
  IF keyword_set(civil) THEN angle = 96.0d
  IF keyword_set(nautical) THEN angle = 102.0d
  IF keyword_set(astronomical) THEN angle = 108.0d
  IF keyword_set(instrumental) THEN angle = 90.0d + (5.0d/6.0d)-2.5

  latRad = ((!dpi/180.0d)*double(lat))   
  sdRad  = ((!dpi/180.0d)*double(solarDec))
  HAarg = (cos(((!dpi/180.0d)*angle))/(cos(latRad)*cos(sdRad))-tan(latRad) * tan(sdRad)) 
  HA = (acos(cos(((!dpi/180.0d)*angle))/(cos(latRad)*cos(sdRad))-tan(latRad) * tan(sdRad))) 
  return,ha
END 

FUNCTION calcHourAngleSunset,lat, solarDec,civil=civil,nautical=nautical,astronomical=astronomical,instrumental=instrumental
;calculate the hour angle of the sun at sunset for the latitude (in radians)
  angle = 90.0d + (5.0d/6.0d)
  IF keyword_set(civil) THEN angle = 96.0d
  IF keyword_set(nautical) THEN angle = 102.0d
  IF keyword_set(astronomical) THEN angle = 108.0d
  IF keyword_set(instrumental) THEN angle = 90.0d + (5.0d/6.0d)-2.5

  latRad = ((!dpi/180.0d)*double(lat))   
  sdRad  = ((!dpi/180.0d)*double(solarDec))
  HAarg = (cos(((!dpi/180.0d)*angle))/(cos(latRad)*cos(sdRad))-tan(latRad) * tan(sdRad))
  HA = (acos(cos(((!dpi/180.0d)*angle))/(cos(latRad)*cos(sdRad))-tan(latRad) * tan(sdRad)))
  return,(-1.0d)*ha
END   

FUNCTION calcSunriseUTC,JD, latitude, longitude,civil=civil,nautical=nautical,astronomical=astronomical,instrumental=instrumental
;calculate time of sunrise for the given day at the given location on earth
;(in minute since 0 UTC)
  t = calcTimeJulianCent(double(jd))
  eqTime = calcEquationOfTime(double(t)) 
  solarDec = calcSunDeclination(double(t)) 
  hourAngle = calcHourAngleSunrise(double(latitude), double(solarDec),civil=civil,nautical=nautical,astronomical=astronomical,instrumental=instrumental) 
  delta = double(longitude) - ((180.0d/!dpi)*hourAngle) 
  timeDiff = 4.0d * delta    
  timeUTC = 720.0d + timeDiff - eqTime
  newt = calcTimeJulianCent(calcJDFromJulianCent(double(t)) + timeUTC/1440.0d) 
  eqTime = calcEquationOfTime(double(newt)) 
  solarDec = calcSunDeclination(double(newt)) 
  hourAngle = calcHourAngleSunrise(double(latitude), double(solarDec),civil=civil,nautical=nautical,astronomical=astronomical,instrumental=instrumental) 
  delta = double(longitude) - ((180.0d/!dpi)*hourAngle)
  timeDiff = 4.0d * delta          
  timeUTC = 720.0d + timeDiff - eqTime 
  return,timeutc
END 

FUNCTION calcSunsetUTC,JD, latitude, longitude,civil=civil,nautical=nautical,astronomical=astronomical,instrumental=instrumental
;calculate time of sunset for the given day at the given location on earth
;(in minute since 0 UTC)
  t = calcTimeJulianCent(double(JD))
  eqTime = calcEquationOfTime(double(t)) 
  solarDec = calcSunDeclination(double(t)) 
  hourAngle = calcHourAngleSunset(double(latitude), double(solarDec),civil=civil,nautical=nautical,astronomical=astronomical,instrumental=instrumental) 
  delta = double(longitude) - ((180.0d/!dpi)*hourAngle) 
  timeDiff = 4.0d * delta     
  timeUTC = 720.0d + timeDiff - eqTime 
  newt = calcTimeJulianCent(calcJDFromJulianCent(double(t)) + timeUTC/1440.0d) 
  eqTime = calcEquationOfTime(double(newt)) 
  solarDec = calcSunDeclination(double(newt)) 
  hourAngle = calcHourAngleSunset(double(latitude), double(solarDec),civil=civil,nautical=nautical,astronomical=astronomical,instrumental=instrumental)
  delta = double(longitude) - ((180.0d/!dpi)*hourAngle) 
  timeDiff = 4.0d * delta         
  timeUTC = 720.0d + timeDiff - eqTime 
  return,timeutc

END   

FUNCTION calcSolNoonUTC,jd, longitude
;calculate time of solar noon the given day at the given location on earth
;(in minute since 0 UTC)
  t = calcTimeJulianCent(double(JD))
  newt = calcTimeJulianCent(calcJDFromJulianCent(double(t)) + 0.5d + double(longitude)/360.0d) 
  eqTime = calcEquationOfTime(double(newt)) 
  solarNoonDec = calcSunDeclination(double(newt)) 
  solNoonUTC = 720.0d + (double(longitude) * 4.0d) - eqTime 
  return,solnoonutc
END 

FUNCTION sunrise,doy,yr,lat=lat,lon=inlon,melarr=melarr,civil=civil,nautical=nautical,astronomical=astronomical, $
                 instrumental=instrumental
;Main sunrise/sunset calculator by Ankur Desai April 2003

;DOY is Day of Year (1 = 1/1/2003 to 365 or 366)
;  Default is Today (from computer clock)
;Yr is the Year - default is this year

;DOY and/or YR can be arrays (one or both)

;Lat and Lon are latitude and longitude (default is Sylvania Flux Tower)
;melarr is just to keep some comptibility with old code
;Do not give a latitude above the (ant)arctic circle in mid winter or summer
;Negative Longitide is west
;Positive latitude is north

;Default time is Official sunset and sunrise (zenith angle of 90.833)
;Set keyword Civil for civil sunset/rise (96 degrees)
;            Nautical is 102 degrees
;            Astronomical is 108 degrees
;            Instrumental is 90.833-2.5

;Output
;returns time (in hours UTC) of [sunrise,sunset,local noon,length of day]
;If DOY and/or YR is an array, then you get a bunch of these in rows

;Example:
;print,sunrise(180,2002,lat=40.0,lon=-90.0)
;   10.564731   25.548302   18.057778   14.983571
; sunrise (UTC)   sunset      noon      day length (hrs) 

lon=inlon ; preserve user input.
IF lon EQ 0. THEN lon=0.0001  ; Quick fudge to get code working for 0 longitude
                              ; (previously gave incorrect result).
IF doy LE 0 OR doy GT 366 THEN MESSAGE, 'doy out of range'			       

  caldat,systime(/julian),mo,day,year
  IF NOT keyword_Set(doy) THEN BEGIN 
    yr = year
    doy = DayofYear(mo,day,yr)
  ENDIF ELSE IF NOT keyword_set(yr) THEN yr = year
    
  IF NOT keyword_set(lat) THEN lat = 46.242017d
  IF NOT keyword_set(lon) THEN lon = -89.347650d

;calculate month and day from Day of Year and then find Julian Day

  nelem = max([n_elements(doy),n_elements(yr)])

  output = dblarr(4,nelem)
  FOR i = 0,nelem-1 DO BEGIN 
    IF n_elements(doy) GT 1 THEN thedoy = doy[i] ELSE thedoy = doy
    IF n_elements(yr) GT 1 THEN theyr = yr[i] ELSE theyr = yr

    md = CalcMDY(thedoy,theyr)
    mo = md[0]
    dy = md[1]    
    jd = calcjd(theyr,mo,dy)

    rise = CalcSunriseUTC(jd,lat,lon*(-1.0d),civil=civil,nautical=nautical,astronomical=astronomical,instrumental=instrumental)/60.0d
    set = CalcSunSetUTC(jd,lat,lon*(-1.0d),civil=civil,nautical=nautical,astronomical=astronomical,instrumental=instrumental)/60.0d
    noon = CalcSolNoonUTC(jd,lon*(-1.0d))/60.0d
    daylen = set-rise
    output[*,i] = [rise,set,noon,daylen]
  ENDFOR

  return,reform(output)
END


;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    SUN2000 computes the sun vector in geocentric inertial (equatorial)
;    coordinates.
;
;  This subroutine computes the Sun vector in geocentric inertial 
;  (equatorial) coodinates.  It uses the model referenced in The 
;  Astronomical Almanac for 1984, Section S (Supplement) and documented
;  for the SeaWiFS Project in "Constants and Parameters for SeaWiFS
;  Mission Operations", in TBD.  The accuracy of the Sun vector is
;  approximately 0.1 arcminute.
;
pro sun2000,iyr,iday,sec,sun,rs
;	implicit real*8 (a-h,o-z)
;	real*4 sun(3),rs
	xk=0.0056932		;Constant of aberration 
	imon=1
;	common nutcm,dpsi,eps,nutime
;	radeg = 180.d0/3.14159265359d0

;   Compute floating point days since Jan 1.5, 2000 
;    Note that the Julian day starts at noon on the specified date
	jd,iyr,imon,iday,jday
	t = jday - 2451545.0d0 + (sec-43200.d0)/86400.d0
	n = n_elements(t)
	sun=fltarr(3,n)

;  Compute solar ephemeris parameters
	ephparms,t,xls,gs,xlm,omega

;  Check if need to compute nutation corrections for this day
;	nt = t
;	if (nt ne nutime) then begin
;	  nutime = nt
	  nutate,t,xls,gs,xlm,omega,dpsi,eps
;	endif

;  Compute planet mean anomalies
;   Venus Mean Anomaly 	
	g2 = 50.40828 + 1.60213022*t 	
	g2 = g2 mod 360.d0

;   Mars Mean Anomaly 		
	g4 = 19.38816 + 0.52402078*t 	
	g4 = g4 mod 360.d0

;  Jupiter Mean Anomaly 
	g5 = 20.35116 + 0.08309121*t 	
	g5 = g5 mod 360.d0

;  Compute solar distance (AU)
	rs = 1.00014-0.01671*cos(gs/!radeg)-0.00014*cos(2.*gs/!radeg)

;  Compute Geometric Solar Longitude 
	dls = 	(6893. - 4.6543463D-4*t)*sin(gs/!radeg) $
     		+ 72.*sin(2.*gs/!radeg) 		$
     		- 7.*cos((gs - g5)/!radeg) 		$
     		+ 6.*sin((xlm - xls)/!radeg) 		$ 
     		+ 5.*sin((4.*gs - 8.*g4 + 3.*g5)/!radeg) $
     		- 5.*cos((2.*gs - 2.*g2)/!radeg) 	$
     		- 4.*sin((gs - g2)/!radeg) 		$
     		+ 4.*cos((4.*gs - 8.*g4 + 3.*g5)/!radeg) $
     		+ 3.*sin((2.*gs - 2.*g2)/!radeg) 	$
     		- 3.*sin(g5/!radeg) 			$
     		- 3.*sin((2.*gs - 2.*g5)/!radeg)

	xlsg = xls + dls/3600.d0

;  Compute Apparent Solar Longitude; includes corrections for nutation 
;   in longitude and velocity aberration
	xlsa = xlsg + dpsi - xk/rs

;   Compute unit Sun vector 
	sun(0,*) = cos(xlsa/!radeg)
	sun(1,*) = sin(xlsa/!radeg)*cos(eps/!radeg)
	sun(2,*) = sin(xlsa/!radeg)*sin(eps/!radeg)

	return
	end


;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    L_SUN computes unit sun vector in geocentric rotating coordinates.
;
;  Computes unit Sun vector in geocentric rotating coodinates, using 
;  subprograms to compute inertial Sun vector and Greenwich hour angle
pro l_sun,iyr,iday,sec,sunr,rs

;	implicit real*8 (a-h,o-z)
;	real*4 sunr(3),su(3),rs
	
;	radeg = 180.d0/3.14159265359d0

;  Get unit Sun vector in geocentric inertial coordinates
	sun2000,iyr,iday,sec,su,rs

;  Get Greenwich mean sideral angle
	day = iday + sec/864.d2 
	gha2000,iyr,day,gha
	ghar = gha/!radeg

;  Transform Sun vector into geocentric rotating frame
	n = n_elements(day)
	sunr = fltarr(3,n)
	sunr(0,*) = su(0,*)*cos(ghar) + su(1,*)*sin(ghar)
	sunr(1,*) = su(1,*)*cos(ghar) - su(0,*)*sin(ghar)
	sunr(2,*) = su(2,*)

	return
end


;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    JD converts a calendar date to a julian day from noon on the calendar day.
;
;    This function converts a calendar date to the corresponding Julian
;    day starting at noon on the calendar date.  The algorithm used is
;    from Van Flandern and Pulkkinen, Ap. J. Supplement Series 41, 
;    November 1979, p. 400.
;
;
;	Arguments
;     
;     	Name    Type 	I/O 	Description
;     	----	---- 	--- 	-----------
;     	i	I*4  	 I 	Year - e.g. 1970
;     	j       I*4  	 I  	Month - (1-12)
;     	k       I*4  	 I  	Day  - (1-31)
;     	jday    I*4  	 O  	Julian day
;
;     external references
;     -------------------
;      none
;
;
;     Written by Frederick S. Patt, GSC, November 4, 1992
;
;
pro jd,i,j,k,jday
	i = long(i)
	j = long(j)
	k = long(k)
      	jday = 367*i - 7*(i+(j+9)/12)/4 + 275*j/9 + k + 1721014

;  This additional calculation is needed only for dates outside of the 
;   period March 1, 1900 to February 28, 2100
;     	jday = jday + 15 - 3*((i+(j-9)/7)/100+1)/4
      	return
      	end


;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    JDAY returns the julian day.
;
pro jday,i,k,jd
j = 1
l = long((j-14)/12)
jd = k - long(32075) + long(1461)*long(i+4800+l)/4 
jd = jd + long(367)*long(j-2-l*12)/12 - 3*((i+4900+l)/100)/4
return
end


;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    GHA2000 computes the greenwich hour angle in degrees for the input time.
;
;  This subroutine computes the Greenwich hour angle in degrees for the
;  input time.  It uses the model referenced in The Astronomical Almanac
;  for 1984, Section S (Supplement) and documented for the SeaWiFS
;  Project in "Constants and Parameters for SeaWiFS Mission Operations",
;  in TBD.  It includes the correction to mean sideral time for nutation
;  as well as precession.
;
pro gha2000,iyr,day,gha
;	implicit real*8 (a-h,o-z)
	imon=1
;	common nutcm,dpsi,eps,nutime
;	radeg = 180.d0/3.14159265359d0

;  Compute days since J2000
	iday = long(day)
	fday = day - iday
	jday,iyr,iday,jd
	t = jd - 2451545.5d0 + fday
	
;  Compute Greenwich Mean Sidereal Time	(degrees)
	gmst = 100.4606184d0 + 0.9856473663d0*t + 2.908d-13*t*t

;  Check if need to compute nutation correction for this day
;	nt = t
;	if (nt ne nutime) then begin
;	  nutime = nt
	  ephparms,t,xls,gs,xlm,omega
	  nutate,t,xls,gs,xlm,omega,dpsi,eps
;	endif

;  Include apparent time correction and time-of-day
	gha = gmst + dpsi*cos(eps/!radeg) + fday*360.d0
	gha = gha mod 360.d0
	neg = where(gha lt 0.d0)
        if (neg(0) gt -1) then gha(neg) = gha(neg) + 360.d0
;
	return
	end

	pro ephparms,t,xls,gs,xlm,omega

;  This subroutine computes ephemeris parameters used by other Mission
;  Operations routines:  the solar mean longitude and mean anomaly, and
;  the lunar mean longitude and mean ascending node.  It uses the model
;  referenced in The Astronomical Almanac for 1984, Section S 
;  (Supplement) and documented for the SeaWiFS Project in "Constants
;  and Parameters for SeaWiFS Mission Operations", in TBD.  These
;  parameters are used to compute the solar longitude and the nutation
;  in longitude and obliquity.
;
;	implicit real*8 (a-h,o-z)
;	radeg = 180.d0/3.14159265359d0

;  Sun Mean Longitude 		
	xls = 280.46592d0 + 0.9856473516d0*t
	xls = xls mod 360.d0
 
;  Sun Mean Anomaly		
	gs = 357.52772d0 + 0.9856002831d0*t 
	gs = gs mod 360.d0

;  Moon Mean Longitude		
	xlm = 218.31643d0 + 13.17639648d0*t 
	xlm = xlm mod 360.d0

;  Ascending Node of Moon's Mean Orbit 	
	omega = 125.04452d0 - 0.0529537648d0*t 
	omega = omega mod 360.d0
	
	return
	end

	pro nutate,t,xls,gs,xlm,omega,dpsi,eps

;  This subroutine computes the nutation in longitude and the obliquity
;  of the ecliptic corrected for nutation.  It uses the model referenced
;  in The Astronomical Almanac for 1984, Section S (Supplement) and 
;  documented for the SeaWiFS Project in "Constants and Parameters for 
;  SeaWiFS Mission Operations", in TBD.  These parameters are used to 
;  compute the apparent time correction to the Greenwich Hour Angle and 
;  for the calculation of the geocentric Sun vector.  The input 
;  ephemeris parameters are computed using subroutine ephparms.  Terms 
;  are included to 0.1 arcsecond.
;
;	implicit real*8 (a-h,o-z)
;	radeg = 180.d0/3.14159265359d0

	
;  Nutation in Longitude
	dpsi = - 17.1996*sin(omega/!radeg) 	$
     	 	+ 0.2062*sin(2.*omega/!radeg)	$
     	     	- 1.3187*sin(2.*xls/!radeg) 	$
     		+ 0.1426*sin(gs/!radeg) 	$
     		- 0.2274*sin(2.*xlm/!radeg)
	
;  Mean Obliquity of the Ecliptic	
	epsm = 23.439291d0 - 3.560d-7*t 

;  Nutation in Obliquity 
	deps = 9.2025*cos(omega/!radeg) + 0.5736*cos(2.*xls/!radeg)

;  True Obliquity of the Ecliptic 
	eps = epsm + deps/3600.d0

	dpsi = dpsi/3600.d0
        return
	end


;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    LONLAT2SOLZ calculates the solar and sensor view angles from geodetic
;    longitude and latitude and observation time.
;
; -----------------------------------------------------------------
; pro lonlat2solz
;
; Calculates solar and sensor view angles from geodetic longitude 
; and latitude, and observation time.
;
; Inputs:
;     lon	- longitude of pixel in degrees
;     lat	- longitude of pixel in degrees
;     year      - observation year
;     day       - observation day of year
;     msec      - observation millisecs of day
;
; Outputs:
;     solz	- solar zenith angle of pixel in degrees
;     sola	- solar azimuth angle of pixel in degrees
;
; Notes:
;     Inputs can be scalar or vector, but vectors must be of
;     equal length.
;
; Written By:
;     Bryan Franz, SAIC GSC, NASA/SIMBIOS Project, April 1998.
;     (with much help from geonav.pro by Fred Patt)
;
; -----------------------------------------------------------------

pro lonlat2solz,lon,lat,year,day,msec,solz,sola

    Re = 6378.137        ; Earth radius in km
    f  = 1/298.257       ; Earth flattening factor

    l_sun,year,day,msec/1000.D0,usun

    n    = n_elements(lon)
    solz = fltarr(n)
    sola = fltarr(n)
    rmat = fltarr(3,3) ; Rotation matrix

    for i=0L,n-1 do begin

      rlon   = lon(i)*!pi/180.
      rlat   = lat(i)*!pi/180.

      ;
      ; First, we must define the axes (in ECR space) of a
      ; pixel-local coordinate system, with the z-axis along
      ; the geodetic pixel normal, x-axis pointing east, and
      ; y-axis pointing north.
      ;
      up    = [cos(rlat)*cos(rlon),cos(rlat)*sin(rlon),sin(rlat)]
      upxy  = sqrt(up(0)*up(0)+up(1)*up(1))
      ea    = [-up(1)/upxy,up(0)/upxy,0.0]
      no    = crossp(up,ea)

      ;
      ; Compute geocentric pixel location vector.
      ;
      phi   = atan(tan(rlat)*(1-f)^2)                ; geocentric latitude
      R     = Re*(1-f)/sqrt(1-(2-f)*f*(cos(phi)^2))  ; dist to Earth surface
      gvec  = R*[cos(phi)*cos(rlon),cos(phi)*sin(rlon),sin(phi)]

      ;
      ; Now we can transform Sun vectors into the local frame.
      ;
      rmat(0,*) = ea
      rmat(1,*) = no
      rmat(2,*) = up
      sunl      = rmat # usun(*,i)

      ;
      ; Compute the solar zenith and azimuth
      ;
      solz(i) = atan(sqrt(sunl(0)*sunl(0)+sunl(1)*sunl(1)),sunl(2)) * !radeg
      if (solz(i) gt 0.05) then begin
          sola(i) = atan(sunl(0),sunl(1)) * !radeg
      endif else begin
          sola(i) = 0.0
      endelse
      if (sola(i) lt 0.0) then $
          sola(i) = sola(i) + 360.0d0

    endfor

    if (n eq 1) then begin
        solz = solz(0)
        sola = sola(0)
    endif

    return
end


FUNCTION jd_to_yyyymmdd, jd, hhmm=hhmm

;------------------------------------------------------------------------------
; BY:  Elizabeth GOOD
; DATE: 25 May 2011
;
; PURPOSE:
; Function to convert a Julian day, as per the IDL function JULDAY, to a 
; YYYYMMDD string.
;
; Works on arrays or scalars.
;
; INPUTS:
; jd - Julian day (i.e. that corresponds to the IDL JULDAY function.
;
; MODIFICATIONS:
; 24 November, Elizabeth Good - modified to optionally return hhmm
;------------------------------------------------------------------------------

; Get month, day and year
CALDAT, jd, month, day, year, hour, minute

; Convert to strings
yyyy='0000'+STRCOMPRESS(year, /REMOVE_ALL)
mm='00'+STRCOMPRESS(month, /REMOVE_ALL)
dd='00'+STRCOMPRESS(day, /REMOVE_ALL)
hh='00'+STRCOMPRESS(hour, /REMOVE_ALL)
mi='00'+STRCOMPRESS(minute, /REMOVE_ALL)

; Trim excess digits
yyyy=STRMID(yyyy, 3, 4, /REVERSE_OFFSET)
mm=STRMID(mm, 1, 2, /REVERSE_OFFSET)
dd=STRMID(dd, 1, 2, /REVERSE_OFFSET)
hh=STRMID(hh, 1, 2, /REVERSE_OFFSET)
mi=STRMID(mi, 1, 2, /REVERSE_OFFSET)

hhmm=hh+mi

RETURN, yyyy+mm+dd

END 


FUNCTION yyyymmdd_to_jd, yyyymmdd, hhmm=hhmm

;------------------------------------------------------------------------------
; BY:  Elizabeth GOOD
; DATE: 27 July 2011
;
; PURPOSE:
; Function to convert yyyymmdd to Julian day, as per the IDL function JULDAY.
;
; Works on arrays or scalars.
;
; INPUTS:
; yyyymmdd - STRING or NUMERIC YYYYMMDD
;
; KEYWORD:
; hhmm - optional hhmm string for hour and minute.
;
; Modified, 02/02/12, Elizabeth Good - added hhmm keyword.
;------------------------------------------------------------------------------

; Input as a string
yyyymmdd_new=STRCOMPRESS(yyyymmdd)

; Check length is appropriate.
IF MIN(STRLEN(yyyymmdd_new)) NE 8 OR MAX(STRLEN(yyyymmdd_new)) NE 8 THEN $
   MESSAGE, 'Input <yyyymmdd> must be 8 characters'

year=FIX(STRMID(yyyymmdd_new, 0, 4))
month=FIX(STRMID(yyyymmdd_new, 4, 2))
day=FIX(STRMID(yyyymmdd_new, 6, 2))

IF KEYWORD_SET(hhmm) THEN BEGIN
   hh=FIX(STRMID(hhmm, 0, 2))
   mm=FIX(STRMID(hhmm, 2, 2))
   RETURN, JULDAY(month, day, year, hh, mm)
ENDIF ELSE $
   RETURN, JULDAY(month, day, year)

END 


FUNCTION NCDF_QUICK_READ, file,     $
                          name,     $
                          varname,  $
                          NO_STR_CONV=no_str_conv, $
                          READFILL=readFill, $
                          FILLVALUE=fillValue, $
			  ERROR=error

;+
; NAME:
;       NCDF_QUICK_READ
; PURPOSE:
;       Quickly read a variable or attribute from a NetCDF file
;       without needing to worry about any of the mechanics of
;       opening and reading NetCDF.
; CALLING SEQUENCE:
;       Data = NCDF_QUICK_READ( File, Name [, VarName] [, /NO_STR_CONV] )
; METHOD:
;       The routine uses the built in IDL commands to read the
;       NetCDF file. The user only needs to specify the minimum
;       information about what they want to read; everything else
;       is worked out from the NetCDF file or the user inputs.
; CODE OWNER:
;       Simon Good.
; INPUTS:
;       File - Name of file to read. Alternatively, this can be
;              the file ID of an already open NetCDF file. In this
;              case the file is left open after the read.
;       Name - The name of the variable or attribute to read.
;       VarName - If reading an attribute of a variable use this
;                 to specify the name of that variable. Do not
;                 specify VarName if reading a variable or a 
;                 global attribute.
;       readFill - Set to tell the program to attempt to read
;                  a fill value.
;       fillValue - Used to return the fill value.
; KEYWORDS:
;       NO_STR_CONV - If this is set strings are returned as
;                     an array of bytes.
; OUTPUTS:
;       The function returns the variable or attribute.
; COMMON BLOCK:
;       None.
; RETURN:
;       The data that have been read.
; MODIFICATION HISTORY:
;       See FCM
;       V1.4 - modified by Lizzie Good to return -1 if 
;              specified variable name does not exist.  Accompanied by
;              new keyword 'error' return value of 1 (0 if no errors).
;              (3 June 2015)
; EXAMPLES:
;       Read the variable 'temperature' from Example.nc:
;         Var = NCDF_QUICK_READ('Example.nc','temperature')
;       Read the attribute '_FillValue' for variable 'temperature':
;         Att = NCDF_QUICK_READ('Example.nc','_FillValue','temperature')
;       Read the global attribute 'History'. This is a string but to
;       have the string returned as an array of bytes, set no_str_conv.
;         Att = NCDF_QUICK_READ('Example.nc','History',/NO_STR_CONV)

; NOTES:
;       
;-
;      $URL: svn://fcm9/idl_svn/desktop/trunk/ukmo/lib/ukmo_template $
;   $Author: itkt $
;     $Date: 2008-11-24 11:24:47 +0000 (Mon, 24 Nov 2008) $
; $Revision: 2545 $
;       $Id: ukmo_template 2545 2008-11-24 11:24:47Z itkt $

@comm_error.inc

COMPILE_OPT IDL2,HIDDEN

quiet  = !QUIET
!QUIET = 1

vartype = SIZE(file,                $ 
               /TYPE)

IF vartype EQ 7 THEN BEGIN
    ; file is a string.
    fid = NCDF_OPEN(file,           $
                    /NOWRITE)
ENDIF ELSE BEGIN
    fid = file
ENDELSE

; Default setting no errors - overwritten later.
error=0

; If there are three parameters then this must be an attribute read.
IF N_PARAMS() EQ 3 THEN BEGIN

    varid = NCDF_VARID(fid,     $
                       varname)

    IF varid LT 0 THEN BEGIN

        ; Check if it exists.  Return with error code if not.
	VarDetails = NCDF_ATTINQ( fid,   $
                                  varname, $
                                  /GLOBAL)
        IF VarDetails.datatype EQ 'UNKNOWN' THEN BEGIN
	   error=1
           IF vartype EQ 7 THEN $
              NCDF_CLOSE, fid
           RETURN, -1
	ENDIF
    ENDIF

    VarDetails = NCDF_ATTINQ( fid,   $
                              varid, $
                              name)
    
    IF VarDetails.datatype EQ 'UNKNOWN' THEN BEGIN
	error=1
        IF vartype EQ 7 THEN $
           NCDF_CLOSE, fid
        RETURN, -1
    ENDIF

    NCDF_ATTGET, fid,           $
                 varid,         $
                 name,          $
                 value

ENDIF ELSE BEGIN

    ; Need to decide if this is a var or a global attribute.
    varid = NCDF_VARID(fid,     $
                       name)

    IF varid LT 0 THEN BEGIN

        ; Check if it exists.  Return with error code if not.
	VarDetails = NCDF_ATTINQ( fid,   $
                                  name,  $
                                  /GLOBAL )
        IF VarDetails.datatype EQ 'UNKNOWN' THEN BEGIN
	   error=1
           IF vartype EQ 7 THEN $
              NCDF_CLOSE, fid
           RETURN, -1
	ENDIF
	
        ; It is a global attribute.
        NCDF_ATTGET, fid,       $
                     name,      $
                     value,     $
                     /GLOBAL


    ENDIF ELSE BEGIN

        NCDF_VARGET, fid,       $
                     varid,     $
                     value

        VarDetails = NCDF_VARINQ( fid, $
                                  varid )


    ENDELSE

    if keyword_set(readFill) then $
        fillValue = ncdf_quick_read(fid, '_FillValue', name)

ENDELSE

; Convert into a string if required.
IF NOT KEYWORD_SET(NO_STR_CONV) THEN BEGIN                           
  IF STRUPCASE(VarDetails.DataType) EQ 'CHAR' THEN BEGIN
    value = STRING( value )
  ENDIF
ENDIF

IF vartype EQ 7 THEN $
  NCDF_CLOSE, fid

!QUIET = quiet

RETURN,value
END
                          

PRO generate_lat_long_grid, lat1, lat2, lon1, lon2, latres, lonres, lat_grid, lon_grid, $
                            grid_points=grid_points

;----------------------------------------------------------------------------------------
; AUTHOR:    Elizabeth Noyes
; DATE:      Jan 06
;
; ABOUT THE PRORGRAM:
;
; Procedure to generate a lat and long grid.  Grid cells have fixed resolution, and grid
; lat1 and long1 must be specified.
;
; lat1 - starting latitude of grid (-90 to 90)
; lon1 - starting longitude of grid (0 to 360)
; latres - resolution of grid in latitude (Degrees)
; lonres - resolution of grid in longitude (Degrees)
; lat_grid - returned latitude grid
; lon_grid - returned longitude grid.
;
; KEYWORDS:
; grid_point - allows user to force first and last grid points to be as defined by lat1, lat2, lon1
;              and lon2.
;
; LIMITATIONS:
;
; Assumes grid cells are same dimensions in terms of latitude and longitude degrees.
;
; VERSION 2: Elizabeth Noyes, July 2007
; changed lat_strip=lat1-(lat_array*res) to lat_strip=lat1+(lat_array*res).  Will also now work
; for any grid range (i.e. no whole-earth).

; VERSION 3: Elizabeth Noyes, October 2007
; Addition of keyword: grid_point, which allows the user to force the program to use the lat1
; and lon1 as the first grid point.  Also bug fix to deal with lat and long grids that go from max
; to min in place of min to max (e.g. like ECMWF latitude grid).

; VERSION 4: Elizabeth Good (formerly Noyes), April 2008
; Modified to permit grids with differing lat/long resolutions.
;
; VERSION 5: Elizabeth Good, July 2011
; Modified to remove loops and use rebin instead.
;----------------------------------------------------------------------------------------

; Check whether grid_points keyword is set: if so, we want to start the grid at lat1 and lon1 and
; end at lat2 and lon2.  Define start offset and an add-on factor accordingly.
IF KEYWORD_SET(grid_points) THEN BEGIN
   lat_offset=0.
   lon_offset=0.
   add_on=1
ENDIF ELSE BEGIN
   lat_offset=latres/2.
   lon_offset=lonres/2.
   add_on=0
ENDELSE


; First generate 'strips' for lat and long values as along each row of pixels
; the latitude will be the same, and for each column of data, the longitude will
; be the same.
nlat=ABS((lat2-lat1)/FLOAT(latres))+add_on
nlon=ABS((lon2-lon1)/FLOAT(lonres))+add_on
IF FIX(nlat) NE nlat OR FIX(nlon) NE nlon THEN $
   MESSAGE, 'grid for given lat/long range not possible for given resolution'

; Make strip arrays for lats and longs.  These are then mosaiced together to form a grid.
lat_array=FINDGEN(nlat)
lon_array=FINDGEN(nlon)
IF lat2 GT lat1 THEN BEGIN
   start_lat=lat1+lat_offset
   lat_strip=start_lat+(lat_array*latres)
ENDIF ELSE BEGIN
   start_lat=lat1-lat_offset
   lat_strip=start_lat-(lat_array*latres)
ENDELSE
IF lon2 GT lon1 THEN BEGIN
   start_lon=lon1+lon_offset
   lon_strip=start_lon+(lon_array*lonres)
ENDIF ELSE BEGIN
   start_lon=lon1-lon_offset
   lon_strip=start_lon-(lon_array*lonres)
ENDELSE


; Now replicate each column of latitudes along the whole grid and each row
; of longitudes along the whole grid.  NOTE: NEED TO REVISE CODE BELOW USING
; REBIN.
lat_grid=REBIN(TRANSPOSE(lat_strip), nlon, nlat)
lon_grid=REBIN(lon_strip, nlon, nlat)

END


FUNCTION num_to_string, number, FORMAT=format, INTEGER=integer

; Quick routine to convert a number to a string with formatting.

IF N_ELEMENTS(format) EQ 0 THEN format='(F6.2)'
IF N_ELEMENTS(integer) EQ 1 THEN format='(I50)'

RETURN, STRCOMPRESS(STRING(number, format=format), /REMOVE_ALL)

END
