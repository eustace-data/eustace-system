; Copied from research/metoffice-sst_to_mat/setup_satellte.pro

function normalise,inarr
  return, inarr/total(inarr*inarr)
end

function distance, x1,y1,x2,y2
  xdiff = x1-x2
  index = where(abs(xdiff) gt 180,ct)
  if ct gt 0 then xdiff[index] = 360-abs(xdiff[index])
  return, sqrt((xdiff)^2 + (y1-y2)^2)
end

pro get_fourier_pentads, neofs, C_ptd, C_daily

if not keyword_set(neofs) then neofs = 4

npts = 365

C_daily = fltarr(npts,neofs+1)
C_daily[*,0] = 1.0
for i =1,neofs,2 do begin
C_daily[*,i] = sin((i+1)*!pi*findgen(npts)/npts)
C_daily[*,i+1] = cos((i+1)*!pi*findgen(npts)/npts)
endfor
for i = 0,neofs do C_daily[*,i] = normalise(C_daily[*,i])

C_ptd = fltarr(73,neofs+1)
for i = 0, neofs do begin
for p = 0,72 do begin
   C_ptd[p,i] = avg(C_daily[p*5:p*5+4,i])
endfor
endfor

return
end

pro get_fourier_dailies, neofs, C_daily1, C_daily2

if not keyword_set(neofs) then neofs = 4

npts = 365

C_daily = fltarr(npts,neofs+1)
C_daily[*,0] = 1.0
for i =1,neofs,2 do begin
C_daily[*,i] = sin((i+1)*!pi*findgen(npts)/npts)
C_daily[*,i+1] = cos((i+1)*!pi*findgen(npts)/npts)
endfor
for i = 0,neofs do C_daily[*,i] = normalise(C_daily[*,i])

C_daily1 = C_daily
C_daily2 = C_daily

return
end

function ymd2dn,yr,m,d, help=hlp
  if (n_params(0) lt 3) or keyword_set(hlp) then begin
     print,' Convert from year, month, day to day number of year.'
     print,' dy = ymd2dn(yr,m,d)'
     print,'   yr = year (like 1988).               in'
     print,'   m = month number (like 11 = Nov).    in'
     print,'   d = day of month (like 5).           in'
     print,'   dy = day number in year (like 310).  out'
     return, -1
  endif
;----  Days before start of each month (non-leap year)  -----
  idays = [0,31,59,90,120,151,181,212,243,273,304,334,366]
;----  Correct for leap year if month ge 3  -------------
  lpyr = (((yr mod 4) eq 0) and ((yr mod 100) ne 0)) $
         or ((yr mod 400) eq 0) and (m ge 3)
  dy = d + idays[m-1] + lpyr
  return, dy
end

FUNCTION leap, year, month = month, days = days
; Function to test if a year is a leap year or not
;returns 1 if it is a leap year and 0 otherwise.
;will also return the number of days in each month as an array (days)
;if
;month keyword is non-zero
  IF (year mod 4 EQ 0) AND (year mod 100 NE 0 OR year mod 400 EQ 0) THEN BEGIN
     result = 1 
  ENDIF ELSE BEGIN
     result = 0
  ENDELSE
  IF KEYWORD_SET(month) THEN BEGIN
     IF result EQ 0 THEN BEGIN
        days = [31,28,31,30,31,30,31,31,30,31,30,31]
     ENDIF ELSE BEGIN
        days = [31,29,31,30,31,30,31,31,30,31,30,31]
     ENDELSE
  ENDIF
  RETURN,result
END


function which_lon_box,lon,res

  if lon gt 180.0 or lon lt -180.0 then stop,'Bad lon',lon

  ;observations on grid boundaries are pushed east or north
  shift_lon = lon + 0.001

  ;work out number of lon and lat boxes
  max_lon_box = 360/res - 1

  lon_box = FLOOR((shift_lon - (-180.0)) / res)

  index = where(lon_box lt 0,ct)
  if ct gt 0 then lon_box(index) = max_lon_box
  index = where(lon_box gt max_lon_box,ct)
  if ct gt 0 then lon_box(index) = 0

  return,lon_box
end


function which_lat_box,lat,res
  
  if lat gt 90.0 or lat lt -90 then stop,'Bad lat',lat

  ;observations on grid boundaries are pushed east or north
  shift_lat = lat + 0.001

  ;work out number of lon and lat boxes
  max_lat_box = 180/res - 1

  lat_box = FLOOR((shift_lat - (-90.0)) / res) 

  index = where(lat_box lt 0,ct)
  if ct gt 0 then lat_box(index) = 0
  index = where(lat_box gt max_lat_box,ct)
  if ct gt 0 then lat_box(index) = max_lat_box

  return,lat_box
end


function which_pentad,year,month,day 

  isleap = leap(year,/month,days=ndays_in_month) ; Get number of days in month and check if year is a leap year

  if day gt ndays_in_month(month-1) then stop,'bad day'

  yd = ymd2dn(year,month,day)

; Shift data from Feb 29th onward back by 1 day if we are in a leap
; year so we compute anomalies relative to 28th Feb 
  IF isleap EQ 1 THEN BEGIN
     if yd gt 59 THEN yd = yd - 1
  ENDIF  
  
  pentad = CEIL(yd / 5.0) - 1   ; Calculate pentad and set to be zero-subscripted

  return,pentad
end

function get_offset,C,mean_params,neofs,lat,lon,day_in_year
  offset = 0.0
  xx = which_lon_box(lon,1)
  yy = which_lat_box(lat,1)
  for i = 0, neofs do offset = offset + C[day_in_year,i] * mean_params[i].data[xx,yy]
  return, offset
end


function get_hires_offset,C,mean_params,neofs,lat,lon,day_in_year
  offset = 0.0
  xx = which_lon_box(lon,0.25)
  yy = which_lat_box(lat,0.25)
  for i = 0, neofs do offset = offset + C[day_in_year,i] * mean_params[i].data[xx,yy]
  return, offset
end

function get_single_hires_offset,C,mean_params,target,neofs,lat,lon,day_in_year
  offset = 0.0
  xx = which_lon_box(lon,0.25)
  yy = which_lat_box(lat,0.25)
  offset = C[day_in_year,target] * mean_params[target].data[xx,yy]
  return, offset
end


