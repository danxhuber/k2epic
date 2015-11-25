;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Wrapper to make a K2 Ecliptic Plane Input Catalog (EPIC)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; dependencies
@makecat.pro
@makeASCII.pro
@loadcolors.pro

pro epic,ra=ra,dec=dec,inrad=inrad,step=instep,epicid=epicid

PREF_SET, 'IDL_PATH', '+./lib/:<IDL_DEFAULT>', /COMMIT

loadcolors

; campaign 0
;ra = 98.296250
;dec = 21.587778
;inrad = 9.0
;step = 1.

; campaign 1
;ra = 173.93958
;dec = 1.4172222
;inrad = 9.
;step = 1.

;campaign 2
;ra = +246.1264
;dec = -22.4473 
;inrad = 9.0
;step = 1.0

; campaign 3
;ra = 336.67
;dec = -11.1
;inrad = 2.0
;step = 1.0

; campaign 4
;ra = 59.0758
;dec = 18.6606
;inrad = 1.0
;step = 1.0

; campaign 5
;ra = 130.158
;dec = 16.8297
;inrad = 2.0
;step = 1.0

if not keyword_set(epicid) then epicid = long(0) else epicid = long(epicid)
if not keyword_set(step) then step=1
if not keyword_set(ra) or not keyword_set(dec) or not keyword_set(inrad) then begin
    print,'error: must specify ra, dec and catalog radius'
    stop
endif

; name for output files
filename = 'ra'+strtrim(string(ra,format='(f5.1)'),2)+'_de'+$
strtrim(string(dec,format='(f5.1)'),2)+'_r'+strtrim(string(inrad,format='(f5.1)'),2)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; call Python code to download input catalogs from Vizier
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print,'------------------------------------------------------------'
print,'downloading catalogs from Vizier ...'
print,'------------------------------------------------------------'
;t0 = systime(1)
;spawn,'python2.7 getcatalogs.py '+strtrim(string(ra),2)+' '+strtrim(string(dec),2)+' '+strtrim(string(inrad),2)+' '+strtrim(string(step),2)
spawn,'python getcatalogs.py '+strtrim(string(ra),2)+' '+strtrim(string(dec),2)+' '+strtrim(string(inrad),2)+' '+strtrim(string(step),2)
;t1 = systime(1)
;print,'done, took',t1-t0,'seconds'


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; perform cross-matching of input catalogs (and optional stellar classification)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print,'------------------------------------------------------------'
print,'performing Xmatching and classification ...'
print,'------------------------------------------------------------'
makecat,k2cat,epicid=epicid,meanra=ra,meande=dec,maxrad=inrad


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; convert output to MAST ascii format
print,'------------------------------------------------------------'
print,'writing SOC and DMC output ...'
print,'------------------------------------------------------------'
makeascii,k2cat,filename=filename,meanra=ra,meande=dec,maxrad=inrad

stop

end
