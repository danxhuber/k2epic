;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; generates K2 input catalogs from Hipparcos, Tycho-2, UCAC-4, 2MASS and SDSS DR9
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; subroutines and functions
@progress.pro
@transformations.pro
@torad.pro
@ucacerrs.pro

pro makecat,k2cat,epicid=epicid,$
meanra=meanra,meande=meande,maxrad=maxrad

; load input catalogs
print,'--------------------------------------'
print,'loading input catalogs'

;window,xsize=900,ysize=1000
!p.charsize=2.5
!p.multi=[0,2,3]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Hipparcos 2007
hip = mrdfits('catalogs/hip.fits',1,head)
nhip = n_elements(hip)
plot,hip._raj2000,hip._dej2000,psym=3,$
xrange=[meanra-maxrad,meanra+maxrad],yrange=[meande-maxrad,meande+maxrad],/xs,/ys,$
title='Hipparcos',xtitle='RA (deg)',ytitle='DEC (deg)'

u = where(hip.e_plx ne hip.e_plx or hip.e_plx lt 0.)
if (u[0] ne -1) then begin
	hip[u].plx = 0.
	hip[u].e_plx = 0.
endif

u = where(hip.e_pmra lt 0. or hip.e_pmra ne hip.e_pmra)
if (u[0] ne -1) then begin
	hip[u].pmra = 0.
	hip[u].e_pmra = 0.
endif

u = where(hip.e_pmde lt 0. or hip.e_pmde ne hip.e_pmde)
if (u[0] ne -1) then begin
	hip[u].pmde = 0.
	hip[u].e_pmde = 0.
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Tycho-2
tycho = mrdfits('catalogs/tycho.fits',1,head)
td = intarr(n_elements(tycho))
ntycho = n_elements(tycho)
plot,tycho._raj2000,tycho._dej2000,psym=3,$
xrange=[meanra-maxrad,meanra+maxrad],yrange=[meande-maxrad,meande+maxrad],/xs,/ys,$
title='Tycho-2',xtitle='RA (deg)',ytitle='DEC (deg)'

; replace NaN's
u = where(tycho.e_btmag lt 0. or tycho.e_btmag ne tycho.e_btmag)
if (u[0] ne -1) then begin
	tycho[u].btmag = 0.
	tycho[u].e_btmag = 0.
endif

u = where(tycho.e_vtmag lt 0. or tycho.e_vtmag ne tycho.e_vtmag)
if (u[0] ne -1) then begin
	tycho[u].vtmag = 0.
	tycho[u].e_vtmag = 0.
endif

u = where(tycho.e_pmra lt 0. or tycho.e_pmra ne tycho.e_pmra)
if (u[0] ne -1) then begin
	tycho[u].pmra = 0.
	tycho[u].e_pmra = 0.
endif

u = where(tycho.e_pmde lt 0. or tycho.e_pmde ne tycho.e_pmde)
if (u[0] ne -1) then begin
	tycho[u].pmde = 0.
	tycho[u].e_pmde = 0.
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; UCAC-4
ucac = mrdfits('catalogs/ucac.fits',1,head)
plot,ucac._raj2000,ucac._dej2000,psym=3,$
xrange=[meanra-maxrad,meanra+maxrad],yrange=[meande-maxrad,meande+maxrad],/xs,/ys,$
title='UCAC-4',xtitle='RA (deg)',ytitle='DEC (deg)'

; replace NaN's
u = where(ucac.e_bmag lt 0.)
if (u[0] ne -1) then begin
	ucac[u].bmag = 0.
	ucac[u].e_bmag = 0.
endif

u = where(ucac.e_vmag lt 0.)
if (u[0] ne -1) then begin
	ucac[u].vmag = 0.
	ucac[u].e_vmag = 0.
endif

u = where(ucac.e_gmag lt 0.)
if (u[0] ne -1) then begin
	ucac[u].gmag = 0.
	ucac[u].e_gmag = 0.
endif

u = where(ucac.e_rmag lt 0.)
if (u[0] ne -1) then begin
	ucac[u].rmag = 0.
	ucac[u].e_rmag = 0.
endif

u = where(ucac.e_imag lt 0.)
if (u[0] ne -1) then begin
	ucac[u].imag = 0.
	ucac[u].e_imag = 0.
endif

u = where(ucac.e_bmag lt 0.)
if (u[0] ne -1) then begin
	ucac[u].bmag = 0.
	ucac[u].e_bmag = 0.
endif

u = where(ucac.pmra ne ucac.pmra)
if (u[0] ne -1) then begin
	ucac[u].pmra = 0.
	ucac[u].e_pmra = 0.
endif

u = where(ucac.pmde ne ucac.pmde)
if (u[0] ne -1) then begin
	ucac[u].pmde = 0.
	ucac[u].e_pmde = 0.
endif


; some UCAC errors are unrealistically small; replace those that are 'formal' and zero 
; with typical values based on fit of real uncertainties
u = where((strmatch(ucac.f_bmag,'-') and ucac.e_bmag eq 1.) or ucac.e_bmag eq 0.)
ucac[u].e_bmag = ucacerrs(ucac[u].bmag,/bmag) 

u = where((strmatch(ucac.f_vmag,'-') and ucac.e_vmag eq 1.) or ucac.e_vmag eq 0.)
ucac[u].e_vmag = ucacerrs(ucac[u].vmag,/vmag) 

u = where((strmatch(ucac.f_gmag,'-') and ucac.e_gmag eq 1.) or ucac.e_gmag eq 0.)
ucac[u].e_gmag = ucacerrs(ucac[u].gmag,/gmag) 

u = where((strmatch(ucac.f_rmag,'-') and ucac.e_rmag eq 1.) or ucac.e_rmag eq 0.)
ucac[u].e_rmag = ucacerrs(ucac[u].rmag,/rmag) 

u = where((strmatch(ucac.f_imag,'-') and ucac.e_imag eq 1.) or ucac.e_imag eq 0.)
ucac[u].e_imag = ucacerrs(ucac[u].imag,/imag) 

;orig = ucac
; convert APASS g'r'i' to Sloan gri
u = where(ucac.imag ne 0. or ucac.gmag ne 0. or ucac.rmag ne 0.)
for q=0.,n_elements(u)-1 do begin
    res = tosloan(ucac[u[q]].gmag,ucac[u[q]].rmag,ucac[u[q]].imag)
    ucac[u[q]].gmag = res[0]
    ucac[u[q]].rmag = res[1]
    ucac[u[q]].imag = res[2]
end

ud = intarr(n_elements(ucac))
nucac = n_elements(ucac)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 2MASS
mass = mrdfits('catalogs/2mass.fits',1,head)

dis = sqrt( (mass._raj2000-meanra)^2D + (mass._dej2000-meande)^2D )
u = where(dis le 9.)
mass = mass[u]

plot,mass._raj2000,mass._dej2000,psym=3,$
xrange=[meanra-maxrad,meanra+maxrad],yrange=[meande-maxrad,meande+maxrad],/xs,/ys,$
title='2MASS',xtitle='RA (deg)',ytitle='DEC (deg)'

jqflg = strmid(mass.qflg,0,1)
; only keep sources that are detected in J-band
u = where(strmatch(jqflg,'A') or strmatch(jqflg,'B') or strmatch(jqflg,'C'))
mass = mass[u]
jqflg = strmid(mass.qflg,0,1)
hqflg = strmid(mass.qflg,1,1)
kqflg = strmid(mass.qflg,2,1)

hgood = replicate(0.,n_elements(hqflg))
u = where(strmatch(hqflg,'A') or strmatch(hqflg,'B') or strmatch(hqflg,'C'))
hgood[u] = 1.
kgood = replicate(0.,n_elements(kqflg))
u = where(strmatch(kqflg,'A') or strmatch(kqflg,'B') or strmatch(kqflg,'C'))
kgood[u] = 1.

md = intarr(n_elements(mass))
nmass = n_elements(mass)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SDSS
nosdss = 0
if file_test('catalogs/sdss.fits') then begin
sdss = mrdfits('catalogs/sdss.fits',1,head)
; eliminate bad SDSS photometry with large errors or bad quality flags
u = where(sdss.e_gmag lt 0.5 and sdss.e_rmag lt 0.5 and sdss.e_imag lt 0.5 and strmatch(sdss.q_mode,'*+*'))
sdss=sdss[u]

plot,sdss._raj2000,sdss._dej2000,psym=3,$
xrange=[meanra-maxrad,meanra+maxrad],yrange=[meande-maxrad,meande+maxrad],/xs,/ys,$
title='SDSS DR9',xtitle='RA (deg)',ytitle='DEC (deg)'

; replace NaNs
u = where(sdss.pmra ne sdss.pmra or sdss.e_pmra lt 0.)
if (u[0] ne -1) then begin
sdss[u].pmra = 0.
sdss[u].e_pmra = 0.
endif
u = where(sdss.pmde ne sdss.pmde or sdss.e_pmde lt 0.)
if (u[0] ne -1) then begin
sdss[u].pmde = 0.
sdss[u].e_pmde = 0.
endif
u = where(sdss.umag ne sdss.umag or sdss.e_umag lt 0.)
if (u[0] ne -1) then begin
sdss[u].umag = 0.
sdss[u].e_umag = 0.
endif
u = where(sdss.gmag ne sdss.gmag or sdss.e_gmag lt 0.)
if (u[0] ne -1) then begin
sdss[u].gmag = 0.
sdss[u].e_gmag = 0.
endif
u = where(sdss.rmag ne sdss.rmag or sdss.e_rmag lt 0.)
if (u[0] ne -1) then begin
sdss[u].rmag = 0.
sdss[u].e_rmag = 0.
endif
u = where(sdss.imag ne sdss.imag or sdss.e_imag lt 0.)
if (u[0] ne -1) then begin
sdss[u].imag = 0.
sdss[u].e_imag = 0.
endif
u = where(sdss.zmag ne sdss.zmag or sdss.e_zmag lt 0.)
if (u[0] ne -1) then begin
sdss[u].zmag = 0.
sdss[u].e_zmag = 0.
endif

sd = intarr(n_elements(sdss))
nsdss = n_elements(sdss)
endif else begin
	nosdss = 1
	nsdss = 0
endelse


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NOMAD
nomad = mrdfits('catalogs/nomad.fits',1,head)
; only keep 2MASS cross-matched NOMAD sources, since USNO-B sources aren't really reliable
u = where(strmatch(nomad.ym,'*M*'))
nomad = nomad[u]
plot,nomad._raj2000,nomad._dej2000,psym=3,$
xrange=[meanra-maxrad,meanra+maxrad],yrange=[meande-maxrad,meande+maxrad],/xs,/ys,$
title='NOMAD proper motions',xtitle='RA (deg)',ytitle='DEC (deg)'

;saveimage,'gifs/'+filename+'.input.gif'

print,''


!p.multi=0

; convert Tycho BtVt into Johnson BV (only if both B and V are valid)
u = where(tycho.btmag ne 0. and tycho.vtmag ne 0.)

readcol,'Tycho_BV_Bessel2000.dat',bvt,vvt,dbv,/silent,skipline=2
bv = interpol(dbv,bvt,tycho[u].btmag-tycho[u].vtmag)+tycho[u].btmag-tycho[u].vtmag
v = interpol(vvt,bvt,tycho[u].btmag-tycho[u].vtmag)+tycho[u].vtmag
b = bv + v
tycho[u].vtmag = v
tycho[u].btmag = b

; max number of entries
nmax = nhip+ntycho+nucac+nmass+nsdss

; master struct
k2cat = replicate($
{obs: $
	{ra:0D, de:0D, bmag:0D, vmag:0D, umag:0D, gmag:0D, rmag:0D, $
	imag:0D, zmag:0D, jmag:0D, hmag:0D, kmag:0D, pmra:0D, pmde:0D, plx:0D, $
	e_bmag:0D, e_vmag:0D, e_umag:0D, e_gmag:0D, e_rmag:0D, $
	e_imag:0D, e_zmag:0D, e_jmag:0D, e_hmag:0D, e_kmag:0D, e_pmra:0D, e_pmde:0D, e_plx:0D, $
	mflg:'', prox:0D}, $
id: $
	{epic:'', hip:'', tycho:'', ucac:'', mass:'', sdss:'', nomad:'', kepflag:'', stparflag:'', $
	class:'' }, $
deriv: $
	{teff:dblarr(4), logg:dblarr(4), feh:dblarr(4), rad:dblarr(4), mass:dblarr(4), $
	rho:dblarr(4), lum:dblarr(4), age:dblarr(4), dis:dblarr(4), ebv:0D, kp:0D }},nmax)

ix=0.
count=0.

k2cat.id.class = replicate('STAR',nmax)

print,'--------------------------------------'
print,'cross-matching input catalogs'

; don't trust 2MASS photometry below this J magnitude
masslim = 5.

; magnitude matching criteria; somewhat arbitrary
maglim = 1.5
maglim2 = 4.

; match radius
range = 3./60./60.;0.001

; loop over Hipparcos stars
;goto,skiphip

; this calcuates the closest Eucledian distance match between catalogs
match_mass=MATCH_2D(hip._raj2000,hip._dej2000,mass._raj2000,mass._dej2000,range)
match_ucac=MATCH_2D(hip._raj2000,hip._dej2000,ucac._raj2000,ucac._dej2000,range)
if not (nosdss) then match_sdss=MATCH_2D(hip._raj2000,hip._dej2000,sdss._raj2000,sdss._dej2000,range)

for i=0.,n_elements(hip)-1 do begin

	t1 = systime(1)

	k2cat[ix].obs.ra = hip[i]._raj2000
	k2cat[ix].obs.de = hip[i]._dej2000

	k2cat[ix].obs.plx = hip[i].plx
	k2cat[ix].obs.e_plx = hip[i].e_plx
	k2cat[ix].obs.pmra = hip[i].pmra
	k2cat[ix].obs.e_pmra = hip[i].e_pmra
	k2cat[ix].obs.pmde = hip[i].pmde
	k2cat[ix].obs.e_pmde = hip[i].e_pmde
	k2cat[ix].id.hip = 'HIP:'+strtrim(string(hip[i].hip),2)
	
	; cross-match Tycho
	u = where(hip[i].hip eq tycho.hip)
	if (u[0] ne -1) then begin
		k2cat[ix].obs.bmag = tycho[u[0]].btmag
		k2cat[ix].obs.e_bmag = tycho[u[0]].e_btmag
		k2cat[ix].obs.vmag = tycho[u[0]].vtmag
		k2cat[ix].obs.e_vmag = tycho[u[0]].e_vtmag
		k2cat[ix].id.tycho = 'TYC:'+strtrim(tycho[u[0]].tyc1,2)+'-'+strtrim(tycho[u[0]].tyc2,2)$
		+'-'+strtrim(tycho[u[0]].tyc3,2)
		tyc = u[[0]]
		td[u] = 1
	endif else tyc = -1
	
	; cross-match 2MASS
	if (match_mass[i] ne -1) then begin
	    pos = match_mass[i]

    	    ; estimate V mag from JHK of the matching object; if there is no V mag to compare with, assume sources match
	    if (k2cat[ix].obs.vmag ne 0. and kgood[pos] and hgood[pos]) then vmatch = jhktov(mass[pos].jmag,mass[pos].hmag,mass[pos].kmag) else $
	    vmatch = k2cat[ix].obs.vmag

	    if (mass[pos].jmag gt masslim and abs(k2cat[ix].obs.vmag-vmatch) lt maglim) then begin
		k2cat[ix].obs.jmag = mass[pos].jmag
		k2cat[ix].obs.e_jmag = mass[pos].e_jmag	
    	    	if (hgood[pos]) then begin
		    k2cat[ix].obs.hmag = mass[pos].hmag
		    k2cat[ix].obs.e_hmag = mass[pos].e_hmag
		endif
		if (kgood[pos]) then begin
		    k2cat[ix].obs.kmag = mass[pos].kmag
		    k2cat[ix].obs.e_kmag = mass[pos].e_kmag	
		endif
		k2cat[ix].id.mass = '2MASS:'+strtrim(string(mass[pos]._2mass),2)
		k2cat[ix].obs.mflg = strtrim(string(mass[pos].qflg),2)+'-'+$
		strtrim(string(mass[pos].rflg),2)+'-'+$
		strtrim(string(mass[pos].bflg),2)+'-'+$
		strtrim(string(mass[pos].cflg),2)+'-'+$
		strtrim(string(fix(mass[pos].xflg)),2)+'-'+$
		strtrim(string(fix(mass[pos].aflg)),2)	
		k2cat[ix].obs.prox = strtrim(string(mass[pos].prox),2)	
		
		md[pos] = 1
		
	    endif
	endif
		
	; cross-match UCAC
	if (match_ucac[i] ne -1) then begin
	    pos = match_ucac[i]
	    
	    ; first check: compare if there is overlapping V-band; if not, assume sources match
	    if (k2cat[ix].obs.vmag ne 0. and ucac[pos].vmag ne 0.) then diff = abs(k2cat[ix].obs.vmag-ucac[pos].vmag) else diff = 0.
	    
	    ; second check: see whether UCAC-Tycho cross-match exists
	    tycmatch = 0
	    if (tyc[0] ne -1) then $
	    if (strmatch(ucac[pos].tycho_2,'*'+strtrim(tycho[tyc].tyc1,2)+'*') $
	    and strmatch(ucac[pos].tycho_2,'*'+strtrim(tycho[tyc].tyc2,2)+'*') and $
	    strmatch(ucac[pos].tycho_2,'*'+strtrim(tycho[tyc].tyc3,2)+'*')) then tycmatch = 1 else tycmatch = 0
	    	    		    
	    ; accept match only if the HIP-TYC TYC-UCAC link matches, or the magnitudes agree
	    if (tycmatch or (diff lt maglim)) then begin
	   ; if (ucac[pos].e_gmag lt 0.1 and ucac[pos].e_rmag lt 0.1 and ucac[pos].e_imag lt 0.1) then begin
		k2cat[ix].obs.gmag = ucac[pos].gmag
		k2cat[ix].obs.e_gmag = ucac[pos].e_gmag/100D
		k2cat[ix].obs.rmag = ucac[pos].rmag
		k2cat[ix].obs.e_rmag = ucac[pos].e_rmag/100D
		k2cat[ix].obs.imag = ucac[pos].imag
		k2cat[ix].obs.e_imag = ucac[pos].e_imag/100D
		
		if (k2cat[ix].obs.vmag eq 0.) then begin
		    k2cat[ix].obs.vmag = ucac[pos].vmag
		    k2cat[ix].obs.e_vmag = ucac[pos].e_vmag/100D
		endif
		
		if (k2cat[ix].obs.bmag eq 0.) then begin
		    k2cat[ix].obs.bmag = ucac[pos].bmag
		    k2cat[ix].obs.e_bmag = ucac[pos].e_bmag/100D
		endif
		
		k2cat[ix].id.ucac = 'UCAC:'+strtrim(string(ucac[pos].ucac4),2)
		if (ucac[pos].leda gt 0.) then k2cat[ix].id.class = 'EXTENDED'
		ud[pos] = 1		
	    endif
	    
	endif
    	
	; match SDSS, but don't record photometry since these will all be saturated
	if not (nosdss) then if (match_sdss[i] ne -1) then begin 
    	    pos = match_sdss[i]
	    
	    ; estimate gmag from BV; if BV doesn;t exist, assume that they match
	    if (k2cat[ix].obs.vmag ne 0. and k2cat[ix].obs.bmag ne 0.) then gmatch = bvtog(k2cat[ix].obs.bmag,k2cat[ix].obs.vmag) else gmatch = sdss[pos].gmag
	
	    if (abs(gmatch-sdss[pos].gmag) lt maglim) then begin
	    	k2cat[ix].id.sdss = 'SDSS:'+strtrim(string(sdss[pos].sdss_id),2)	    
	    	;if (sdss[pos].cl eq 3) then k2cat[ix].id.class = 'GALAXY'
	    	sd[pos] = 1	
	    endif
	endif	
	
	count++
	ix++
		
	t2 = systime(1)
	;print,t2-t1
	progress,'cross-matching Hipparcos',i,nhip
	
endfor
progress,'cross-matching Hipparcos',i,nhip
print,'used ',count,' of ',n_elements(hip),' Hipparcos stars'
print,''
skiphip:

count=0.
	
;goto,skiptyc
; loop over Tycho stars

match_mass=MATCH_2D(tycho._raj2000,tycho._dej2000,mass._raj2000,mass._dej2000,range)
match_ucac=MATCH_2D(tycho._raj2000,tycho._dej2000,ucac._raj2000,ucac._dej2000,range)
if not (nosdss) then  match_sdss=MATCH_2D(tycho._raj2000,tycho._dej2000,sdss._raj2000,sdss._dej2000,range)

;dumarr1 = dblarr(n_elements(tycho))
;dumarr2 = dblarr(n_elements(tycho))

for i=0.,n_elements(tycho)-1 do begin
		
	if (td[i] eq 1) then continue
	;print,i
	;t1 = systime(1)	
	
	k2cat[ix].obs.ra = tycho[i]._raj2000
	k2cat[ix].obs.de = tycho[i]._dej2000

	k2cat[ix].obs.bmag = tycho[i].btmag
	k2cat[ix].obs.e_bmag = tycho[i].e_btmag
	k2cat[ix].obs.vmag = tycho[i].vtmag
	k2cat[ix].obs.e_vmag = tycho[i].e_vtmag

	k2cat[ix].obs.pmra = tycho[i].pmra
	k2cat[ix].obs.e_pmra = tycho[i].e_pmra
	k2cat[ix].obs.pmde = tycho[i].pmde
	k2cat[ix].obs.e_pmde = tycho[i].e_pmde
	k2cat[ix].id.tycho = 'TYC:'+strtrim(tycho[i].tyc1,2)+'-'+strtrim(tycho[i].tyc2,2)+'-'+strtrim(tycho[i].tyc3,2)

	; cross-match 2MASS
	if (match_mass[i] ne -1) then begin
	    pos = match_mass[i]
	    if (k2cat[ix].obs.vmag ne 0. and hgood[pos] and kgood[pos]) then vmatch = jhktov(mass[pos].jmag,mass[pos].hmag,mass[pos].kmag) else $
	    vmatch = k2cat[ix].obs.vmag
	   
	 ;  dumarr1[i] = k2cat[ix].obs.vmag-vmatch
	   
	    ; Hipparcos stars are bright, so discard those that might have messed up 2MASS photometry
	   if (mass[pos].jmag gt masslim and abs(k2cat[ix].obs.vmag-vmatch) lt maglim) then begin
		k2cat[ix].obs.jmag = mass[pos].jmag
		k2cat[ix].obs.e_jmag = mass[pos].e_jmag	
    	    	if (hgood[pos]) then begin
		    k2cat[ix].obs.hmag = mass[pos].hmag
		    k2cat[ix].obs.e_hmag = mass[pos].e_hmag
		endif
		if (kgood[pos]) then begin
		    k2cat[ix].obs.kmag = mass[pos].kmag
		    k2cat[ix].obs.e_kmag = mass[pos].e_kmag	
		endif
		k2cat[ix].id.mass = '2MASS:'+strtrim(string(mass[pos]._2mass),2)
		k2cat[ix].obs.mflg = strtrim(string(mass[pos].qflg),2)+'-'+$
		strtrim(string(mass[pos].rflg),2)+'-'+$
		strtrim(string(mass[pos].bflg),2)+'-'+$
		strtrim(string(mass[pos].cflg),2)+'-'+$
		strtrim(string(fix(mass[pos].xflg)),2)+'-'+$
		strtrim(string(fix(mass[pos].aflg)),2)	
		k2cat[ix].obs.prox = strtrim(string(mass[pos].prox),2)			
		
		md[pos] = 1
	    endif
	endif
	
	; cross-match UCAC
	if (match_ucac[i] ne -1) then begin
	    pos = match_ucac[i]	    
	    if (strmatch(ucac[pos].tycho_2,'*'+strtrim(tycho[i].tyc1,2)+'*') $
	    and strmatch(ucac[pos].tycho_2,'*'+strtrim(tycho[i].tyc2,2)+'*') and strmatch(ucac[pos].tycho_2,'*'+strtrim(tycho[i].tyc3,2)+'*')) then begin
	   ; if (ucac[pos].e_gmag lt 0.1 and ucac[pos].e_rmag lt 0.1 and ucac[pos].e_imag lt 0.1) then begin
		k2cat[ix].obs.gmag = ucac[pos].gmag
		k2cat[ix].obs.e_gmag = ucac[pos].e_gmag/100D
		k2cat[ix].obs.rmag = ucac[pos].rmag
		k2cat[ix].obs.e_rmag = ucac[pos].e_rmag/100D
		k2cat[ix].obs.imag = ucac[pos].imag
		k2cat[ix].obs.e_imag = ucac[pos].e_imag/100D
		k2cat[ix].id.ucac = 'UCAC:'+strtrim(string(ucac[pos].ucac4),2)
		if (ucac[pos].leda gt 0.) then k2cat[ix].id.class = 'EXTENDED'
		ud[pos] = 1		
	    endif
	endif
    	
	;if strmatch(k2cat[ix].id.tycho,'*272*358*1*') then stop
	
	; match SDSS, but don't record photometry since these will all be saturated
	if not (nosdss) then if (match_sdss[i] ne -1) then begin 
    	    pos = match_sdss[i]
	    if (k2cat[ix].obs.bmag ne 0. and k2cat[ix].obs.vmag ne 0.) then gmatch = bvtog(k2cat[ix].obs.bmag,k2cat[ix].obs.vmag) else gmatch = sdss[pos].gmag
	    if (abs(gmatch-sdss[pos].gmag) lt maglim) then begin
	    	k2cat[ix].id.sdss = 'SDSS:'+strtrim(string(sdss[pos].sdss_id),2)	    
	    	;if (sdss[pos].cl eq 3) then k2cat[ix].id.class = 'GALAXY'
	    	sd[pos] = 1	
	    endif
	endif
	
	count++
	ix++

	;t2 = systime(1)
	;print,t2-t1
	;an=''
	;read,an
	progress,'cross-matching Tycho',i,ntycho
	
endfor
progress,'cross-matching Tycho',i,ntycho
print,'used ',count,' of ',n_elements(tycho),' Tycho stars'
print,''
skiptyc:
count=0.


if not (nosdss) then match_sdss=MATCH_2D(ucac._raj2000,ucac._dej2000,sdss._raj2000,sdss._dej2000,range)
match_mass=MATCH_2D(ucac._raj2000,ucac._dej2000,mass._raj2000,mass._dej2000,range)

;openw,5,'compphot.dat'

dumarr = dblarr(n_elements(ucac))

;goto,skipucac
; loop over UCAC stars
for i=0.,n_elements(ucac)-1 do begin
		
	if (ud[i] eq 1) then continue
	
	t1 = systime(1)	
	
	k2cat[ix].obs.ra = ucac[i]._raj2000
	k2cat[ix].obs.de = ucac[i]._dej2000

    k2cat[ix].obs.bmag = ucac[i].bmag
	k2cat[ix].obs.e_bmag = ucac[i].e_bmag/100D
	k2cat[ix].obs.vmag = ucac[i].vmag
	k2cat[ix].obs.e_vmag = ucac[i].e_vmag/100D

	k2cat[ix].obs.gmag = ucac[i].gmag
	k2cat[ix].obs.e_gmag = ucac[i].e_gmag/100D
	k2cat[ix].obs.rmag = ucac[i].rmag
	k2cat[ix].obs.e_rmag = ucac[i].e_rmag/100D
	k2cat[ix].obs.imag = ucac[i].imag
	k2cat[ix].obs.e_imag = ucac[i].e_imag/100D

	k2cat[ix].obs.pmra = ucac[i].pmra
	k2cat[ix].obs.e_pmra = ucac[i].e_pmra
	k2cat[ix].obs.pmde = ucac[i].pmde
	k2cat[ix].obs.e_pmde = ucac[i].e_pmde
	k2cat[ix].id.ucac = 'UCAC:'+strtrim(string(ucac[i].ucac4),2)
	
	; UCAC galaxy cross-match
	if (ucac[i].leda gt 0.) then k2cat[ix].id.class = 'EXTENDED'
	
	; cross-match 2MASS
	if (match_mass[i] ne -1) then begin
	    pos = match_mass[i]
	    if ((ucac[i]._2MKEY eq mass[pos].CNTR)) then begin
		k2cat[ix].obs.jmag = mass[pos].jmag
		k2cat[ix].obs.e_jmag = mass[pos].e_jmag	
    	if (hgood[pos]) then begin
		    k2cat[ix].obs.hmag = mass[pos].hmag
		    k2cat[ix].obs.e_hmag = mass[pos].e_hmag
		endif
		if (kgood[pos]) then begin
		    k2cat[ix].obs.kmag = mass[pos].kmag
		    k2cat[ix].obs.e_kmag = mass[pos].e_kmag	
		endif
		k2cat[ix].id.mass = '2MASS:'+strtrim(string(mass[pos]._2mass),2)
		k2cat[ix].obs.mflg = strtrim(string(mass[pos].qflg),2)+'-'+$
		strtrim(string(mass[pos].rflg),2)+'-'+$
		strtrim(string(mass[pos].bflg),2)+'-'+$
		strtrim(string(mass[pos].cflg),2)+'-'+$
		strtrim(string(fix(mass[pos].xflg)),2)+'-'+$
		strtrim(string(fix(mass[pos].aflg)),2)	
		k2cat[ix].obs.prox = strtrim(string(mass[pos].prox),2)		
		
		md[pos] = 1
	    endif
	endif
    	
	; keep UCAC gri photometry even if there is a SDSS match, since SDSS might be saturated
	if not (nosdss) then if (match_sdss[i] ne -1) then begin 
    	    pos = match_sdss[i]
	    ; if UCAC gmag is 0, try proxy from BV mag
	    gmatch = ucac[i].gmag
	    if (gmatch eq 0. and ucac[i].vmag ne 0. and ucac[i].bmag ne 0.) then gmatch = bvtog(k2cat[ix].obs.bmag,k2cat[ix].obs.vmag)
	    if (gmatch eq 0.) then gmatch = sdss[pos].gmag
	;    dumarr[i] = gmatch-sdss[pos].gmag
	    if (abs(gmatch-sdss[pos].gmag) lt maglim) then begin
	    	k2cat[ix].id.sdss = 'SDSS:'+strtrim(string(sdss[pos].sdss_id),2)	    
	    	if (sdss[pos].cl eq 3) then k2cat[ix].id.class = 'EXTENDED'
		
		;printf,5,ucac[i].gmag,ucac[i].rmag,ucac[i].imag,sdss[pos].gmag,sdss[pos].rmag,sdss[pos].imag
		;printf,5,ucac[i].gmag,k2cat[ix].obs.de,$
		;ucac[i].pmde,ucac[i].e_pmde,ucac[i].pmra,ucac[i].e_pmra,$
		;sdss[pos].pmde,sdss[pos].e_pmde,sdss[pos].pmra,sdss[pos].e_pmra,$
		;format='(d10.3,d10.6,d10.4,d10.4,d10.4,d10.4,d10.4,d10.4,d10.4,d10.4)'
		
	    	sd[pos] = 1	
	    endif
	endif
	
	count++
	ix++

	;progress,'cross-matching UCAC',i,nucac
	
	;t2 = systime(1)
	;print,t2-t1
	;an=''
	;read,an
endfor

progress,'cross-matching UCAC',i,nucac
print,'used ',count,' of ',n_elements(ucac),' UCAC stars'
print,''

;close,5

skipucac:

count=0.

if not (nosdss) then  match_sdss=MATCH_2D(mass._raj2000,mass._dej2000,sdss._raj2000,sdss._dej2000,range)

;dumc=0.

dumarr = dblarr(n_elements(mass))

;ctemp=0.
;ctemp2=0.
; loop over 2MASS stars
;goto,skipmass
for i=0.,n_elements(mass)-1 do begin	
	if (md[i] eq 1) then continue
	
	; skip if there are bad photometry flags
	if (mass[i].jmag gt masslim) then begin
	
	t1 = systime(1)	
	
	k2cat[ix].obs.ra = mass[i]._raj2000
	k2cat[ix].obs.de = mass[i]._dej2000
	k2cat[ix].id.mass = '2MASS:'+strtrim(string(mass[i]._2mass),2)
	k2cat[ix].obs.mflg = strtrim(string(mass[i].qflg),2)+'-'+$
		strtrim(string(mass[i].rflg),2)+'-'+$
		strtrim(string(mass[i].bflg),2)+'-'+$
		strtrim(string(mass[i].cflg),2)+'-'+$
		strtrim(string(fix(mass[i].xflg)),2)+'-'+$
		strtrim(string(fix(mass[i].aflg)),2)	
	k2cat[ix].obs.prox = strtrim(string(mass[i].prox),2)	
	
    k2cat[ix].obs.jmag = mass[i].jmag
	k2cat[ix].obs.e_jmag = mass[i].e_jmag
	if (hgood[i]) then begin
	    k2cat[ix].obs.hmag = mass[i].hmag
	    k2cat[ix].obs.e_hmag = mass[i].e_hmag
	endif
	if (kgood[i]) then begin
	    k2cat[ix].obs.kmag = mass[i].kmag
	    k2cat[ix].obs.e_kmag = mass[i].e_kmag	
    	endif
	md[i] = 1
	
	if (mass[i].xflg eq 1) then k2cat[ix].id.class = 'EXTENDED'

    	;if strmatch(mass[i]._2mass,'*11193805-0256024*') then stop

    	if not (nosdss) then if (match_sdss[i] ne -1) then begin 
	    pos = match_sdss[i]
;	    ctemp++
	    if (hgood[i] and kgood[i]) then gmatch = jhktog(mass[i].jmag,mass[i].hmag,mass[i].kmag) else begin
		 gmatch = sdss[pos].gmag
;		 ctemp2++
	    endelse
	    
	   ; if (mass[i].xflg ne 1 and sdss[pos].cl ne 3 and abs(gmatch-sdss[pos].gmag) gt maglim) then dumc++
    	 ;   if (mass[i].xflg ne 1 and sdss[pos].cl ne 3) then dumarr[i] = gmatch-sdss[pos].gmag

    	;    if (dumarr[i] lt -5. or dumarr[i] gt 5.) then stop

	    ; if it's a galaxy, the color transformation match isn't valid
	    ;gmatch = sdss[pos].gmag
	    if (abs(gmatch-sdss[pos].gmag) lt maglim2 or mass[i].xflg eq 1 or sdss[pos].cl eq 3) then begin
	
	    k2cat[ix].id.sdss = 'SDSS:'+strtrim(string(sdss[pos].sdss_id),2)
	    
	    k2cat[ix].obs.gmag = sdss[pos].gmag
	    k2cat[ix].obs.e_gmag = sdss[pos].e_gmag
	    k2cat[ix].obs.rmag = sdss[pos].rmag
	    k2cat[ix].obs.e_rmag = sdss[pos].e_rmag
	    k2cat[ix].obs.imag = sdss[pos].imag
	    k2cat[ix].obs.e_imag = sdss[pos].e_imag
    	    
	    k2cat[ix].obs.umag = sdss[pos].umag
	    k2cat[ix].obs.e_umag = sdss[pos].e_umag
	    k2cat[ix].obs.zmag = sdss[pos].zmag
	    k2cat[ix].obs.e_zmag = sdss[pos].e_zmag
	    
	    k2cat[ix].obs.pmra = sdss[pos].pmra
	    k2cat[ix].obs.e_pmra = sdss[pos].e_pmra
	    k2cat[ix].obs.pmde = sdss[pos].pmde
	    k2cat[ix].obs.e_pmde = sdss[pos].e_pmde
	    
	    if (sdss[pos].cl eq 3) then k2cat[ix].id.class = 'EXTENDED'
	        
	    sd[pos] = 1	
	    
	    endif ;else begin
	    ;	print,mass[i]._2mass,mass[i]._raj2000,mass[i]._dej2000
	;	print,sdss[pos].sdss_id,gmatch,sdss[pos].gmag
	;	stop
	  ;  endelse
	    
	endif
	
	progress,'cross-matching 2MASS',i,nmass
	count++
	ix++
	
	endif

endfor
;print,ctemp
;print,ctemp2
;stop
;print,dumc

progress,'cross-matching 2MASS',i,nmass
print,'used ',count,' of ',n_elements(mass),' 2MASS stars'
print,''
skipmass:
count=0.
	
; loop over SDSS stars
if not (nosdss) then begin
for i=0.,n_elements(sdss)-1 do begin
	
	if (sd[i] eq 1) then continue
	
	k2cat[ix].obs.ra = sdss[i]._raj2000
	k2cat[ix].obs.de = sdss[i]._dej2000

	k2cat[ix].id.sdss = 'SDSS:'+strtrim(string(sdss[i].sdss_id),2)
	
	k2cat[ix].obs.gmag = sdss[i].gmag
	k2cat[ix].obs.e_gmag = sdss[i].e_gmag
	k2cat[ix].obs.rmag = sdss[i].rmag
	k2cat[ix].obs.e_rmag = sdss[i].e_rmag
	k2cat[ix].obs.imag = sdss[i].imag
	k2cat[ix].obs.e_imag = sdss[i].e_imag
    	
	k2cat[ix].obs.umag = sdss[i].umag
	k2cat[ix].obs.e_umag = sdss[i].e_umag
	k2cat[ix].obs.zmag = sdss[i].zmag
	k2cat[ix].obs.e_zmag = sdss[i].e_zmag
	
	k2cat[ix].obs.pmra = sdss[i].pmra
	k2cat[ix].obs.e_pmra = sdss[i].e_pmra
	k2cat[ix].obs.pmde = sdss[i].pmde
	k2cat[ix].obs.e_pmde = sdss[i].e_pmde
    	
	if (sdss[i].cl eq 3) then k2cat[ix].id.class = 'EXTENDED'

	progress,'cross-matching SDSS',i,nsdss
	count++
	ix++
	
endfor

progress,'cross-matching SDSS',i,nsdss
print,'used ',count,' of ',n_elements(sdss),' SDSS stars'
endif

print,''
print,'total number of unique objects:',long(ix)
k2cat = k2cat[0:ix-1]
print,''


; add in NOMAD proper motions
match_nomad=MATCH_2D(k2cat.obs.ra,k2cat.obs.de,nomad._raj2000,nomad._dej2000,range)

ad=0.
for i=0.,n_elements(k2cat)-1 do begin
	; if there is a proper motion already, don't override
	if (k2cat[i].obs.pmra ne 0. or k2cat[i].obs.pmde ne 0.) then continue
	if (match_nomad[i] eq -1) then continue
	pos = match_nomad[i]
	if (abs(k2cat[i].obs.jmag - nomad[pos].jmag) gt 0.005) then continue
	k2cat[i].obs.pmra = nomad[pos].pmra
	k2cat[i].obs.e_pmra = nomad[pos].e_pmra
	k2cat[i].obs.pmde = nomad[pos].pmde
	k2cat[i].obs.e_pmde = nomad[pos].e_pmde
	k2cat[i].id.nomad = 'NOMAD:'+nomad[pos].nomad1

	progress,'adding in NOMAD proper motions',i,ix
	ad++
endfor

progress,'adding in NOMAD proper motions',i,ix
print,'added ',ad,' proper motions'
print,''

crazymag=0.

; calculate Kepler magnitudes
for i=0.,ix-1 do begin

    gmag = 0.
    rmag = 0.
    imag = 0.
    kepmag = 0.

    kepmag = getkepmag(k2cat[i].obs.gmag,k2cat[i].obs.rmag,k2cat[i].obs.imag)

    ; if gri isn't available, use other colors
    if (kepmag eq 0.) then begin
	if (k2cat[i].obs.bmag ne 0. and k2cat[i].obs.vmag ne 0.) then begin
		grmag = 1.124*(k2cat[i].obs.bmag-k2cat[i].obs.vmag) - 0.252
    	    	;rimag = 1.040*(k2cat[i].obs.bmag-k2cat[i].obs.vmag) - 0.224
    	    	gmag = k2cat[i].obs.vmag + 0.634*(k2cat[i].obs.bmag-k2cat[i].obs.vmag) - 0.108
		rmag = gmag-grmag
		imag = 0.   ;rmag-rimag
		k2cat[i].id.kepflag='BV'
		kepmag = getkepmag(gmag,rmag,imag)
	endif else begin
		
		; if within color range of calibration, use polynomials by Howell et al. (2012)		
		if (k2cat[i].obs.jmag ne 0. and  k2cat[i].obs.hmag ne 0. and k2cat[i].obs.kmag ne 0.) then begin
	    	    x = k2cat[i].obs.jmag-k2cat[i].obs.kmag
	    	
		    if (k2cat[i].obs.jmag-k2cat[i].obs.hmag gt 0.75 and k2cat[i].obs.hmag-k2cat[i].obs.kmag gt 0.1 and $
		    x gt -0.2 and x lt 1.2) then begin
                        kepmag = 0.42443603 + 3.7937617*x - 2.3267277*x^2D + 1.4602553*x^3D + k2cat[i].obs.kmag
                    endif else begin
		    if (x gt -0.2 and x lt 1.0) then $
                        kepmag = (0.314377 + 3.85667*x + 3.176111*x^2D - 25.3126*x^3D + $
                        40.7221*x^4D - 19.2112*x^5D) + k2cat[i].obs.kmag
                    endelse
		    
		    k2cat[i].id.kepflag='JHK'
		    
		endif
		
		; otherwise use Howell et al. (2012) proxy based on J-band only
		if (kepmag eq 0. and k2cat[i].obs.jmag ne 0.) then begin
	   	    
		    if (k2cat[i].obs.jmag gt 10. and k2cat[i].obs.jmag lt 16.7) then kepmag = -398.04666+149.08127*k2cat[i].obs.jmag-$
		    	21.952130*k2cat[i].obs.jmag^2D + 1.5968619*k2cat[i].obs.jmag^3D - 0.057478947*k2cat[i].obs.jmag^4D +$
			0.00082033223*k2cat[i].obs.jmag^5D + k2cat[i].obs.jmag
			
		    if (k2cat[i].obs.jmag gt 16.7) then kepmag = 0.1918 + 0.08156*k2cat[i].obs.jmag + k2cat[i].obs.jmag
		    
		    ; this should be pretty rare ... if a star is brighter than J<10 and still has no kepmag, 
		    ; just assign J-band mag + a 1.7 mag offset estimated from Fig. 18 in Howell et al. (2012)
		    if (kepmag eq 0.) then begin
		    	crazymag++
			kepmag = k2cat[i].obs.jmag+1.7
		    endif
		    
		    ;grmag = 1.951*(k2cat[i].obs.jmag-k2cat[i].obs.hmag)+1.199*(k2cat[i].obs.hmag-k2cat[i].obs.kmag)-0.230
            	    ;rimag = 0.991*(k2cat[i].obs.jmag-k2cat[i].obs.hmag)+0.792*(k2cat[i].obs.hmag-k2cat[i].obs.kmag)-0.210
            	    ;gmagt = 1.379*grmag+1.702*rimag+0.518+k2cat[i].obs.jmag
            	    ;rmagt = -grmag+gmagt
            	    ;imagt = -rimag+rmagt
            	    ;gmag = gmagt
            	    ;rmag = rmagt
            	    ;imag = imagt
		    ;kepmag = getkepmag(gmag,rmag,imag)
		    
		    k2cat[i].id.kepflag='J'
		    
		endif    
	endelse
    endif else k2cat[i].id.kepflag='gri'
    k2cat[i].deriv.kp=kepmag
    
    progress,'calculating Kepmags',i,ix
endfor

progress,'calculating Kepmags',i,ix

print,crazymag

; removes entries with no photometry (no kepmags)
u = where(k2cat.deriv.kp ne 0.)
k2cat = k2cat[u]

print,n_elements(k2cat),' unique sources with Kepler magnitudes'

; assign IDs
ix = n_elements(k2cat)
epic = long(epicid)+indgen(ix,/long)
sorted = sort(k2cat.obs.de)
k2cat[sorted].id.epic = epic
sorted = sort(k2cat.id.epic)
k2cat = k2cat[sorted]

;stop
;save,file='output/'+filename+'.idl',k2cat

end
