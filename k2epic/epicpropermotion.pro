; takes an EPIC catalog and checks for high proper motion stars using DSS images
; each images shows: original positions of all stars in the EPIC (red circles), 
; EPIC position of the target that's being inspected (blue square), proper motion 
; corrected position of the target that's being inspected (green circle); 
; if the EPIC proper motions are correct, the green circle should overlap with a star of 
; roughly the correct brightness
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; REQUIRED INPUT:
; epic:          path to EPIC SOC delivery file (*soc.mrg)
; lctargetfile:  path to LC target file
; sctargetfile:  path to SC target file
;
; OPTIONAL INPUT:
; excludefile:   path to a file containing EPIC IDs that should be excluded (i.e. high
;                proper motion stars delivered by Tom based on proposer input)
; pmlim:         proper motion limit in pixels (default=1)
; imsize:        size of image for visual checking (default=3 arcminutes)
; outfile:       output file name (default=propermotions_Cx_epic.dat)
;
; OUTPUT:	- output file containing EPIC ID, proper motion in RA (as/yr) and 
;                 proper motion in DEC (as/yr)
;               - screenshot for each visually checked image (in screenshots/ directory)
; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro epicpropermotion,epic=epic,lctargetfile=lctargetfile,sctargetfile=sctargetfile,$
pmlim=pmlim,excludefile=excludefile,imsize=imsize,outfile=outfile

; due to long paths, defining the input here may be easier
;epic='~/science/K2/EPIC/deliveries/14252_01_epic_c3_soc/d14252_01_epic_c3_soc_comprehensive.mrg'
;lctargetfile='~/science/K2/EPIC/C23/C3/c3_nov2014/go_lc_campaign3_oct2014.txt'
;sctargetfile='~/science/K2/EPIC/C23/C3/c3_nov2014/go_sc_campaign3_oct2014.txt'
;excludefile='~/science/K2/EPIC/C23/C3/c3_nov2014/propermotion_c3.txt'


if not keyword_set(epic) then begin 
	print,'an EPIC input catalog must be provided.'
	stop
endif

if not keyword_set(lctargetfile) then begin 
	print,'a LC target file must be provided.'
	stop
endif

if not keyword_set(sctargetfile) then begin 
	print,'a SC target file must be provided.'
	stop
endif

if not keyword_set(pmlim) then pmlim=1.
if not keyword_set(imsize) then imsize=3.
if not keyword_set(outfile) then outfile='propermotions_Cx_epic.dat'

; fudge the number in the following line if circles are too big or too small
ss=10.

; read the EPIC SOC delivery
n=0.
s=''
openr,1,epic
while not eof(1) do begin
    readf,1,s
    n++
endwhile
close,1

pmra_in=dblarr(n)
pmde_in=dblarr(n)
ra_in=dblarr(n)
de_in=dblarr(n)
kp_in=dblarr(n)
epic_in=lonarr(n)
str_in=strarr(n)

i=0.
openr,1,epic
while not eof(1) do begin
    readf,1,s	
    t=strsplit(s,'|',/extract,/preserve_null)
    ra_in[i]=double(t[0])*15D
    de_in[i]=double(t[1])
    pmra_in[i]=double(t[2])
    pmde_in[i]=double(t[3])
    kp_in[i]=double(t[14]) 
    epic_in[i] = long(t[15])
    str_in[i]=s      
    i++
endwhile
close,1
ra_all=ra_in
de_all=de_in
kp_all=kp_in


; only consider targets that went into target management
readcol,lctargetfile,lctargets,format='L'
readcol,sctargetfile,sctargets,format='L'
alltargets=[lctargets,sctargets]
;alltargets=alltargets[unique(alltargets)]
alltargets=unique(alltargets,/sort)

; exclude proper motions provided by Tom (user-provided)
if keyword_set(excludefile) then begin
	readcol,excludefile,epicpm,format='L'
	match2,alltargets,epicpm,ix,iy
	u=where(ix eq -1)
	alltargets=alltargets[u]
endif

match,epic_in,alltargets,ix,iy
epic_in=epic_in[ix]
ra_in=ra_in[ix]
de_in=de_in[ix]
pmra_in=pmra_in[ix]
pmde_in=pmde_in[ix]
kp_in=kp_in[ix]
str_in=str_in[ix]

; calculate proper motion since JD2000 in units of Kepler pixels
pmo=sqrt(pmde_in^2D + (pmra_in^2D * cos(de_in)*cos(de_in)))   
thisyear=float((strsplit(systime(),/extract))[4])
pmo=pmo*(thisyear-2000.)/4.

s = where(pmo gt pmlim)
print,strtrim(string(n_elements(s)),2),' stars with proper motion above',pmlim,$
' pixels selected. Hit Enter to continue'
an=''
read,an

; loop over those stars and check whether EPIC proper motions make sense
openw,5,outfile
for i=0.,n_elements(s)-1 do begin

	; get the image
	cenra = ra_in[s[i]]
	cende = de_in[s[i]]
	cenkp = kp_in[s[i]]
	QueryDSS,[cenra,cende], Im, Hdr, ImSize=imsize
	spawn,'rm test.fits'
	MWRFITS, Im, 'test.fits', hdr 

	dateobs=double(strmid(fxpar(Hdr,'DATE-OBS'),0,4))
	
	; for some reason this is what needs to be done on Dan's Mac to get dateobs
	;temp=fxpar(Hdr,'DATE-OBS')
	;t=strsplit(temp,'/',/extract)
	;dateobs=float(strtrim(string(19),2)+(t[2]))
	
	; correct proper motions
	tstep=2000.-dateobs	
	cenracor = ra_in[s[i]]-tstep*(pmra_in[s[i]]/3600.)
	cendecor = de_in[s[i]]-tstep*(pmde_in[s[i]]/3600.)

	print,'---'
	print,'EPIC:',epic_in[s[i]]
	print,'RA,DEC,Kp:',cenra,cende,cenkp
	print,'PMRA,PMDEC:',pmra_in[s[i]],pmde_in[s[i]]
	print,'correcting',tstep,'years of proper motion'
	print,'corrected RA & DEC:',cenracor,cendecor

	dis = sqrt( (ra_all-cenra)^2D + (de_all-cende)^2D )
	u = where(dis lt 0.3)
	ra = ra_all[u]
	de = de_all[u]
	kp = kp_all[u]
	maxb = min(kp)
	sizes = replicate(ss,n_elements(ra))
	sizes = sizes*(maxb/kp)^2.

	openw,1,'ds9.reg'
	printf,1,'global color=red, font="helvetica 6 normal roman"'
	printf,1,'fk5'
	for q=0.,n_elements(ra)-1 do begin
		printf,1,'circle   ',ra[q],de[q],sizes[q],'" # text={'+strtrim(string(kp[q],format='(d5.1)'),2)+'}'
	endfor

	sizes = ss*(maxb/cenkp)^2.

	printf,1,'global color=blue, font="helvetica 6 normal roman"'
	printf,1,'fk5'
	printf,1,'box   ',ra_in[s[i]],de_in[s[i]],' ',strtrim(string(2.*sizes),2)+'" ',strtrim(string(2.*sizes),2)+'" # text={'+strtrim(string(cenkp,format='(d5.1)'),2)+'}',$
	format='(A,D,D,A,A,A,A)'

	printf,1,'global color=green, font="helvetica 6 normal roman"'
	printf,1,'fk5'
	printf,1,'circle   ',cenracor,cendecor,sizes,'" # text={'+strtrim(string(cenkp,format='(d5.1)'),2)+'}'

	close,1

	filename='screenshots/'+strtrim(string(epic_in[s[i]]),2)+'.png'
	print,filename

    	spawn,'ds9 test.fits -regions ds9.reg -width 750 -height 750 -zoom to fit -wcs degrees -saveimage png '+filename,unit=dummy,pid=pid
	wait,2.5

	print,'log this as correct proper motion (1-yes,0-no)'
	an=0
	read,an
	
	spawn,'kill '+strtrim(string(pid),2)

	if an eq 1 then begin
		printf,5,epic_in[s[i]],pmra_in[s[i]],pmde_in[s[i]]
	endif

endfor
close,/all
stop

end
