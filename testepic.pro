; takes an EPIC catalog and checks coverage using DSS images; the input must be an ascii 
; output file generated by epic.pro (the same files which are delivered to the SOC and DMC)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; input:      path to input file (ascii output of epic.pro, either *soc.dat or *dmc.dat)
; inra/indec: input RA & DEC for test image; if not specified, the code will pick random 
;             coordinates near the center of the field
; imsize:     size of image in arcminutes (default=15', max=60')
; kplim:      faint limit to include when overlaying the catalog (default=30)
; syms:       scale factor for circle size around stars (default=10)
;
; examples: testepic,input='output/ra59.1_de18.7_r1.0.dmc.dat',inra=58.5,indec=18.3,imsize=10.
;           testepic,input='output/ra59.1_de18.7_r1.0.soc.dat',imsize=30.,kplim=16.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro testepic,input=input,inra=inra,indec=indec,imsize=imsize,kplim=kplim,syms=syms

if not keyword_set(input) then begin
    print,'error: input file needed.'
    stop
endif

if strmatch(input,'*soc*') then soc=1 else soc=0
if strmatch(input,'*dmc*') then dmc=1 else dmc=0

; expected number of columns
if (dmc) then ncol=66 else ncol=41

if not keyword_set(imsize) then ims=15. else ims=imsize
if not keyword_set(kplim) then kplim=30. else kplim=kplim
if not keyword_set(syms) then ss=10. else ss=syms

file=input

n=0.
s=''
openr,1,file
while not eof(1) do begin
    readf,1,s
    n++
endwhile
close,1

ra_in=dblarr(n)
de_in=dblarr(n)
kp_in=dblarr(n)
epic_in=lonarr(n)
str_in=strarr(n)

if (soc) then begin
    i=0.
    j=long(0)
    openr,1,file
    while not eof(1) do begin
    	readf,1,s	
    	t=strsplit(s,'|',/extract,/preserve_null)
    	ra_in[i]=double(t[0])*15D	; SOC RA is in hours!
    	de_in[i]=double(t[1])
    	kp_in[i]=double(t[14]) 
    	epic_in[i] = long(t[15])
    	str_in[i]=s 
    	if (n_elements(t) eq ncol) then j++        
    	i++
    endwhile
    close,1    
endif

if (dmc) then begin
    i=0.
    j=long(0)
    openr,1,file
    while not eof(1) do begin
    	readf,1,s	
    	t=strsplit(s,'|',/extract,/preserve_null)
    	ra_in[i]=double(t[9])
    	de_in[i]=double(t[10])
    	kp_in[i]=double(t[45]) 
    	epic_in[0] = long(t[0])
    	str_in[i]=s   
    	if (n_elements(t) eq ncol) then j++   
    	i++
    endwhile
    close,1    
endif

; some sanity checks
print,'------------------------------------'
print,'Sanity checks:'
print,'------------------------------------'
print,'total number of catalog entries:',n_elements(ra_in)
print,'number of entries with Kepmag>0:',n_elements(where(kp_in gt 0.))
print,'number of entries with',ncol,'  colums:',j
print,'min Kp mag:',min(kp_in,pos),'  at RA=',ra_in[pos],'  DEC=',de_in[pos]
print,'max Kp mag:',max(kp_in)

an=''
read,an,prompt='press enter to continue'

; use input coordinates when provided; otherwise, pick at random
if keyword_set(inra) and keyword_set(indec) then begin
    centera=inra
    centerd=indec
    ran=0.
endif else begin
    centera=mean(ra_in)
    centerd=mean(de_in)
    ran=((max(ra_in)-min(ra_in))/2.)*0.3
endelse

; no reason to loop if there is only one set of input coordinates
if keyword_set(inra) and keyword_set(indec) then nit=0 else nit=1e3

for i=0.,nit do begin

    cenra = centera+randomn(seed)*ran
    cende = centerd+randomn(seed)*ran

    print,'grabbing image centered on:',cenra,cende

    QueryDSS,[cenra,cende], Im, Hdr,imsize=ims
    spawn,'rm test.fits'
    MWRFITS, Im, 'test.fits', hdr 

    dis = sqrt( (ra_in-cenra)^2D + (de_in-cende)^2D )
    u = where(dis lt ims/60. and kp_in lt kplim)
    ra = ra_in[u]
    de = de_in[u]
    kp = kp_in[u]
    str = str_in[u]
    maxb = min(kp)
    sizes = replicate(ss,n_elements(ra))
    sizes = sizes*(maxb/kp)^2.

    openw,1,'ds9.reg'
    printf,1,'global color=red, font="helvetica 6 normal roman"'
    printf,1,'fk5'
    for q=0.,n_elements(ra)-1 do begin
	printf,1,'circle   ',ra[q],de[q],sizes[q],'" # text={'+strtrim(string(kp[q],format='(d5.1)'),2)+'}'
    endfor

    close,1

    spawn,'ds9 test.fits -regions ds9.reg -width 750 -height 750 -scale mode 99.5 -zoom to fit -wcs degrees';,unit=dummy
    
    ; some options for saving the plots as png files
    ; filename='t.png'
    ;spawn,'ds9 test.fits -regions ds9.reg -width 750 -height 750 -zoom to fit -wcs degrees -scale mode 99.5 -invert -saveimage png '+filename,unit=dummy
    ;spawn,'ds9 test.fits -regions ds9.reg -width 750 -height 750 -zoom to fit -wcs degrees - saveimage png C4_pleiades_kplt12.png ',unit=dummy

    if (nit gt 0.) then begin
    	an=''
    	read,an,prompt='press enter for next image'
	endif

endfor

end
