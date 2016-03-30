;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; writes EPIC catalog to ASCII format
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

@fixstring.pro

pro makeascii,k2cat,filename=filename,meanra=meanra,meande=meande,maxrad=maxrad

device,set_font='Helvetica Bold',/tt_font
;tvlct,0,0,0,255
;tvlct,255,255,255,0

;window,xsize=900,ysize=1000
!p.charsize=3.0
!p.multi=[0,2,3]
u = where(strmatch(k2cat.id.hip,'*HIP*'))
if u[0] ne -1 then plot,k2cat[u].obs.ra,k2cat[u].obs.de,psym=3,$
xrange=[meanra-maxrad,meanra+maxrad],yrange=[meande-maxrad,meande+maxrad],/xs,/ys,$
title='Hipparcos',xtitle='RA (deg)',ytitle='DEC (deg)',FONT=1 else plot,[0],[0]

u = where(strmatch(k2cat.id.tycho,'*TYC*'))
if u[0] ne -1 then plot,k2cat[u].obs.ra,k2cat[u].obs.de,psym=3,$
xrange=[meanra-maxrad,meanra+maxrad],yrange=[meande-maxrad,meande+maxrad],/xs,/ys,$
title='Tycho-2',xtitle='RA (deg)',ytitle='DEC (deg)',FONT=1 else plot,[0],[0]

u = where(strmatch(k2cat.id.ucac,'*UCAC*'))
if u[0] ne -1 then plot,k2cat[u].obs.ra,k2cat[u].obs.de,psym=3,$
xrange=[meanra-maxrad,meanra+maxrad],yrange=[meande-maxrad,meande+maxrad],/xs,/ys,$
title='UCAC-4',xtitle='RA (deg)',ytitle='DEC (deg)',FONT=1 else plot,[0],[0]

u = where(strmatch(k2cat.id.mass,'*MASS*'))
if u[0] ne -1 then plot,k2cat[u].obs.ra,k2cat[u].obs.de,psym=3,$
xrange=[meanra-maxrad,meanra+maxrad],yrange=[meande-maxrad,meande+maxrad],/xs,/ys,$
title='2MASS',xtitle='RA (deg)',ytitle='DEC (deg)',FONT=1 else plot,[0],[0]

u = where(strmatch(k2cat.id.sdss,'*SDSS*'))
if u[0] ne -1 then plot,k2cat[u].obs.ra,k2cat[u].obs.de,psym=3,$
xrange=[meanra-maxrad,meanra+maxrad],yrange=[meande-maxrad,meande+maxrad],/xs,/ys,$
title='SDSS DR9',xtitle='RA (deg)',ytitle='DEC (deg)',FONT=1 else plot,[0],[0]

u = where(strmatch(k2cat.id.nomad,'*NOMAD*'))
if u[0] ne -1 then plot,k2cat[u].obs.ra,k2cat[u].obs.de,psym=3,$
xrange=[meanra-maxrad,meanra+maxrad],yrange=[meande-maxrad,meande+maxrad],/xs,/ys,$
title='NOMAD proper motions',xtitle='RA (deg)',ytitle='DEC (deg)',FONT=1 else plot,[0],[0]
!p.multi=0

;saveimage,'gifs/'+filename+'.output.gif'

u = where(strmatch(k2cat.id.class,'STAR'))
u2 = where(strmatch(k2cat.id.class,'EXTENDED'))
;k2cat[u2].id.class = 'EXTENDED'

glxid = replicate(0,n_elements(k2cat))
glxid[u2] = 1

;ix = n_elements(k2cat)
;epic = long(201000001)+indgen(ix,/long)
;k2cat.id.epic = epic

; stellar properties have not yet been validated, so set them all to zero for this delivery
k2cat[*].deriv.teff[3] = 0.
k2cat[*].deriv.logg[3] = 0.
k2cat[*].deriv.feh[3] = 0.
k2cat[*].deriv.rad[3] = 0.
k2cat[*].deriv.mass[3] = 0.
k2cat[*].deriv.rho[3] = 0.
k2cat[*].deriv.lum[3] = 0.
k2cat[*].deriv.dis[3] = 0.
k2cat[*].deriv.ebv = 0.

print,'final catalog has',n_elements(k2cat),'total sources'
print,'final catalog has',n_elements(u),'sources classified as stars'
print,'final catalog has',n_elements(u2),'sources classified as galaxies'

print,'writing ASCII output'

openw,3,'output/'+filename+'.dmc.dat'
openw,2,'output/'+filename+'.soc.dat'
openw,1,'output/'+filename+'.dmc.dat.readme'
;openw,2,'delivery/1435_01_epic_field1_soc/d1435_01_epic_field1_soc.mrg'
;openw,3,'delivery/1435_02_epic_field1_dmc/d1435_02_epic_field1_dmc.mrg'

printf,1,'# K2 Input Catalog'
printf,1,'# file generated on '+systime(0)
printf,1,'---------------------------------------------------------------------------------'
printf,1,'### Identifiers & Flags'
printf,1,'#1  ID          I10           ... K2 Input Catalog Identifier'
printf,1,'#2  HIP         I8            ... Hipparcos Identifier'
printf,1,'#3  TYC         A15           ... Tycho2 Identifier'
printf,1,'#4  UCAC        A15           ... UCAC4 Identifier'
printf,1,'#5  2MASS       A20           ... 2MASS Identifier'
printf,1,'#6  SDSS        A20           ... SDSS DR9 Identifier'
printf,1,'#7  Objtype     A10           ... Object Type'
printf,1,'#8  Kepflag     A5            ... Kepler Magnitude Flag'
printf,1,'#9  StpropFlag  A5            ... Stellar Properties Flag'
printf,1,'---------------------------------------------------------------------------------'
printf,1,'### Observables'
printf,1,'#10 RA          D10.6         ... Right Ascension JD2000 (Deg)'
printf,1,'#11 DEC         D10.6         ... Declination JD2000 (Deg)'
printf,1,'#12 pmRA        D10.3         ... Proper motion in RA (mas/yr)'
printf,1,'#13 e_pmRA      D10.3         ... Uncertainty in Proper motion in RA (mas/yr)'
printf,1,'#14 pmDEC       D10.3         ... Proper motion in DEC (mas/yr)'
printf,1,'#15 e_pmDEC     D10.3         ... Uncertainty in Proper motion in DEC (mas/yr)'
printf,1,'#16 plx         D10.3         ... Hipparcos parallax (mas)'
printf,1,'#17 e_plx       D10.3         ... Uncertainty in parallax (mas)'
printf,1,'#18 Bmag        D6.3          ... Johnson B Magnitude (mag)'
printf,1,'#19 e_Bmag      D6.3          ... Uncertainty in Johnson B Magnitude (mag)'
printf,1,'#20 Vmag        D6.3          ... Johnson V Magnitude (mag)'
printf,1,'#21 e_Vmag      D6.3          ... Uncertainty in Johnson V Magnitude (mag)'
printf,1,'#22 umag        D6.3          ... Sloan u Magnitude (mag)'
printf,1,'#23 e_umag      D6.3          ... Uncertainty in Sloan u Magnitude (mag)'
printf,1,'#24 gmag        D6.3          ... Sloan g Magnitude (mag)'
printf,1,'#25 e_gmag      D6.3          ... Uncertainty in Sloan g Magnitude (mag)'
printf,1,'#26 rmag        D6.3          ... Sloan r Magnitude (mag)'
printf,1,'#27 e_rmag      D6.3          ... Uncertainty in Sloan r Magnitude (mag)'
printf,1,'#28 imag        D6.3          ... Sloan i Magnitude (mag)'
printf,1,'#29 e_imag      D6.3          ... Uncertainty in Sloan i Magnitude (mag)'
printf,1,'#30 zmag        D6.3          ... Sloan z Magnitude (mag)'
printf,1,'#31 e_zmag      D6.3          ... Uncertainty in Sloan z Magnitude (mag)'
printf,1,'#32 Jmag        D6.3          ... 2MASS J Magnitude (mag)'
printf,1,'#33 e_Jmag      D6.3          ... Uncertainty in 2MASS J Magnitude (mag)'
printf,1,'#34 Hmag        D6.3          ... 2MASS H Magnitude (mag)'
printf,1,'#35 e_Hmag      D6.3          ... Uncertainty in 2MASS H Magnitude (mag)'
printf,1,'#36 Kmag        D6.3          ... 2MASS K Magnitude (mag)'
printf,1,'#37 e_Kmag      D6.3          ... Uncertainty in 2MASS K Magnitude (mag)'
printf,1,'#38 W1mag       D6.3          ... WISE W1 Magnitude (mag)'
printf,1,'#39 e_W1mag     D6.3          ... Uncertainty in WISE W1 Magnitude (mag)'
printf,1,'#40 W2mag       D6.3          ... WISE W2 Magnitude (mag)'
printf,1,'#41 e_W2mag     D6.3          ... Uncertainty in WISE W2 Magnitude (mag)'
printf,1,'#42 W3mag       D6.3          ... WISE W3 Magnitude (mag)'
printf,1,'#43 e_W3mag     D6.3          ... Uncertainty in WISE W3 Magnitude (mag)'
printf,1,'#44 W4mag       D6.3          ... WISE W4 Magnitude (mag)'
printf,1,'#45 e_W4mag     D6.3          ... Uncertainty in WISE W4 Magnitude (mag)'
printf,1,'---------------------------------------------------------------------------------'
printf,1,'### Derived Properties'
printf,1,'#46 Kp          D6.3          ... Kepler Magnitude (mag)'
printf,1,'#47 Teff        D6.0          ... Effective Temperature (K)'
printf,1,'#48 ep_Teff     D6.0          ... Upper Uncertainty in Effective Temperature (K)'
printf,1,'#49 em_Teff     D6.0          ... Lower Uncertainty in Effective Temperature (K)'
printf,1,'#50 logg        D6.3          ... log Surface Gravity (cgs)'
printf,1,'#51 ep_logg     D6.3          ... Upper Uncertainty in log Surface Gravity (cgs)'
printf,1,'#52 em_logg     D6.3          ... Lower Uncertainty in log Surface Gravity (cgs)'
printf,1,'#53 [Fe/H]      D6.3          ... Metallicity (dex)'
printf,1,'#54 ep_[Fe/H]   D6.3          ... Upper Uncertainty in Metallicity (dex)'
printf,1,'#55 em_[Fe/H]   D6.3          ... Lower Uncertainty in Metallicity (dex)'
printf,1,'#56 Rad         D8.3          ... Stellar Radius (solar units)'
printf,1,'#57 ep_Rad      D8.3          ... Upper Uncertainty in Stellar Radius (solar units)'
printf,1,'#58 em_Rad      D8.3          ... Lower Uncertainty in Stellar Radius (solar units)'
printf,1,'#59 Mass        D8.3          ... Stellar Mass (solar units)'
printf,1,'#60 ep_Mass     D8.3          ... Upper Uncertainty in Stellar Mass (solar units)'
printf,1,'#61 em_Mass     D8.3          ... Lower Uncertainty in Stellar Mass (solar units)'
printf,1,'#62 rho         E10.3         ... Stellar Density (solar units)'
printf,1,'#63 ep_rho      E10.3         ... Upper Uncertainty in Stellar Density (solar units)'
printf,1,'#64 em_rho      E10.3         ... Lower Uncertainty in Stellar Density (solar units)'
printf,1,'#65 Lum         E10.3         ... Stellar Luminosity (solar units)'
printf,1,'#66 ep_Lum      E10.3         ... Upper Uncertainty in Stellar Luminosity (solar units)'
printf,1,'#67 em_Lum      E10.3         ... Lower Uncertainty in Stellar Luminosity (solar units)'
printf,1,'#68 d           E9.3          ... Distance (pc)'
printf,1,'#69 ep_d        E9.3          ... Upper Uncertainty in Distance (parsec)'
printf,1,'#70 em_d        E9.3          ... Lower Uncertainty in Distance (parsec)'
printf,1,'#71 E(B-V)      D7.4          ... Extinction (mag)'
printf,1,'#72 ep_E(B-V)   D6.4          ... Upper Uncertainty in Extinction (mag)'
printf,1,'#73 em_E(B-V)   D6.4          ... Upper Uncertainty in Extinction (mag)'
printf,1,'---------------------------------------------------------------------------------'
printf,1,'### Additional Information'
printf,1,'#74 NOMAD       A20           ... NOMAD1 Identifier (for sources with supplemented NOMAD1 proper motions)'
printf,1,'#75 Mflg        A20           ... 2MASS Photometry Flags (Qflg-Rflg-Bflg-Cflg-Xflg-Aflg)'
printf,1,'#76 prox        D6.3          ... 2MASS nearest neighbor distance (arcsec)'
printf,1,'---------------------------------------------------------------------------------'
close,1

;printf,1,'   EPIC    |    HIP   |       TYC       |      UCAC       |       2MASS          |         SDSS         |   Objtype  | Kepflag | StparFlag |      RA    |   DEC      |  Bmag   |  e_Bmag |  Vmag   |  e_Vmag |  umag   |  e_umag |   gmag  |  e_gmag |  rmag   |  e_rmag |  imag   |  e_imag |  zmag   |  e_zmag |  jmag   |  e_jmag |  hmag   |  e_hmag |  kmag   |  e_kmag |    Kp   |    Teff  |   e_Teff |   logg   |  e_logg  |  [Fe/H]  |  e_[Fe/H] |     Rad   |    e_Rad   |    Mass    |    e_Mass  |   rho      |  e_rho     |    Lum     |    e_Lum   |      d     |    e_d     |   E(B-V)'
;printf,1,'ID|HIP|TYC|UCAC|2MASS|SDSS|NOMAD|Objtype|Kepflag|StpropFlag|RA|DEC|pmRA|e_pmRA|pmDEC|e_pmDEC|'+$
;'plx|e_plx|Bmag|e_Bmag|Vmag|e_Vmag|umag|e_umag|gmag|e_gmag|rmag|e_rmag|imag|e_imag|zmag|e_zmag|jmag'+$
;'|e_jmag|hmag|e_hmag|kmag|e_kmag|Mflg|prox|w1mag|e_w1mag|w2mag|e_w2mag|w3mag|e_w3mag|w4mag|e_w4mag|Kp|Teff|e_Teff|'+$
;'logg|e_logg|[Fe/H]|e_[Fe/H]|Rad|e_Rad|Mass|e_Mass|rho|e_rho|Lum|e_Lum|d|e_d|E(B-V)'

idstr = '(I10,A3,A8,A3,A15,A3,A15,A3,A20,A3,A20,A3,A10,A3,A7,A3,A9,A3,'
radec = 'd10.6,A3,d10.6,A3,'
mags = 'd7.3,A3,d7.3,A3,d7.3,A3,d7.3,A3,d7.3,A3,d7.3,A3,d7.3,A3,d7.3,A3,d7.3,A3,d7.3,A3,'+$
	   'd7.3,A3,d7.3,A3,A30,A3,d7.3,A3,d7.3,A3,d7.3,A3,d7.3,A3,d7.3,A3,d7.3,A3,d7.3,A3,d7.3,A3,d7.3,A3,'
derv = 'd7.3,A3,d8.0,A3,d8.0,A3,d8.3,A3,d8.3,A3,d8.3,A3,d8.3,A3,E10.3,A3,E10.3,A3,E10.3,A3,E10.3,A3,E10.3,A3,E10.3,A3,E10.3,A3,E10.3,A3,E10.3,A3,E10.3,A3,d8.3)'

for i=0.,n_elements(k2cat)-1 do begin

t = strsplit(k2cat[i].id.hip,':',/extract)
t2 = strsplit(k2cat[i].id.tycho,':',/extract)
t3 = strsplit(k2cat[i].id.ucac,':',/extract)
t4 = strsplit(k2cat[i].id.mass,':',/extract)
t5 = strsplit(k2cat[i].id.sdss,':',/extract)
t6 = strsplit(k2cat[i].id.nomad,':',/extract)

if n_elements(t) eq 1 then t = ['',''] else t[1] = long(t[1])
if n_elements(t2) eq 1 then t2 = ['','']
if n_elements(t3) eq 1 then t3 = ['','']
if n_elements(t4) eq 1 then t4 = ['','']
if n_elements(t5) eq 1 then t5 = ['','']
if n_elements(t6) eq 1 then t6 = ['','']

teffe = 0.
logge = 0.
fehe = 0.
rade = 0.
masse = 0.
rhoe = 0.
lume = 0.
de = 0.

if (k2cat[i].obs.bmag eq 0.) then p1=''+'|'+''+'|' else p1 = strtrim(string(k2cat[i].obs.bmag,format='(d6.3)'),2)+'|'+strtrim(string(k2cat[i].obs.e_bmag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.vmag eq 0.) then p2=''+'|'+''+'|' else p2 = strtrim(string(k2cat[i].obs.vmag,format='(d6.3)'),2)+'|'+strtrim(string(k2cat[i].obs.e_vmag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.umag eq 0.) then p3=''+'|'+''+'|' else p3 = strtrim(string(k2cat[i].obs.umag,format='(d6.3)'),2)+'|'+strtrim(string(k2cat[i].obs.e_umag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.gmag eq 0.) then p4=''+'|'+''+'|' else p4 = strtrim(string(k2cat[i].obs.gmag,format='(d6.3)'),2)+'|'+strtrim(string(k2cat[i].obs.e_gmag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.rmag eq 0.) then p5=''+'|'+''+'|' else p5 = strtrim(string(k2cat[i].obs.rmag,format='(d6.3)'),2)+'|'+strtrim(string(k2cat[i].obs.e_rmag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.imag eq 0.) then p6=''+'|'+''+'|' else p6 = strtrim(string(k2cat[i].obs.imag,format='(d6.3)'),2)+'|'+strtrim(string(k2cat[i].obs.e_imag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.zmag eq 0.) then p7=''+'|'+''+'|' else p7 = strtrim(string(k2cat[i].obs.zmag,format='(d6.3)'),2)+'|'+strtrim(string(k2cat[i].obs.e_zmag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.jmag eq 0.) then p8=''+'|'+''+'|' else p8 = strtrim(string(k2cat[i].obs.jmag,format='(d6.3)'),2)+'|'+strtrim(string(k2cat[i].obs.e_jmag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.hmag eq 0.) then p9=''+'|'+''+'|' else p9 = strtrim(string(k2cat[i].obs.hmag,format='(d6.3)'),2)+'|'+strtrim(string(k2cat[i].obs.e_hmag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.kmag eq 0.) then p10=''+'|'+''+'|' else p10 = strtrim(string(k2cat[i].obs.kmag,format='(d6.3)'),2)+'|'+strtrim(string(k2cat[i].obs.e_kmag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.pmra eq 0. and k2cat[i].obs.e_pmra eq 0.) then p11=''+'|'+''+'|' else p11 = strtrim(string(k2cat[i].obs.pmra,format='(d10.3)'),2)+'|'+strtrim(string(k2cat[i].obs.e_pmra,format='(d10.3)'),2)+'|'
if (k2cat[i].obs.pmde eq 0. and k2cat[i].obs.e_pmde eq 0.) then p12=''+'|'+''+'|' else p12 = strtrim(string(k2cat[i].obs.pmde,format='(d10.3)'),2)+'|'+strtrim(string(k2cat[i].obs.e_pmde,format='(d10.3)'),2)+'|'
if (k2cat[i].obs.plx eq 0. and k2cat[i].obs.e_plx eq 0.) then p13=''+'|'+''+'|' else p13 = strtrim(string(k2cat[i].obs.plx,format='(d10.3)'),2)+'|'+strtrim(string(k2cat[i].obs.e_plx,format='(d10.3)'),2)+'|'
; additional 2MASS columns added after C1
if (strlen(k2cat[i].obs.mflg) lt 2) then m0=''+'|' else m0 = k2cat[i].obs.mflg+'|'
if (k2cat[i].obs.prox eq 0.) then m1='' else m1 = strtrim(string(k2cat[i].obs.prox,format='(d8.3)'),2);+'|'



; SOC delivery requires parallaxes/proper motions in arcsecs
if (k2cat[i].obs.bmag eq 0.) then p1s=''+'|' else p1s = strtrim(string(k2cat[i].obs.bmag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.vmag eq 0.) then p2s=''+'|' else p2s = strtrim(string(k2cat[i].obs.vmag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.umag eq 0.) then p3s=''+'|' else p3s = strtrim(string(k2cat[i].obs.umag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.gmag eq 0.) then p4s=''+'|' else p4s = strtrim(string(k2cat[i].obs.gmag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.rmag eq 0.) then p5s=''+'|' else p5s = strtrim(string(k2cat[i].obs.rmag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.imag eq 0.) then p6s=''+'|' else p6s = strtrim(string(k2cat[i].obs.imag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.zmag eq 0.) then p7s=''+'|' else p7s = strtrim(string(k2cat[i].obs.zmag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.jmag eq 0.) then p8s=''+'|' else p8s = strtrim(string(k2cat[i].obs.jmag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.hmag eq 0.) then p9s=''+'|' else p9s = strtrim(string(k2cat[i].obs.hmag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.kmag eq 0.) then p10s=''+'|' else p10s = strtrim(string(k2cat[i].obs.kmag,format='(d6.3)'),2)+'|'
if (k2cat[i].obs.pmra eq 0. and k2cat[i].obs.e_pmra eq 0.) then p11s=''+'|' else p11s = strtrim(string(k2cat[i].obs.pmra/1000D,format='(d12.6)'),2)+'|'
if (k2cat[i].obs.pmde eq 0. and k2cat[i].obs.e_pmde eq 0.) then p12s=''+'|' else p12s = strtrim(string(k2cat[i].obs.pmde/1000D,format='(d12.6)'),2)+'|'
if (k2cat[i].obs.plx eq 0. and k2cat[i].obs.e_plx eq 0.) then p13s=''+'|' else p13s = strtrim(string(k2cat[i].obs.plx/1000D,format='(d12.6)'),2)+'|'

; stellar properties are delivered separately, so these are all left empty
stpropstr='|||||||||||||||||||||||||||'

;if (k2cat[i].deriv.teff[3] eq 0.) then d1=''+'|'+''+'|' else d1 = strtrim(string(k2cat[i].deriv.teff[3],format='(d6.0)'),2)+'|'+strtrim(string(teffe,format='(d6.0)'),2)+'|'
;if (k2cat[i].deriv.logg[3] eq 0.) then d2=''+'|'+''+'|' else d2 = strtrim(string(k2cat[i].deriv.logg[3],format='(d6.3)'),2)+'|'+strtrim(string(logge,format='(d6.3)'),2)+'|'
;if (k2cat[i].deriv.feh[3] eq 0.) then d3=''+'|'+''+'|' else d3 = strtrim(string(k2cat[i].deriv.feh[3],format='(d6.3)'),2)+'|'+strtrim(string(fehe,format='(d6.3)'),2)+'|'
;if (k2cat[i].deriv.rad[3] eq 0.) then d4=''+'|'+''+'|' else d4 = strtrim(string(k2cat[i].deriv.rad[3],format='(d8.3)'),2)+'|'+strtrim(string(rade,format='(d8.3)'),2)+'|'
;if (k2cat[i].deriv.mass[3] eq 0.) then d5=''+'|'+''+'|' else d5 = strtrim(string(k2cat[i].deriv.mass[3],format='(d8.3)'),2)+'|'+strtrim(string(masse,format='(d8.3)'),2)+'|'
;if (k2cat[i].deriv.rho[3] eq 0.) then d6=''+'|'+''+'|' else d6 = strtrim(string(k2cat[i].deriv.rho[3],format='(e10.3)'),2)+'|'+strtrim(string(rhoe,format='(e10.3)'),2)+'|'
;if (k2cat[i].deriv.lum[3] eq 0.) then d7=''+'|'+''+'|' else d7 = strtrim(string(k2cat[i].deriv.lum[3],format='(e10.3)'),2)+'|'+strtrim(string(lume,format='(e10.3)'),2)+'|'
;if (k2cat[i].deriv.dis[3] eq 0.) then d8=''+'|'+''+'|' else d8 = strtrim(string(k2cat[i].deriv.dis[3],format='(d8.1)'),2)+'|'+strtrim(string(de,format='(d8.1)'),2)+'|'
;if (k2cat[i].deriv.ebv eq 0.) then d9=''+'|' else d9 = strtrim(string(k2cat[i].deriv.ebv,format='(d6.3)'),2)+'|'

; DMC delivery
str = strtrim(string(k2cat[i].id.epic),2)+'|'+strtrim(fixstring(t[1]),2)+$
'|'+strtrim(fixstring(t2[1]),2)+'|'+strtrim(fixstring(t3[1]),2)+'|'+$
strtrim(fixstring(t4[1]),2)+'|'+strtrim(fixstring(t5[1]),2)+'|'+$
k2cat[i].id.class+'|'+k2cat[i].id.kepflag+'|'+k2cat[i].id.stparflag+'|'+$
strtrim(string(k2cat[i].obs.ra,format='(d10.6)'),2)+'|'+$
strtrim(string(k2cat[i].obs.de,format='(d10.6)'),2)+'|'+$
p11+p12+p13+p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+$
; the next 4 lines are placeholders for WISE
''+'|'+''+'|'+$
''+'|'+''+'|'+$
''+'|'+''+'|'+$
''+'|'+''+'|'+$
strtrim(string(k2cat[i].deriv.kp,format='(d6.3)'),2)+'|'+$
stpropstr+strtrim(t6[1],2)+'|'+m0+m1

; SOC delivery requires RA in hours
; SOC delivery, RA, DEC, Kepmag, KeplerID only
str2 = strtrim(string(k2cat[i].obs.ra/15D,format='(d10.6)'),2)+'|'+strtrim(string(k2cat[i].obs.de,format='(d10.6)'),2)+'|||||||||||||'+$
strtrim(string(k2cat[i].deriv.kp,format='(d6.3)'),2)+'|'+strtrim(string(k2cat[i].id.epic),2)+'|||||||||||||||||||||||||'

; SOC delivery
str3 = strtrim(string(k2cat[i].obs.ra/15D,format='(d10.6)'),2)+'|'+strtrim(string(k2cat[i].obs.de,format='(d10.6)'),2)+'|'+$
p11s+p12s+p3s+p4s+p5s+p6s+p7s+'||'+p8s+p9s+p10s+$
strtrim(string(k2cat[i].deriv.kp,format='(d6.3)'),2)+'|'+strtrim(string(k2cat[i].id.epic),2)+'|||||'+strtrim(string(glxid[i]),2)+$
'||||||||||||||'+p13s+'|||||'

printf,3,str
printf,2,str3


endfor
close,/all

;stop

skipascii:

; make FITS output

print,'writing FITS output'
temp=k2cat
n=n_elements(temp)
k2cat = replicate({ra:0D, de:0D, bmag:0D, vmag:0D, umag:0D, gmag:0D, rmag:0D, $
	imag:0D, zmag:0D, jmag:0D, hmag:0D, kmag:0D, pmra:0D, pmde:0D, plx:0D, $
	e_bmag:0D, e_vmag:0D, e_umag:0D, e_gmag:0D, e_rmag:0D, $
	e_imag:0D, e_zmag:0D, e_jmag:0D, e_hmag:0D, e_kmag:0D, e_pmra:0D, e_pmde:0D, e_plx:0D, $
	epic:'', hip:'', tycho:'', ucac:'', twomass:'', sdss:'', kepflag:'', stparflag:'', objclass:'', $
    	teff:0D, logg:0D, feh:0D, rad:0D, mass:0D, rho:0D, lum:0D, dis:0D,$
	e_teff:0D, e_logg:0D, e_feh:0D, e_rad:0D, e_mass:0D, e_rho:0D, e_lum:0D, e_dis:0D, $
	ebv:0D, kp:0D },n)

k2cat.ra = temp.obs.ra
k2cat.de = temp.obs.de
k2cat.bmag = temp.obs.bmag
k2cat.vmag = temp.obs.vmag
k2cat.umag = temp.obs.umag
k2cat.rmag = temp.obs.rmag
k2cat.imag = temp.obs.imag
k2cat.zmag = temp.obs.zmag
k2cat.jmag = temp.obs.jmag
k2cat.hmag = temp.obs.hmag
k2cat.kmag = temp.obs.kmag
k2cat.pmra = temp.obs.pmra
k2cat.pmde = temp.obs.pmde
k2cat.plx = temp.obs.plx
k2cat.e_bmag = temp.obs.e_bmag
k2cat.e_vmag = temp.obs.e_vmag
k2cat.e_umag = temp.obs.e_umag
k2cat.e_rmag = temp.obs.e_rmag
k2cat.e_imag = temp.obs.e_imag
k2cat.e_zmag = temp.obs.e_zmag
k2cat.e_jmag = temp.obs.e_jmag
k2cat.e_hmag = temp.obs.e_hmag
k2cat.e_kmag = temp.obs.e_kmag
k2cat.e_pmra = temp.obs.e_pmra
k2cat.e_pmde = temp.obs.e_pmde
k2cat.e_plx = temp.obs.e_plx

k2cat.epic = temp.id.epic
k2cat.hip = temp.id.hip
k2cat.tycho = temp.id.tycho
k2cat.ucac = temp.id.ucac
k2cat.twomass = temp.id.mass
k2cat.sdss = temp.id.sdss
k2cat.kepflag = temp.id.kepflag
k2cat.stparflag = temp.id.stparflag
k2cat.objclass = temp.id.class

k2cat.kp = temp.deriv.kp

k2cat[*].e_teff =  0.
k2cat[*].e_logg =  0.
k2cat[*].e_feh =  0.
k2cat[*].e_rad =  0.
k2cat[*].e_rho =  0.
k2cat[*].e_lum =  0.
k2cat[*].e_dis =  0.
k2cat[*].e_mass =  0.

mwrfits,k2cat,'output/'+filename+'.fits',/create
;save,file=filename+'.dmc.idl',k2cat

end
