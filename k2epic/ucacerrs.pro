; this function returns typical errors of UCAC photometry based on power law fits to a 
; sample with reliable uncertainties
function ucacerrs,inmag,bmag=bmag,vmag=vmag,gmag=gmag,rmag=rmag,imag=imag

inmag = float(inmag)

if keyword_set(bmag) then begin
    res = 4.+((inmag-10.)*0.21)^4.5
    u = where(inmag lt 10.)
    if (u[0] ne -1) then res[u] = 4.
endif
    
if keyword_set(vmag) then begin
    res = 3.+((inmag-10.)*0.22)^4.5
    u = where(inmag lt 10.)
    if (u[0] ne -1) then res[u] = 3.
endif

if keyword_set(gmag) then begin
    res = 3.+((inmag-10.)*0.2)^4.2
    u = where(inmag lt 10.)
    if (u[0] ne -1) then res[u] = 3.
endif

if keyword_set(rmag) then begin
    res = 4.+((inmag-10.)*0.22)^4.5
    u = where(inmag lt 10.)
    if (u[0] ne -1) then res[u] = 4.
endif

if keyword_set(imag) then begin
    res = 6.+((inmag-10.)*0.26)^4.5
    u = where(inmag lt 10.)
    if (u[0] ne -1) then res[u] = 6.
endif

return,res
end
