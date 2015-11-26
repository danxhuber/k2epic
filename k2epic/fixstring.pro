; fix screwed up strings in GDL
function fixstring,in
temp=byte(in)
u=where(temp ne 0)
if (u[0] ne -1) then return,string(temp[u]) else return,string(temp)
end