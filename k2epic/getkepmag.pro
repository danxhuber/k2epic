; returns a Kepler magnitude based on the gri->Kp equations in Brown et al. (2011)
function getkepmag,g,r,i

if (g ne 0. and r ne 0. and i eq 0.) $
	then if (g-r le 0.8) then return,0.1*g+0.9*r else return,0.2*g+0.8*r
	
if (g ne 0. and r eq 0. and i ne 0.) $
	then if (g-i le -0.05) then return,0.55*g+0.45*i else return,0.3*g+0.7*i	
	
if (g eq 0. and r ne 0. and i ne 0.) $
	then if (r-i le 0.673) then return,0.65*r+0.35*i else return,1.2*r-0.2*i	

if (g ne 0. and r ne 0. and i ne 0.) $
	then if (g-r le 0.3) then return,0.25*g+0.75*r else return,0.3*g+0.7*i
	
end
