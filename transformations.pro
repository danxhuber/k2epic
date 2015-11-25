; various color transformations; see documentation for details
function jhktov,j,h,k
bv = 1.622*(j-h)+0.912*(h-k)+0.044
ri = 0.954*(j-h)+0.593*(h-k)+0.025
v = 1.21*bv+1.295*ri-0.046+j
return,v
end

function bvtog,b,v
gr = 1.124*(b-v) - 0.252
g = v + 0.634*(b-v) - 0.108
return,g
end

function jhktog,j,h,k
gr = 1.951*(j-h)+1.199*(h-k)-0.230
ri = 0.991*(j-h)+0.792*(h-k)-0.210
g = 1.379*gr+1.702*ri+0.518+j
return,g
end

function tosloan,gp,rp,ip
;u = up -           ((up-gp)-1.39)
if (gp ne 0. and rp ne 0.) then g = gp - (-0.06) * ((gp-rp)-(0.53)) else g = 0.
if (rp ne 0. and ip ne 0.) then r = rp - (-0.035)* ((rp-ip)-(0.21)) else r = 0.
if (rp ne 0. and ip ne 0.) then i = ip - (-0.041)* ((rp-ip)-(0.21)) else i = 0.
;z = zp - (0.03)  * ((ip-zp)-(0.09))
return,[g,r,i]
end
