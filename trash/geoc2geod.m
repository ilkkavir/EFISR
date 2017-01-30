function [lat_out,lon_out,r_out] = geoc2geod(lat,lon,r,RE)
%
% [lat_out,lon_out,r_out] = geoc2geoe(lat,lon,r,RE)
%
%  Conversion from geocentric to geodetic coordinates. Modified
%  from the c-function geoc2geoe in igrflib.c
%
% INPUT:
%  lat geocentric latitude (deg)
%  lon geocentric longitude (deg)
%  alt geocentric altitude (radial distance - RE) (km)
%  RE  Earh radius (km)
%
% OUTPUT:
%  lat_out geodetic latitude (deg)
%  lon_out geodetic longitude (deg)
%  r_out   geodetic altitude (km)
%
% IV 2016
%

DTOR = pi/180;
a = 6378.1370; %            /* semi-major axis */
f = 1./298.257223563; %    /* flattening */
b = a*(1. -f);       %      /* semi-minor axis */
ee = (2. - f) * f;
e4 = ee*ee;
aa = a*a;

theta = (90. - lat)*DTOR;
phi   = lon * DTOR;

st = sin(theta);
ct = cos(theta);
sp = sin(phi);
cp = cos(phi);

x = (r+RE) * st * cp;
y = (r+RE) * st * sp;
z = (r+RE) * ct;

k0i   = 1. - ee;
pp    = x*x + y*y;
zeta  = k0i*z*z/aa;
rho   = (pp/aa + zeta - e4)/6.;
s     = e4*zeta*pp/(4.*aa);
rho3  = rho*rho*rho;
t     = power(rho3 + s + sqrt(s*(s+2*rho3)), 1./3.);
u     = rho + t + rho*rho/t;
v     = sqrt(u*u + e4*zeta);
w     = ee*(u + v - zeta)/(2.*v);
kappa = 1. + ee*(sqrt(u+v+w*w) + w)/(u + v);

lat_out = atan2(z*kappa,sqrt(pp))/DTOR;
lon_out = lon;
r_out = sqrt(pp + z*z*kappa*kappa)/ee * (1./kappa - k0i);


end
