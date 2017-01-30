function [lat_out,lon_out,r_out] = geod2geoc(lat,lon,alt,RE)
%
% [lat_out,lon_out,r_out] = geod2geoc(lat,lon,alt,RE)
%
%  Conversion from geodetic to geocentric coordinates. Modified
%  from the c-function geod2geoc in igrflib.c
%
% INPUT:
%  lat geodetic latitude (deg)
%  lon geodetic longitude (deg)
%  alt geodetic altitude (km)
%  RE  Earh radius (km)
%
% OUTPUT:
%  lat_out geocentric latitude (deg)
%  lon_out geocentric longitude (deg)
%  r_out   geocentric altitude (radial distance - RE) (km)
%
% IV 2016
%

DTOR = pi/180;

a = 6378.1370;%/* semi-major axis */
f = 1./298.257223563; %/* flattening */
b = a*(1. -f);%/* semi-minor axis */
a2 = a*a;
b2 = b*b;
theta = (90. -lat)*DTOR;%/* colatitude in radians   */
st = sin(theta);
ct = cos(theta);
one = a2*st*st;
two = b2*ct*ct;
three = one + two;
rho = sqrt(three);%/* [km] */
r = sqrt(alt*(alt+2*rho) + (a2*one + b2*two)/three);%    /* [km] ...*/
cd = (alt+rho)/r;
sd = (a2-b2)/rho *ct*st/r;

r_out = r - RE;
lat_out = 90 - acos(ct*cd - st*sd)/DTOR;
lon_out = lon;

end