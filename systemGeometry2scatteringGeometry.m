function [xyz,llhg,llhm,k] = systemGeometry2scatteringGeometry(llhT,llhR,azelR,r,time)
%
% Convert the measurement system geometry into location of the
% measurement volume and scattering wave vector direction.
%
% INPUT:
%  llhT  latitude, longitude, height of the transmitter (deg,deg,m)
%  llhR  latitude, longitude, height of the receiver (deg,deg,m)
%  azelR azimuth, elevation of the receiver (deg)
%  r     distance from the receiver to the target (m)
%  time  time for which the aacgm coordinates are calculated
%        (matlab datetime)
%
% OUTPUT:
%  xyz   scattering volume location in cartesian ECEF coordinates (m)
%  llhg  latitude, longitude, height of the scattering volume in
%        geodetic coordinates (deg,deg,m)
%  llhm  latitude, longitude, height of the scattering volume in
%        altitude adjusted corrected geomagnetic (aacgm)
%        coordinates  (deg,deg,m). Height is radial height (radial
%        distance - RE)
%
%  k     scattering wave vector direction (unit vector) in
%        cartesian ECEF coordinates
%
% IV 2016
%

DRT = pi/180;

% the reference ellipsoid
wgs84 = wgs84Ellipsoid('meters');

% site locations in ECEF
[xT,yT,zT] = geodetic2ecef(wgs84,llhT(1),llhT(2),llhT(3));
[xR,yR,zR] = geodetic2ecef(wgs84,llhR(1),llhR(2),llhR(3));

% receiver beam pointing in local north-east-up
[beamN,beamE,beamU] = sph2cart(DRT*azelR(1),DRT*azelR(2),r);


% target location in ECEF
[xS,yS,zS] = enu2ecef(beamE,beamN,beamU,llhR(1),llhR(2),llhR(3), ...
                         wgs84,'degrees');


% receiver beam in ECEF
xyzRB = [xS,yS,zS] - [xR,yR,zR];


% transmitter beamin ECEF
xyzTB = [xS,yS,zS] - [xT,yT,zT];

% unit vectors in beam directions, away from the antennas
xyzTB0 = xyzTB./sqrt(sum(xyzTB.^2));
xyzRB0 = xyzRB./sqrt(sum(xyzRB.^2));

% the scattering wave vector direction
% minus sign because this vector must be towards the radar(s)
k = -(xyzTB0 + xyzRB0);
k = k./sqrt(sum(k.^2));

% target location in geodetic lat, lon, h
[latS,lonS,hS] = ecef2geodetic(wgs84,xS,yS,zS);
llhg = [latS,lonS,hS];

% target location in geomagnetic lat, lon, h
[mlatS,mlonS,mhS] = geodetic2aacgm(latS,lonS,hS/1000,time);
llhm = [mlatS,mlonS,mhS*1000];

% cartesian ECEF coordinates of the target
xyz = [xS,yS,zS];


end