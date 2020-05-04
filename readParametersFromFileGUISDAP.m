function param = readParametersFromFileGUISDAP( fpath )
%function [Ne,NeStd,Ti,TiStd,Tr,TrStd,Coll,CollStd,Vi,ViStd,Comp,CompStd,status,chisqr,ts,te,mlt,llhT,llhR,azelR,r,h,phi,site,ecefS,llhgS,llhmS,kS,B] = readParametersFromFileGUISDAP( fpath )
%
% [Ne,NeStd,Ti,TiStd,Tr,TrStd,Coll,CollStd,Vi,Vistd,Comp,CompStd,status,
% chisqr,ts,te,mlt,llhT,llhR,azelR,r,h,phi,site,ecefS,llhgS,llhmS,kS,B]
% = readParametersFromFileGUISDAP( fpath )
%
% Read plasma parameters and related information from GUISDAP fit
% result files. This function reads exactly one file
%
% INPUT:
%  fpath  path to the GUISDAP output file
%
% OUTPUT:
%
%  param a struct with fields:
%
%   Ne        Electron density (m^-3)
%   NeStd     Ne standard deviation (m^-3)
%   Ti        Ion temperature (K)
%   TiStd     Ti standard deviation (K)
%   Tr        Temperature ratio (electron/ion) (dimensionless)
%   TrStd     Tr standard deviation
%   Coll      Ion-neutral collision frequency (s^-1)
%   CollStd   Coll standard deviation (s^-1)
%   Vi        Line-of-sight ion velocity (ms^-1)
%   ViStd     Vi standard deviation (ms^-1)
%   status    fit status (0 for success)
%   chisqr    chi-squared
%   ts        integration start time (unix time)
%   te        integration end time (unix time)
%   mlt       magnetic local time at integration end time (hours)
%   llhT      latitude, longitude, height of the transmitter (deg,deg,m)
%   llhR      latitude, longitude, height of the receiver (deg,deg,m)
%   azelR     azimuth and elevation angles of the receiver (deg,deg)
%   r         range (m)
%   h         height (m)
%   phi       scattering angle (deg)
%   site      EISCAT site name (char)
%   ecefS     location of the scattering volume(s) in cartesian ecef
%             (earth-centred earth-fixed) coordinates (m)
%   llhgS     geodetic coordinates of the scattering volume(s). (deg,deg,m)
%   llhmS     altitude-adjusted corrected geomagnetic (v2)
%             coordinates of the scattering volume (deg,deg,m). The
%             last coordinate is geocentric height (radial distance
%             - Re)
%   kS        scattering wave vector direction(s) (unit vector(s)) in
%             cartesian ecef coordinates
%   B         magnetic field vector in local North-East-Down
%            coordinates (nT)
%
%
% IV 2016
%


% read the data
load(fpath);

param = struct();

% parameters and their standard deviations
param.Ne = r_param(:,1);
param.Ti = r_param(:,2);
param.Tr = r_param(:,3);
param.Coll = r_param(:,4);
param.Vi = r_param(:,5);
param.Comp = r_param(:,6);

param.NeStd = r_error(:,1);
param.TiStd = r_error(:,2);
param.TrStd = r_error(:,3);
param.CollStd = r_error(:,4);
param.ViStd = r_error(:,5);
param.CompStd = r_error(:,6);


% fit status
param.status = r_status;

% chi-squared
param.chisqr = r_res(:,1);

% number of velocity points
nv = length(param.Vi);

% integration start and end time as unix time
param.ts = repmat( posixtime(datetime(r_time(1,1),r_time(1,2),r_time(1,3),r_time(1,4), ...
             r_time(1,5),r_time(1,6))) , nv , 1 );
param.te = repmat( posixtime(datetime(r_time(2,1),r_time(2,2),r_time(2,3),r_time(2,4), ...
             r_time(2,5),r_time(2,6))) , nv , 1 );

% transmitter location (lat,lon,height)
param.llhT = repmat( r_XMITloc , nv , 1 );

% receiver location (lat,lon,height)
param.llhR = repmat( r_RECloc , nv , 1 );

% azimuth and elevation of the receiver antenna
param.azelR = repmat( [ r_az , r_el ] , nv , 1 );

% range (distance from target to receiver? check this!)
param.r = r_range.*1000;

% height as calculated by GUISDAP
param.h = r_h.*1000;

% scattering angle calculated by GUISDAP
param.phi = r_SCangle.*180./pi;

% site name
param.site = repmat( name_site , nv , 1 );

% scattering volume locations etc.
param.ecefS =  zeros(nv,3);
param.llhgS =  zeros(nv,3);
param.llhmS =  zeros(nv,3);
param.kS    =  zeros(nv,3);
param.B     =  zeros(nv,3);
param.mlt   =  zeros(nv,1);

% should calculate only unique values like in Daedalus/calculateCoordinatesMadrigal.m!
for hh = 1:nv

    [param.ecefS(hh,:),param.llhgS(hh,:),param.llhmS(hh,:),param.kS(hh,:)] = ...
        systemGeometry2scatteringGeometry(param.llhT(hh,:),param.llhR(hh,:), ...
                                          param.azelR(hh,:),param.r(hh), ...
                                          datetime(param.ts(hh,:), ...
                                                   'ConvertFrom', ...
                                                   'posixtime'));
    % magnetic field direction from igrf
    % maximum altitude of igrf is 600 km...
    if param.llhgS(hh,3) <= 600000
        try
            param.B(hh,:) = igrfmagm( param.llhgS(hh,3) , param.llhgS(hh,1) , param.llhgS(hh,2), ...
                                      decyear(r_time(2,:)));
        catch
            param.B(hh,:) = igrfmagm( param.llhgS(hh,3) , param.llhgS(hh,1) , param.llhgS(hh,2), ...
                                      2020);
    end


    % magnetic local time
    param.mlt(hh) = magneticLocalTime(datetime(param.ts(hh,:),'ConvertFrom','posixtime'),param.llhmS(hh,2));

end


end