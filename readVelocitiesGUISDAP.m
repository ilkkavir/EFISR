function [v,verr,status,chisqr,ts,te,mlt,llhT,llhR,azelR,r,h,phi,site,ecefS,llhgS,llhmS,kS,B] = readVelocitiesGUISDAP( varargin )
%
% [v,verr,status,chisqr,ts,te,mlt,llhT,llhR,azelR,r,h,phi,site,ecefS,llhgS,llhmS,kS,B]
% = readVelocitiesGUISDAP( fpath [, fpath2 , fpath3 , ...] )
%
% Read ion velocities and related information from GUISDAP fit
% results.
%
% INPUT:
%  fpath  path(s) to GUISDAP output files or director/y/ies. Any
%         combination of individual files and directories.
%
% OUTPUT:
%  v         velocity estimates (m/s)
%  verr      standard deviations of the velocity estimates (m/s)
%  status    fit status (0 for success)
%  chisqr    chi-squared
%  ts        integration start time (unix time)
%  te        integration end time (unix time)
%  mlt       magnetic local time at integration end time (hours)
%  llhT      latitude, longitude, height of the transmitter (deg,deg,m)
%  llhR      latitude, longitude, height of the receiver (deg,deg,m)
%  azelR     azimuth and elevation angles of the receiver (deg,deg)
%  r         range (m)
%  h         height (m)
%  phi       scattering angle (deg)
%  site      EISCAT site name (char)
%  ecefS     location of the scattering volume(s) in cartesian ecef
%            (earth-centred earth-fixed) coordinates (m)
%  llhgS     geodetic coordinates of the scattering volume(s). (deg,deg,m)
%  llhmS     altitude-adjusted corrected geomagnetic (v2)
%            coordinates of the scattering volume (deg,deg,m). The
%            last coordinate is geocentric height (radial distance
%            - Re)
%  kS        scattering wave vector direction(s) (unit vector(s)) in
%            cartesian ecef coordinates
%  B         magnetic field vector in local North-East-Down
%            coordinates (nT)
%
%
%
% The function returns all data, including failed iterations. Use
% e.g. the output variables "status" and "chisquare" for removing
% failed iterations and unrealistic results.
%
%
%
%
% IV 2016
%

v = [];
verr = [];
status = [];
chisqr = [];
ts = [];
te = [];
mlt = [];
llhT = [];
llhR = [];
azelR = [];
r = [];
h = [];
phi = [];
site = [];
ecefS = [];
llhgS = [];
llhmS = [];
kS = [];
B = [];

fpath = listmatfiles(varargin{:});

% fpath will be a struct if any mat files were found
if isstruct(fpath)

    nf = length(fpath);

    for ff = 1:nf
        readThisFile = false;
        if ~fpath(ff).isdir
            fname = fpath(ff).name;
            if length(fname)>=12
                if (fname(end-3:end)) == '.mat'
                    tstr = str2num(fpath(ff).name(end-11:end-4));
                    if ~isempty(tstr)
                        if 0 <= tstr <= 31622400
                            readThisFile = true;
                        end
                    end
                end
            end
        end
        if readThisFile
            %Ne,NeStd,Ti,TiStd,Tr,TrStd,Coll,CollStd,Vi,Vistd,Comp,CompStd,status,chisqr,ts,te,mlt,llhT,llhR,azelR,r,h,phi,site,ecefS,llhgS,llhmS,kS,B
            disp(fpath(ff).name)
            paramlist = readParametersFromFileGUISDAP(fpath(ff).name);
            % should pre-allocate!!! the re-allocation takes forever...
            v = [ v ; paramlist.Vi ];
            verr = [ verr ; paramlist.ViStd ];
            status = [ status ; paramlist.status ];
            chisqr = [chisqr ; paramlist.chisqr ];
            ts = [ ts ; paramlist.ts ];
            te = [ te ; paramlist.te ];
            mlt = [ mlt ; paramlist.mlt ];
            llhT = [ llhT ; paramlist.llhT ];
            llhR = [ llhR ; paramlist.llhR ];
            azelR = [ azelR ; paramlist.azelR ];
            r = [ r ; paramlist.r ];
            h = [ h ; paramlist.h ];
            phi = [ phi ; paramlist.phi ];
            site = [ site ; paramlist.site ];
            ecefS = [ ecefS ; paramlist.ecefS ];
            llhgS = [ llhgS ; paramlist.llhgS ];
            llhmS = [ llhmS ; paramlist.llhmS ];
            kS = [ kS ; paramlist.kS ];
            B = [ B ; paramlist.B ];
        end
    end
end

end








