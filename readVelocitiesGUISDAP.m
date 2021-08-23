function [v,verr,status,chisqr,ts,te,mlt,llhT,llhR,azelR,r,h,phi,site,ecefS,llhgS,llhmS,kS,B] = readVelocitiesGUISDAP( hlims , varargin )
%
% [v,verr,status,chisqr,ts,te,mlt,llhT,llhR,azelR,r,h,phi,site,ecefS,llhgS,llhmS,kS,B]
% = readVelocitiesGUISDAP( fpath [, fpath2 , fpath3 , ...] )
%
% Read ion velocities and related information from GUISDAP fit
% results.
%
% INPUT:
%  hlims  altitude limits ([hmin,hmax], in km) of data needed for the analysis
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

% read all mat files
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
            %            disp(fpath(ff).name)
            %fprintf('.')
            paramlist = readParametersFromFileGUISDAP(fpath(ff).name);
            % could pre-allocate, but this is fast compared to the
            % function call above
            v = [ v ; paramlist.Vi ];
            verr = [ verr ; paramlist.ViStd ];
            status = [ status ; paramlist.status ];
            chisqr = [chisqr ; paramlist.chisqr ];
            ts = [ ts ; paramlist.ts ];
            te = [ te ; paramlist.te ];
            %            mlt = [ mlt ; paramlist.mlt ];
            llhT = [ llhT ; paramlist.llhT ];
            llhR = [ llhR ; paramlist.llhR ];
            azelR = [ azelR ; paramlist.azelR ];
            r = [ r ; paramlist.r ];
            h = [ h ; paramlist.h ];
            phi = [ phi ; paramlist.phi ];
            site = [ site ; paramlist.site ];
            % ecefS = [ ecefS ; paramlist.ecefS ];
            % llhgS = [ llhgS ; paramlist.llhgS ];
            % llhmS = [ llhmS ; paramlist.llhmS ];
            % kS = [ kS ; paramlist.kS ];
            % B = [ B ; paramlist.B ];
        end
    end
end
%fprintf('\n')


%nmatdata = length(r);
% read all hdf5 files
fpath = listhdf5files(varargin{:});

% fpath will be a struct if any hdf5 files were found
if isstruct(fpath)
    
    nf = length(fpath);
    
    for ff = 1:nf
        readThisFile = false;
        if ~fpath(ff).isdir
            fname = fpath(ff).name;
            if (fname(end-4:end)) == '.hdf5'
                readThisFile = true;
            end
        end
        if readThisFile
            %Ne,NeStd,Ti,TiStd,Tr,TrStd,Coll,CollStd,Vi,Vistd,Comp,CompStd,status,chisqr,ts,te,mlt,llhT,llhR,azelR,r,h,phi,site,ecefS,llhgS,llhmS,kS,B
            %            disp(fpath(ff).name)
            
            [data, metadata] = load_EISCAT_hdf5(fpath(ff).name);
            
            nt = length(data.datenum_1(:));
            nv = size(data.v_i_los,1);
            
            for tt = 1:nt
                % in hdf5 files velocity is positive away
                v = [ v ; -data.v_i_los(:,tt)];
                verr = [ verr ; data.v_i_los_err(:,tt) ];
                status = [ status ; data.status(:,tt) ];
                chisqr = [chisqr ; data.residual(:,tt) ];
                
                ts = [ ts ; repmat(posixtime(datetime(data.datenum_1(tt),'convertfrom','datenum')),nv,1) ];
                te = [ te ; repmat(posixtime(datetime(data.datenum_2(tt),'convertfrom','datenum')),nv,1) ];
                
                llhT = [ llhT ; repmat([data.r_XMITloc1 data.r_XMITloc2 data.r_XMITloc3],nv,1) ];
                llhR = [ llhR ; repmat([data.r_RECloc1 data.r_RECloc2 data.r_RECloc3],nv,1) ];
                azelR = [ azelR ; repmat([data.az(tt) data.el(tt)],nv,1) ];
                r = [ r ; data.range(:,tt)*1000 ];
                h = [ h ; data.height(:,tt)*1000 ];
                phi = [ phi ; repmat(data.r_SCangle.*180./pi,nv,1)];
                site = [ site ; repmat(metadata.name_site,nv,1) ];
                
            end                    
            
        end
    end
end

% the data may contain NaN values that cause problems to us
inan = isnan(v) | isnan(verr) | isnan(status) | isnan(chisqr) | isnan(ts) | isnan(te) | any(isnan(llhT'))' | any(isnan(llhR'))' | any(isnan(azelR'))' | isnan(r) | isnan(h) | isnan(phi) | isnan(site);

iheight = (h/1000 >= hlims(1)) & (h/1000 <= hlims(2));

ikeep  = iheight & ~inan;

% v(inan) = [];
% verr(inan) = [];
% status(inan) = [];
% chisqr(inan) = [];
% ts(inan) = [];
% te(inan) = [];
% llhT(inan,:) = [];
% llhR(inan,:) = [];
% azelR(inan,:) = [];
% r(inan) = [];
% h(inan) = [];
% phi(inan) = [];
% site(inan) = [];

v = v(ikeep);
verr = verr(ikeep);
status = status(ikeep);
chisqr = chisqr(ikeep);
ts = ts(ikeep);
te = te(ikeep);
llhT = llhT(ikeep,:);
llhR = llhR(ikeep,:);
azelR = azelR(ikeep,:);
r = r(ikeep);
h = h(ikeep);
phi = phi(ikeep);
site = site(ikeep);

if length(r)>0

    %    posmerged = round( [llhT(:,1)*100 llhT(:,2)*100 llhT(:,3)/1000 llhR(:,1)*100 llhR(:,2)*100 llhR(:,3)/1000 azelR*10 r/1000] );
    posmerged = [llhT(:,1) llhT(:,2) llhT(:,3) llhR(:,1) llhR(:,2) llhR(:,3) azelR r];

    posunique = unique(posmerged,'rows');

    nunique = size(posunique,1);
    
    for iu = 1:nunique
        % [ecefStmp,llhgStmp,llhmStmp,kStmp] = systemGeometry2scatteringGeometry(...
        %     posunique(iu,1:3).*[.01 .01 1000],...
        %     posunique(iu,4:6).*[.01 .01 1000],...
        %     posunique(iu,7:8).*.1,posunique(iu,9).*1000,...
        %     datetime(ts(1),'convertfrom','posixtime'));
        [ecefStmp,llhgStmp,llhmStmp,kStmp] = systemGeometry2scatteringGeometry(...
            posunique(iu,1:3),...
            posunique(iu,4:6),...
            posunique(iu,7:8),posunique(iu,9),...
            datetime(ts(1),'convertfrom','posixtime'));
        Btmp = [0 0 0];
        if llhgStmp(3) <= 600000
            try
                Btmp = igrfmagm(llhgStmp(3),llhgStmp(1),llhgStmp(2),...
                                   decyear(datetime(te(1),'convertfrom','posixtime')));
            catch
                Btmp = igrfmagm(llhgStmp(3),llhgStmp(1),llhgStmp(2),2020);
            end
            
        end

        for ii = 1:length(r)
            if all(posmerged(ii,:)==posunique(iu,:))
                ecefS(ii,:) = ecefStmp;
                llhgS(ii,:) = llhgStmp;
                llhmS(ii,:) = llhmStmp;
                kS(ii,:) = kStmp;
                B(ii,:) = Btmp;
                mlt(ii) = magneticLocalTime(datetime(ts(ii),'convertfrom','posixtime'),llhmS(ii,2));
            end
        end
        
    end
    
    % for ii = (nmatdata+1):length(r)
    %     fprintf('\r %i / %i   %i',ii,nmatdata,length(r))
    %     [ecefS(ii,:),llhgS(ii,:),llhmS(ii,:),kS(ii,:)] = systemGeometry2scatteringGeometry(llhT(ii,:),...
    %                                                       llhR(ii,:),azelR(ii,:),r(ii),...
    %                                                       datetime(ts(ii),'convertfrom','posixtime'));
    %     if llhgS(ii,3) <= 600000
    %         try
    %             B(ii,:) = igrfmagm(llhgS(ii,3),llhgS(ii,1),llhgS(ii,2),...
    %                                decyear(datetime(te(ii),'convertfrom','posixtime')));
    %         catch
    %             B(ii,:) = igrfmagm(llhgS(ii,3),llhgS(ii,1),llhgS(ii,2),2020);
    %         end
            
    %         mlt(ii) = magneticLocalTime(datetime(ts(ii),'convertfrom','posixtime'),llhmS(ii,2));
    %     end
    % end
end
%fprintf('\n')

end

