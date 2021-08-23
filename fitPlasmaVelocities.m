function velocity = fitPlasmaVelocities( tres , startTime , gateType ...
                                         , gateLims , maxDiff , ...
                                         ViBzero , ViBEzero , ViBNzero ...
                                         , varargin)
%
% function velocity = fitPlasmaVelocities( tres , gateType ,
%                     gateLims ,  maxDiff , ViPar0 , dpath [,
%                     dpath2 , ...])
%
%
% Fit plasma velocity vectors to incoherent scatter
% measurements. The data can be either from a multistatic radar
% system, from a scanning monostatic system, two or more
% independent radars, or from any combination of these.
%
% INPUT:
%   tres       time resolution [s], OR a vector of unixtimes, OR
%              number of slices to integrate as a negative integer
%   startTime  analysis start time as unix time
%   gateType   type of gating, 'h' or 'mlat'
%   gateLims   gate limits, km if gateType='h', altitude limits and degrees if gateType='mlat'
%              The first two elements are the lower and upper bound and the others are gate limits in degrees of magnetic latitude

%   maxDiff    tolerance for common volume selection. degrees in
%              geomagnetic coordinates
%   ViBzero    logical, 0 for normal fit, 1 to force parallel
%              velocity to zero
%   ViBEzero   logical, 0 for normal fit, 1 to force (magnetic) east
%              velocity to zero
%   ViBNzero   logical, 0 for normal fit, 1 to force (magnetic) north
%              velocity to zero
%   dpath      data path. arbitrary number of paths to GUISDAP
%              output directories and/or individual files
%
% OUTPUT:
%
%  velocity  a struct with fields
%
%   vel       nGate x nTime x 3 array of velocity vectors (m/s) in
%             local NED-coordinates
%   velcov    nGate x nTime x 3 x 3 array of error covariance
%             matrices (m^2/s^2)
%   chisqr    chi-squared of the velocity fit
%   mlat      nGate x nTime array of geomagnetic latitudes (deg)
%   mlon      nGate x nTime array of geomagnetic longitudes (deg)
%   glat      nGate x nTime array of geodetic (wgs84) latitudes (deg)
%   glon      nGate x nTime array of geodetic (wgs84) longitudes (deg)
%   height    nGate x nTime array of heights (km)
%   time      nGate x nTime array of times (unix time)
%   mlt       nGate x nTime array of magnetic local times (hours)
%   Bned      nGate x nTime x 3 array of magnetic field vectors. nT
%             in local cartesian north-east-down coordinates.
%   tlims     nGate x nTime x 2 array of integration time limits (unix time)
%
% IV 2016
%

% read the velocity data
if gateType=='h'
    hlims = [min(gateLims) max(gateLims)];
elseif gateType=='mlat'
    hlims = gateLims(1:2);
else
    error(['Unknown gate type ' type ]);
end

[ v , verr , status , chisqr , ts , te , mlt , llhT , llhR , azelR ...
  , r , h , phi , site , ecefS , llhgS , llhmS , kS , B ] = ...
    readVelocitiesGUISDAP(hlims,varargin{:});

% time-slices, indTime contains a time-slice index for each data point
[ indTime , nTime ] = integrationLimitsTime( ts , te , tres , startTime );

% gates, indGate contains a gate index for each data point
[ indGate , nGate ] = integrationLimitsGate( llhgS , llhmS , gateType ...
                                             , gateLims );
% allocate output arrays
vel = NaN( nGate , nTime , 3 );
velcov = NaN( nGate , nTime , 3 , 3 );
mlat = NaN( nGate ,  nTime );
mlon = NaN( nGate , nTime );
glat = NaN( nGate ,  nTime );
glon = NaN( nGate , nTime );
Bned = NaN( nGate , nTime , 3);
height = NaN( nGate , nTime );
time = NaN( nGate , nTime);
tlims = NaN( nGate , nTime , 2 );
mltime = NaN( nGate , nTime );
chisqrVi = NaN( nGate , nTime );

% loop over time slices
for iT = 1:nTime

    % loop over gates
    for iG = 1:nGate

        % pick data from the correct time slice and gate, do not
        % use unsuccessful fits
        indGT =  ( indTime == iT ) & ( indGate == iG ) & ~status;

        % limit the data to magnetic latitudes and longitudes
        % covered by all sites
        indCV = integrationLimitsCommonVolume( llhmS , site , maxDiff ...
                                               , indGT );
        if ~any(indCV)
            continue
        end
        % pick the velocity data, standard deviations and k-vectors
        VGT = v(indCV);
        VerrGT  = verr(indCV);
        kGT = kS(indCV,:);
        % magnetic field directions
        BGT = B(indCV,:);
        % ecef coordinates
        ecefSGT = ecefS(indCV,:);
        % geodetic latitude, longitude, height
        llhgSGT = llhgS(indCV,:);

        % conversion to local geomagnetic coordinates
        kBGT = ecef2localMagnetic( kGT , ecefSGT , llhgSGT , BGT );

        % if parallel velocity is forced to zero
        if ViBzero
            % a zero-meeasurement
            VGT = [ VGT ; 0 ];
            % a small standard deviation (1 mm/s)
            VerrGT = [ VerrGT ; 1e-3 ];
            % magnetic field direction, we are in geomagnetic
            % coordinates now!
            kBGT = [ kBGT ; [0 0 1] ];
        end

        % if eastward velocity is forced to zero
        if ViBEzero
            % a zero-meeasurement
            VGT = [ VGT ; 0 ];
            % a small standard deviation (1 mm/s)
            VerrGT = [ VerrGT ; 1e-3 ];
            % magnetic field direction, we are in geomagnetic
            % coordinates now!
            kBGT = [ kBGT ; [0 1 0] ];
        end
        % if northward velocity is forced to zero
        if ViBNzero
            % a zero-meeasurement
            VGT = [ VGT ; 0 ];
            % a small standard deviation (1 mm/s)
            VerrGT = [ VerrGT ; 1e-3 ];
            % magnetic field direction, we are in geomagnetic
            % coordinates now!
            kBGT = [ kBGT ; [1 0 0] ];
        end

        % velocity in magnetic coordinates
        [ vel(iG,iT,:) , velcov(iG,iT,:,:) ] = plasmaVelocityVector( ...
            VGT , VerrGT , kBGT );

        % project the velocities to k-vector directions
        velprojs = kBGT * squeeze( vel( iG , iT , : ) );

        % chi-squared
        chisqrVi( iG , iT ) = sum( ((velprojs - VGT ).^2) ./ VerrGT.^2) ...
            / length(velprojs);

        % coordinates
        mlat(iG,iT) = mean( llhmS(indCV,1) );
        mlon(iG,iT) = mean( llhmS(indCV,2) );
        glat(iG,iT) = mean( llhgS(indCV,1) );
        glon(iG,iT) = mean( llhgS(indCV,2) );
        height(iG,iT) = mean( llhgS(indCV,3) )/1000;

        % magnetic field
        Bned(iG,iT,1) = mean(B(indCV,1));
        Bned(iG,iT,2) = mean(B(indCV,2));
        Bned(iG,iT,3) = mean(B(indCV,3));

        % times
        time(iG,iT) = (mean(te(indCV)) + mean(ts(indCV)))/2;
        tlims(iG,iT,1) = min(ts(indCV));
        tlims(iG,iT,2) = max(te(indCV));
        
        % be careful with mlt around the mlt midnight
        mlttmp = mlt(indCV);
        maxmlt = max(mlttmp);
        minmlt = min(mlttmp);
        if maxmlt - minmlt > 20
            mlttmp(mlttmp<12) = mlttmp(mlttmp<12) + 24;
        end
        mltime(iG,iT) =  mod(mean(mlttmp),24);
    end

end

velocity = struct();
velocity.vel = vel;
velocity.velcov = velcov;
velocity.mlat = mlat;
velocity.mlon = mlon;
velocity.glat = glat;
velocity.glon = glon;
velocity.height = height;
velocity.time = time;
velocity.mlt = mltime;
velocity.Bned = Bned;
velocity.chisqrVi = chisqrVi;
velocity.tlims = tlims;



end