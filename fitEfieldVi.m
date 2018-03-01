function EfVi = fitEfieldVi( tres , startTime , gateType , gateLims ...
                             , maxDiff , ViBzero , ViBEzero , ViBNzero ...
                             , varargin)
%
% EfVi = fitEfieldVi( tres , starttime , gateType , gaetLims , maxDiff , ViBzero , varargin)
%
% INPUT:
%
%   tres       time resolution [s], OR a vector of unixtimes, OR
%              number of slices to integrate as a negative integer
%   startTime  analysis start time as unix time, use negative value
%              to start from first data point
%   gateType   type of gating, 'h' or 'mlat'
%   gateLims   gate limits, km if gateType='h', degrees if gateType='mlat'
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
%
% OUTPUT:
%
%  velocity  a struct with fields
%
%   vel       nGate x nTime x 3 array of velocity vectors (m/s)
%   velcov    nGate x nTime x 3 x 3 array of error covariance
%             matrices (m^2/s^2)
%   mlat      nGate x nTime array of geomagnetic latitudes (deg)
%   mlon      nGate x nTime array of geomagnetic longitudes (deg)
%   glat      nGate x nTime array of geodetic (wgs84) latitudes (deg)
%   glon      nGate x nTime array of geodetic (wgs84) longitudes (deg)
%   height    nGate x nTime array of heights (km)
%   time      nGate x nTime array of times (unix time)
%   tlims     nGate x nTime x 2 array of integration time limits (unix time)
%   mlt       nGate x nTime array of magnetic local times (hours)
%   Bned      nGate x nTime x 3 array of magnetic field vectors. nT
%             in local cartesian north-east-down coordinates.
%   E      nGate x nTime x 2 array of  electric field estimates
%   Ecov   nGate x nTime x 2 x 2 error covariance matrix of the
%          electric field estimates
%
% IV 2017
%

% velocities
vel = fitPlasmaVelocities( tres , startTime , gateType , gateLims , ...
                           maxDiff , ViBzero , ViBEzero , ViBNzero , varargin{:} );

% electric field
EfVi = plasmaVelocity2ElectricField( vel );

end