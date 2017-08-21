function EV  = plasmaVelocity2ElectricField( vel )
%
%
% EV = plasmaVelocity2ElectricField( vel )
%
% Electric field from vector velocities
%
% INPUT:
%   vel    an output struct from fitPlasmaVelocities
%
% OUTPUT:
%   EV     a struct with fields
%
%   E      nGate x nTime x 2 array of  electric field estimates
%   Ecov   nGate x nTime x 2 x 2 error covariance matrix of the
%          electric field estimates
%   ...    all contents of the input struct vel
%
%  IV 2017
%

% dimensions
[nGate , nTime] = size(vel.time);

% allocate the electric field and covariance arrays
E = NaN(nGate,nTime,2);
E = NaN(nGate,nTime,2);


for iG = 1:nGate
    for iT = 1:nTime
        Bsqr = sum(vel.Bned(iG,iT,:).^2)*1e-18;
        E(iG,iT,:) = squeeze(vel.vel(iG,iT,[2,1])).*[-1 1]'*sqrt(Bsqr);
        Ecov(iG,iT,:,:) = vel.velcov(iG,iT,[2,1],[2,1])*Bsqr;
    end
end

EV = vel;
EV.E = E;
EV.Ecov = Ecov;

end