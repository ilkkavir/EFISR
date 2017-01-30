function [ indGate , nGate ] = integrationLimitsGate( llhg , llhm , ...
                                                  type , limits )
%
% [indGate,nGate] = integrationLimitsGate(llhg,llhm,type,limits)
%
% gate indices for data points
%
% INPUT:
%  llhg    geodetic latitude, longitude, height
%  llhm    magnetic latitude, longitude, height
%  type    gating type, 'h' or 'mlat'
%  limits  gate limits, height in km if type=='h', magnetic
%          latitudes in degress if type=='mlat'
%
%  OUTPUT:
%   indGate  a vector of gate indices for each data point
%   nGate    number of gates
%
% IV 2017
%

% number of gates
nGate = length(limits) - 1;

% number of data points
nd = size(llhg,1);

% allocate the index vector, use -1 for points outside the
% integration limits
indGate = -ones(nd,1);

% pick the correct data vector
if type=='h'
    dGate = llhg(:,3)./1000; % the limits are given in km
elseif type=='mlat'
    dGate = llhm(:,1);
else
    error(['Unknown gate type ' type ]);
end

% the indices
for iG=1:nGate
    indGate( (dGate >= limits(iG) ) & (dGate < limits(iG+1)) ) = iG;
end

end
