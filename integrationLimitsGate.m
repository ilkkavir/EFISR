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
%  limits  gate limits, height in km if type=='h', height limits and magnetic
%          latitudes in degress if type=='mlat'. The first two elements are the lower and upper bound
%          and the others are gate limits in degrees of magnetic latitude
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
    hlims = [-Inf Inf];
elseif type=='mlat'
    dGate = llhm(:,1);
    nGate = nGate - 2; % the first two elements are the altitude
                       % limits
    hlims = limits(1:2);
    limits = limits(3:end);
else
    error(['Unknown gate type ' type ]);
end
hGate = llhg(:,3)./1000; % heights 

% the indices
for iG=1:nGate
    indGate( (dGate >= limits(iG) ) & (dGate < limits(iG+1)) & hGate ...
             > hlims(1) & hGate < hlims(2) ) = iG;
end

end
