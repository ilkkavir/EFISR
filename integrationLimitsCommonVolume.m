function [ indCV ] = integrationLimitsCommonVolume( llhm , site , ...
                                                  maxDiff , indGT )
%
% [ indCV ] = integrationLimitsCommonVolume( llhm , site , maxDiff )
%
% pick data points which are close to common volume of all sites
%
%  INPUT:
%    llhm    magnetic latitude, longitude, and height
%    site    site index vector
%    maxDiff maximum difference in magnetic latitude and longitude,
%            degrees
%    indGT    indices of data points to check
%
%
%  OUTPUT:
%   indCV    vector of data indices
%
% IV 2017
%

% pick the relevant part of data
llhmGT = llhm(indGT,:);
siteGT = site(indGT);

% sites
sites = unique(siteGT);

% number of sites
nsites = length(sites);

% the velocity estimation does not make sense with one site...
if nsites < 2
    indCV = [];
    return
end

% magnetic latitude coverage for each site
mlatMin = -Inf(nsites,1);
mlatMax = Inf(nsites,1);
mlonMin = -Inf(nsites,1);
mlonMax = Inf(nsites,1);

for iS = 1:nsites
    iSite = siteGT == sites(iS);
    mlatMin(iS) = min(llhmGT(iSite,1));
    mlatMax(iS) = max(llhmGT(iSite,1));
    mlonMin(iS) = min(llhmGT(iSite,2));
    mlonMax(iS) = max(llhmGT(iSite,2));
end

% our final limits in latidue and longitude
minlat = max(mlatMin) - maxDiff;
maxlat = min(mlatMax) + maxDiff;
minlon = max(mlonMin) - maxDiff;
maxlon = min(mlonMax) + maxDiff;

% the final index vector
indCV = indGT & (llhm(:,1)>=minlat) & (llhm(:,1)<=maxlat) & ...
        (llhm(:,2)>=minlon) & (llhm(:,2)<=maxlon);

end