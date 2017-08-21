function [ indTime , nTime ] = integrationLimitsTime( t , tres , ...
                                                  startTime)
%
%
% calculate indices for integration periods. This is very simple
% now, I may implement other optionsn in future
%
%
% INPUT:
%  t          time points
%  tres       time resolution
%  startTime  beginning  of first time slice
%
% OUTPUT:
%  indTime    time slice indices for each data point
%  nTime      number of time slices
%
% IV 2017
%

% start from beginning of data, if starttime is not specified
if startTime < 0 
    startTime = t(1);
end

% create the time-slice limits
tlims = startTime:tres:max(t);

% number of slices
nTime = length(tlims) - 1;

% number of data points
nd = length(t);

% allocate the index vector, use -1 for points outside the
% integration limits
indTime = -ones(nd,1);

% the indices
for iT=1:nTime
    indTime( t >= tlims(iT) & t < tlims(iT+1) ) = iT;
end

end