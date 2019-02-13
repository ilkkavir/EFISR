function [ indTime , nTime ] = integrationLimitsTime( ts, te , tres , ...
                                                  startTime)
%
%
% calculate indices for integration periods. This is very simple
% now, I may implement other optionsn in future
%
%
% INPUT:
%  ts         integration start times
%  te         integration end times
%  tres       time resolution [s], OR a vector of unixtimes, OR
%             number of slices to integrate as a negative integer
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
    %    startTime = round(te(1));
    startTime = round(min(te));
end

if tres(1)>0
    
    % create the time-slice limits
    if length(tres==1)
        tlims = startTime:tres:(max(te)+tres);
    else
        tlims = tres;
    end
    
    % number of slices
    nTime = length(tlims) - 1;
    
    % number of data points
    nd = length(te);
    
    % allocate the index vector, use -1 for points outside the
    % integration limits
    indTime = -ones(nd,1);
    
    % the indices
    for iT=1:nTime
        indTime( te >= tlims(iT) & te < tlims(iT+1) ) = iT;
    end
    
elseif tres(1)<0

    % arrange the integration periods with start time
    tse  = [ts(:)';te(:)']';
    [~,its] = sort(tse(:,1));
    tse = tse(its,:);

    % number of data points
    nd = length(te);

    % initialize the index vector
    indTime = -ones(nd,1);

    % find the period indices for the sorted data
    iT = 1;
    ii = 1;
    while iT <= nd
        iMask = tse(:,1) >= tse(iT,1) & tse(:,1) < tse(iT,2);
        indTime(iMask) = ii;
        ii = ii+1;
        iT = iT + sum(iMask);
    end

    % put the indices in the correct order
    indTime(its) = indTime;

    % the indices are ordered with increasing time, the time
    % integration reduces to rounding
    indTime = ceil(indTime/abs(tres(1)));

    % number of time slices
    nTime = length(unique(indTime));

else
    error(['Invalide tres ',num2str(tres)])
end
end