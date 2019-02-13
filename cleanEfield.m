function EfViClean = cleanEfield( EfVi , varargin )
%
% Clean failed fits from electric field data based on chi-square,
% E-N component correlation, and standard deviation.
%
% INPUT:
%    EfVi     an output list from fitEfield
%
%    optional  parameters as name-value pairs
%      gates     indices of height/mlat gates to clean, default:
%                clean all gates
%      stdlim    standard deviation limit [mV/m], points with std larger
%                than stdlim are not plotted. Default 200
%      chisqrlim chi-square limit, points with chi-square larger than
%                chisqrlim are not plotted, default 10
%      corrlim   upper limit for E-field east and north component
%                correlation, default 0..95
%  
% OUTPUT:
%    EfViClean  a cleaned version of EfVi
%
% IV 2019
%

% parse the inputs
p = inputParser;

defaultGates = 1:size( EfVi.E , 1 );
checkGates = @(x) ( isnumeric(x) & all(x>0) & all(x<=size(EfVi.E,1)));

defaultStdlim = 200;
checkStdlim = @(x) (isnumeric(x) & length(x)==1);

defaultChisqrlim = 10;
checkChisqrlim = @(x) (isnumeric(x) & length(x)==1);

defaultCorrlim = 0.95;
checkCorrlim = @(x) (isnumeric(x) & length(x)==1);

addRequired( p , 'EfVi' , @isstruct ); % data is always required
addParameter( p , 'gates' , defaultGates , checkGates );
addParameter( p , 'stdlim' , defaultStdlim , checkStdlim);
addParameter( p , 'chisqrlim' , defaultChisqrlim , checkChisqrlim);
addParameter( p , 'corrlim' , defaultCorrlim , checkCorrlim);

parse(p,EfVi,varargin{:});

EfViClean = EfVi;

% clean the gates one-by-one to make this easier to understand
for gnum = p.Results.gates

    % chi-square
    chisqrMask = EfVi.chisqrVi(gnum,:) > p.Results.chisqrlim;

    % east-north correlation
    corrEN = EfVi.Ecov(gnum,:,1,2) ./ ( sqrt(EfVi.Ecov(gnum,:,1,1)) .* sqrt(EfVi.Ecov(gnum,:,2,2)));
    corrMask = abs(corrEN) > p.Results.corrlim;

    % standard deviation
    Estdnorth = sqrt( EfVi.Ecov( gnum , : , 1 , 1 ) ) * 1000;
    Estdeast  = sqrt( EfVi.Ecov( gnum , : , 2 , 2 ) ) * 1000;
    stdMask = Estdnorth > p.Results.stdlim | Estdeast > p.Results.stdlim;

    irem = chisqrMask | corrMask | stdMask;

    % remove also individual good points with failed fits on both sides
    iirem = irem(3:end) & irem(1:end-2);
    irem(find(iirem)+1) = true;
    if (irem(2))
        irem(1) = true;
    end
    if (irem(end-1))
        irem(end)=true;
    end

    EfViClean.E( gnum , irem , : ) = NaN;
    EfViClean.Ecov( gnum , irem , : , : ) = NaN;
    EfViClean.vel( gnum , irem , : ) = NaN;
    EfViClean.velcov( gnum , irem , : , : ) = NaN;

end