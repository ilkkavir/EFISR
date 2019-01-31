function fighandle = plotEfieldMlat( EfVi , varargin )
%
% plotEfieldMlat( EfVi )
%
% plot electric field components as a pseudocolor plot as function
% of time and magnetic latitutde
%
% INPUT:
%  EfVi  an output list from fitEfieldVi
%  other parameters as name-value pairs:
%    nlim      y-axis limits for the northward electric field [mV/m], default
%              [-100 100]
%    elim      y-axis limits for the eastward electric field [mV/m],
%              default [-100 100]
%    mlim      y-axis limits [degrees], default: the range of
%              magnetic latitudes in EfVi
%    stdlim    standard deviation limit [mV/m], points with std larger
%              than stdlim are not plotted. Default 200
%    chisqrlim chi-square limit, points with chi-square larger than
%              chisqrlim are not plotted, default 10
%    starttime     x-axis start time as a matlab datetime
%    endtime     x-axis end tiem as a matlab datetime
%
% OUTPUT:
%  fighandle matlab figure handles to the figures
%
% IV 2017, 2018
%

% parse the inputs
p = inputParser;

defaultNlim = [-1 1]*100;
checkNlim = @(x) (isnumeric(x) & length(x)==2 );

defaultElim = [-1 1]*100;
checkElim = @(x) (isnumeric(x) & length(x)==2 );

defaultMlim = [NaN NaN];
checkMlim = @(x) (isnumeric(x) & length(x)==2 );

defaultStdlim = 200;
checkStdlim = @(x) (isnumeric(x) & length(x)==1);

defaultChisqrlim = 10;
checkChisqrlim = @(x) (isnumeric(x) & length(x)==1);

defaultstarttime = NaN;
checkStarttime = @(x) (isdatetime(x));

defaultendtime = NaN;
checkEndtime = @(x) (isdatetime(x));


addRequired( p , 'EfVi' , @isstruct ); % data is always required
addParameter( p , 'nlim' , defaultNlim , checkNlim );
addParameter( p , 'elim' , defaultElim , checkElim );
addParameter( p , 'mlim' , defaultMlim , checkMlim );
addParameter( p , 'stdlim' , defaultStdlim , checkStdlim);
addParameter( p , 'chisqrlim' , defaultChisqrlim , checkChisqrlim);
addParameter( p , 'starttime' , defaultstarttime , checkStarttime );
addParameter( p , 'endtime' , defaultendtime , checkEndtime );

parse(p,EfVi,varargin{:});

elim = p.Results.elim;
nlim = p.Results.nlim;
mlim = p.Results.mlim;

if any(isnan(mlim))
    mlim = [min(EfVi.mlat(:)) max(EfVi.mlat(:))];
end

% the electric field components (mV/m)
Enorth = squeeze( EfVi.E( : , : , 1 ) * 1000);
Eeast =  squeeze(EfVi.E( : , : , 2 ) * 1000);

% remote points with large chi-squared
chisqrMask = EfVi.chisqrVi > p.Results.chisqrlim;
Enorth(chisqrMask) = NaN;
Eeast(chisqrMask) = NaN;

% standard deviations from the covariance matrices
Estdnorth = squeeze( sqrt( EfVi.Ecov( : , : , 1 , 1 ) ) * 1000 );
Estdeast  = squeeze( sqrt( EfVi.Ecov( : , : , 2 , 2 ) ) * 1000 );

% remove points with large std in either component
irem = Estdnorth > p.Results.stdlim | Estdeast > p.Results.stdlim;

Enorth(irem) = NaN;
Eeast(irem)  = NaN;
Estdnorth(irem) = NaN;
Estdeast(irem)  = NaN;

% time as datenum
tt = datenum(datetime(mean(EfVi.time,'omitnan'),'convertfrom','posixtime'));

% the magnetic latitudes
mlat = mean(EfVi.mlat','omitnan');

% x axis limits
starttime = datenum(p.Results.starttime);
endtime = datenum(p.Results.endtime);
if isnan(starttime)
    starttime = min(tt);
end
if isnan(endtime)
    endtime = max(tt);
end
if isnan(starttime)|isnan(endtime)
    fighandle = [];
    return
end
if starttime == endtime
    fighandle = [];
    return
end

% add one row and column to compensate for loosing them in
% pcolor...
[d1 d2] = size(Eeast);
Eeast(d1+1,:) = NaN;
Eeast(:,d2+1) = NaN;
Enorth(d1+1,:) = NaN;
Enorth(:,d2+1) = NaN;
tt = [tt 2*tt(end)-tt(end-1)];
mlat = [mlat 2*mlat(end)-mlat(end-1)];

% open a figure
fighandle = figure;
set( fighandle , 'Position' , [10 10 690 360] , 'PaperPositionMode' ...
                 , 'Auto' );


colAxes = get(gcf,'defaultAxesColorOrder');


% upper panels for the actual fields
% East component on the first panel
h1 = subplot(2,1,1);
pcolor(tt,mlat,Eeast);shading flat;colormap jet;
cb1 = colorbar;
ylim(mlim);
xlim([starttime endtime])
ylabel('aacgm\_v2 lat [deg]')
caxis(elim)

% North component on the second panel
h2 = subplot(2,1,2);
pcolor(tt,mlat,Enorth);shading flat;colormap jet;
cb2 = colorbar;
ylim(mlim);
xlim([starttime endtime])
ylabel('aacgm\_v2 lat [deg]')
caxis(nlim)
xlabel('UT')


xlabel(['UTC  ',datestr(datetime(starttime,'convertfrom','datenum'),'yyyy-mm-dd')])
datetick(h2,'x',13,'keeplimits')

% force update before manipulating the layout
drawnow


% figure settings...
set( h1 , 'XTickLabel' , '' ) % remove x-axis ticks from the upper panel
pos1 = get( h1 , 'Position' ); % position of the first panel
pos2 = get( h2 , 'Position' ); % position of the second panel

% remove empty space from between subplots
set( h1 , 'Position' , [ pos1(1:2)-[0,.02] [1,1.15].*pos1(3:4)] );
set( h2 , 'Position' , [ pos2(1:2)+[0,.02] [1,1.15].*pos1(3:4)] );
set( h1 , 'XTick' , get(h2,'XTick'))

set(h1,'fontsize',12)
set(h2,'fontsize',12)
set(cb1,'fontsize',12)
set(cb2,'fontsize',12)

ylabel(cb1,'E_{east} [mVm^{-1}]')
ylabel(cb2,'E_{north} [mVm^{-1}]')

set(h1,'TickDir','both')
set(h2,'TickDir','both')

drawnow