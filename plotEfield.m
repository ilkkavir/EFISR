function fighandle = plotEfield( EfVi , varargin )
%
% plotEfield( EfVi )
%
% Plot electric field components as function of time.
%
% INPUT:
%  EfVi  an output list from fitEfieldVi
%  other parameters as name-value pairs:
%    gates indices of height/mlat gates to plot, default: plot all
%    nlim      y-axis limits for the northward electric field [mV/m], default
%              [-100 100]
%    elim      y-axis limits for the eastward electric field [mV/m],
%              default [-100 100]
%    glim      y-axis limits for the geographic coordinates plot
%              [degrees], default [0 180]
%    mlim      y-axis limits for the magnetic coordinatess plot
%              [degrees], default [0 180]
%    stdlim    standard deviation limit [mV/m], points with std larger
%              than stdlim are not plotted. Default 200
%    chisqrlim chi-square limit, points with chi-square larger than
%              chisqrlim are not plotted, default 10
%    corrlim   upper limit for E-field east and north component
%              correlation, default 0..95
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

defaultGates = 1:size( EfVi.E , 1 );
checkGates = @(x) ( isnumeric(x) & all(x>0) & all(x<=size(EfVi.E,1)));

defaultNlim = [-1 1]*100;
checkNlim = @(x) (isnumeric(x) & length(x)==2 );

defaultElim = [-1 1]*100;
checkElim = @(x) (isnumeric(x) & length(x)==2 );

defaultGlim = [0 180];
checkGlim = @(x) (isnumeric(x) & length(x)==2 );

defaultMlim = [0 180];
checkMlim = @(x) (isnumeric(x) & length(x)==2 );

defaultStdlim = 200;
checkStdlim = @(x) (isnumeric(x) & length(x)==1);

defaultChisqrlim = 10;
checkChisqrlim = @(x) (isnumeric(x) & length(x)==1);

defaultCorrlim = 0.95;
checkCorrlim = @(x) (isnumeric(x) & length(x)==1);

defaultstarttime = NaN;
checkStarttime = @(x) (isdatetime(x));

defaultendtime = NaN;
checkEndtime = @(x) (isdatetime(x));


addRequired( p , 'EfVi' , @isstruct ); % data is always required
addParameter( p , 'gates' , defaultGates , checkGates );
addParameter( p , 'nlim' , defaultNlim , checkNlim );
addParameter( p , 'elim' , defaultElim , checkElim );
addParameter( p , 'glim' , defaultGlim , checkGlim );
addParameter( p , 'mlim' , defaultMlim , checkMlim );
addParameter( p , 'stdlim' , defaultStdlim , checkStdlim);
addParameter( p , 'chisqrlim' , defaultChisqrlim , checkChisqrlim);
addParameter( p , 'corrlim' , defaultCorrlim , checkCorrlim);
addParameter( p , 'starttime' , defaultstarttime , checkStarttime );
addParameter( p , 'endtime' , defaultendtime , checkEndtime );

parse(p,EfVi,varargin{:});


% use this function recursively to plot all the gates
gatenumfound = 0;
if length(p.Results.gates) > 1
    fighandle = [];
    for gatenum = p.Results.gates
        varlist = varargin;
        for iarg = 1:nargin-1
            if all(isstr(varlist{iarg}))
                if strcmp(varlist{iarg},'gates')
                    varlist{iarg+1} = gatenum;
                    gatenumfound = 1;
                    break;
                end
            end
        end
        if ~gatenumfound
            varlist{nargin} = 'gates';
            varlist{nargin+1} = gatenum;
        end
        fighandle = [ fighandle ; plotEfield(EfVi,varlist{:}) ];
    end
    return;
end

% plot one gate, we have a unique p.results.gates at this point
gnum = p.Results.gates;


EfViClean = cleanEfield( EfVi , 'stdlim' , p.Results.stdlim , ...
                         'chisqrlim' , p.Results.chisqrlim , 'corrlim' ...
                         , p.Results.corrlim );
Enorth = EfViClean.E( gnum , : , 1 ) * 1000;
Eeast =  EfViClean.E( gnum , : , 2 ) * 1000;
Estdnorth = sqrt( EfViClean.Ecov( gnum , : , 1 , 1 ) ) * 1000;
Estdeast  = sqrt( EfViClean.Ecov( gnum , : , 2 , 2 ) ) * 1000;

% time as datetime
tt = datetime(EfViClean.time(gnum,:),'convertfrom','posixtime');

% x axis limits
starttime = p.Results.starttime;
endtime = p.Results.endtime;
if ~isdatetime(starttime)
    starttime = min(tt);
end
if ~isdatetime(endtime)
    endtime = max(tt);
end
if isnat(starttime)|isnat(endtime)
    fighandle = [];
    return
end
if starttime == endtime
    fighandle = [];
    return
end

% for error bars (matlab errorbar-function did not prduce
% satisfactory results...)
tterr = [tt;tt;[tt(2:end) tt(end)]];
errnorth = [ Enorth - Estdnorth; Enorth + Estdnorth ; Enorth*NaN];
erreast = [ Eeast - Estdeast; Eeast + Estdeast ; Eeast*NaN];

% open a figure
fighandle = figure;
set( fighandle , 'Position' , [10 10 690 520] , 'PaperPositionMode' ...
                 , 'Auto' );


colAxes = get(gcf,'defaultAxesColorOrder');


% upper panels for the actual fields

% East component on the first panel
h1 = subplot(3,1,1);
%yyaxis('left')
plot([min(tt) max(tt)]+[-1 1]*1e4,[0 0],'-k','LineWidth',1)
hold on
plot(tterr,erreast,'-','Color','red','LineWidth',1.2);
plot(tt,Eeast,'LineWidth',1.5,'Color','black');
ylim(p.Results.elim);
xlim([starttime endtime])
ylabel('E_{east} [mV/m]')
grid on

% North component on the second panel
h2 = subplot(3,1,2);
%yyaxis('right')
plot([min(tt) max(tt)]+[-1 1]*1e4,[0 0],'-k','LineWidth',1)
hold on
plot(tterr,errnorth,'-','Color','red','LineWidth',1.2);
plot(tt,Enorth,'LineWidth',1.5,'Color','black');
ylim(p.Results.nlim);
xlim([starttime endtime])
ylabel('E_{north} [mV/m]')
grid on

% coordinates on the lower panel
h3 = subplot(3,1,3);

% geodetic coordinates on the left axis
yyaxis('left');
plot(tt,EfViClean.glat(gnum,:),'LineWidth',1.5);
hold on
plot(tt,EfViClean.glon(gnum,:),'LineWidth',1.5);
ylim( p.Results.glim )
xlim([starttime endtime])
ylabel('Degrees')

% aacgm_v2 coordinates on the right axis
yyaxis('right')
plot(tt,EfViClean.mlat(gnum,:),'LineWidth',1.5);
hold on
plot(tt,EfViClean.mlon(gnum,:),'LineWidth',1.5);
ylim( p.Results.mlim )
xlim([starttime endtime])
ylabel('Degrees')


legend('WGS84 lat','WGS84 lon','aacgm\_v2 lat','aacgm\_v2 lon' )
xlabel(['UTC  ',datestr(datetime(starttime,'convertfrom','posixtime'),'yyyy-mm-dd')])
datetick(h3,'x',13,'keeplimits')
grid on

% force update before manipulating the layout
drawnow


% figure settings...
set( h1 , 'XTickLabel' , '' ) % remove x-axis ticks from the upper panels
set( h2 , 'XTickLabel' , '' ) % remove x-axis ticks from the upper panels
pos1 = get( h1 , 'Position' ); % position of the first panel
pos2 = get( h2 , 'Position' ); % position of the second panel
pos3 = get( h3 , 'Position' ); % position of the third panel

% remove empty space from between subplots
set( h1 , 'Position' , [ pos1(1:2)-[0,.07] [1,1.15].*pos1(3:4)] );
set( h2 , 'Position' , [ pos2(1:2)-[0,.07]./2 [1,1.15].*pos1(3:4)] );
set( h3 , 'Position' , [ pos3(1:2) [1,1.15].*pos1(3:4)] );
set( h1 , 'XTick' , get(h3,'XTick'))
set( h2 , 'XTick' , get(h3,'XTick'))

% MLT ticks on the top panel
% the interpolation fails when the x axis is longer than the data vector
%dnum = datetime(EfViClean.time,'convertfrom','posixtime');
%mlt = EfViClean.mlt;
%inan =  isnan(mlt);
%dnum = dnum(~inan);
%mlt = mlt(~inan);
%for kk=2:length(mlt)
%    while (mlt(kk)<mlt(kk-1))
%        mlt(kk)=mlt(kk)+24;
%    end
%end
%mltticks = interp1(dnum,mlt,get(h3,'XTick'),'linear','extrap');
%mltticks = mod(mltticks,24);
utticks = get(h3,'XTick');
mltticks = NaN(length(utticks),1);
for kk=1:length(utticks)
    [mlat,mlon,mh] = geodetic2aacgm(EfViClean.glat(gnum,kk),EfViClean.glon(gnum,kk),EfViClean.height(gnum,kk),utticks(kk));
    mltticks(kk) = magneticLocalTime(utticks(kk),mlon);
end
set(h1,'XAxisLocation','top','XTickLabel', cellstr(num2str(mltticks,'%5.2f')))
xlabel(h1,'MLT (aacgm\_v2)')

set(h1,'fontsize',12)
set(h2,'fontsize',12)
set(h3,'fontsize',12)
linkaxes([h1 h2 h3],'x')

drawnow