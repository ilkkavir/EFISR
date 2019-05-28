function fname = writeEfieldASCII( ef , varargin )
% saveEfield( ef , ...)
% 
% Save the electric field and ion velocity data in an ASCII file with an
% automatically generated file name. 
%
% INPUT:
%  ef   an output list from fitEfVi
%  ...  optional parameters as name value pairs
%     path path to an existing directory, where teh results are
%          saved, default '.'
%
% OUTPUT:
%   none
%
% IV 2019
%


% parse the inputs
p = inputParser;
defaultSavepath = '.';
checkSavepath = @(x) exist(x,'dir');
addRequired(p,'ef',@isstruct);
addParameter(p,'path',defaultSavepath,checkSavepath);
parse(p,ef,varargin{:})

% the file name
fname = fullfile(p.Results.path,[datestr(datetime(round(min(ef.time(:))),'ConvertFrom','posixtime'),'yyyymmddTHHMMss'),'-',datestr(datetime(round(max(ef.time(:))),'ConvertFrom','posixtime'),'yyyymmddTHHMMss'),'_Efield_Vi.dat']);

% header
headstr1 = '      Integration period            Location              Ion velocity     Electric field';
headstr2 = ' YYMM DDHH MMSS YYMM DDHH MMSS   LAT   LON Height   North    East    Down   North    East';

% the data matrix
nd = prod(size(ef.time));
dmat = NaN(nd*2,20);
dline = 1;
for tt=1:size(ef.time,2)
    for gg=1:size(ef.time,1)
        tts = datetime(ef.tlims(gg,tt,1),'convertfrom','posixtime');
        tte = datetime(ef.tlims(gg,tt,2),'convertfrom','posixtime');
        dmat(dline,1) = mod(tts.Year,100);
        dmat(dline,2) = tts.Month;
        dmat(dline,3) = tts.Day;
        dmat(dline,4) = tts.Hour;
        dmat(dline,5) = tts.Minute;
        dmat(dline,6) = round(tts.Second);
        dmat(dline,7) = mod(tte.Year,100);
        dmat(dline,8) = tte.Month;
        dmat(dline,9) = tte.Day;
        dmat(dline,10) = tte.Hour;
        dmat(dline,11) = tte.Minute;
        dmat(dline,12) = round(tte.Second);
        dmat(dline,13) = ef.glat(gg,tt);
        dmat(dline,14) = ef.glon(gg,tt);
        dmat(dline,15) = ef.height(gg,tt);
        dmat(dline,16) = ef.vel(gg,tt,1);
        dmat(dline,17) = ef.vel(gg,tt,2);
        dmat(dline,18) = ef.vel(gg,tt,3);
        dmat(dline,19) = ef.E(gg,tt,1)*1000;
        dmat(dline,20) = ef.E(gg,tt,2)*1000;

        dmat(dline+1,16) = sqrt(ef.velcov(gg,tt,1,1));
        dmat(dline+1,17) = sqrt(ef.velcov(gg,tt,2,2));
        dmat(dline+1,18) = sqrt(ef.velcov(gg,tt,3,3));
        dmat(dline+1,19) = sqrt(ef.Ecov(gg,tt,1,1))*1000;
        dmat(dline+1,20) = sqrt(ef.Ecov(gg,tt,2,2))*1000;

        dline = dline + 2;

    end
end

% replace NaN with -9999
dmat(isnan(dmat)) = -9999;

% write to file
fid = fopen(fname,'w');
fprintf(fid,'%s\n',headstr1);
fprintf(fid,'%s\n',headstr2);
% loop to allow printing the blank space on every second line
for dline=1:2:size(dmat,1)
    % do not write lines with missing timestamps (completely
    % missing data)
    if dmat(dline,1) > -9999
        fprintf(fid,[' %2u%2u %2u%2u %2u%2u %2u%2u %2u%2u %2u%2u%7.2f' ...
                     '%6.2f%6.1f%8.1f%8.1f%8.1f%8.1f%8.1f\n'],dmat(dline,1:20));
        fprintf(fid,['                                                 ' ...
                     '%8.1f%8.1f%8.1f%8.1f%8.1f\n'],dmat(dline+1,16:20));
    end
end

fclose(fid);

end