function [flist] = listmatfiles( varargin )
%
% given an arbitrary number of directories, concatenate the output
% structs of dir('*.mat'). The "names" fields are modified to
% contain full absolute paths to the files
%
% INPUT:
%  directory paths as strings
%
% OUTPUT:
%  flist a matlab struc similar to output of dir, but with absolute
%        paths and outputs of all input directories concatenated.
%
% IV 2016
%


nd  = length(varargin);

flist = [];

if nd > 0
    for dd = 1:nd
        ddpath = varargin{dd};
        % is this an absolute path?
        if ddpath(1) ~= filesep
            % no, make it absolute
            ddpath = fullfile(pwd,ddpath);
        end
        tmplist = [];
        % is this a directory?
        if isdir(ddpath)
            % list mat files
            tmplist = dir(fullfile(ddpath,'*.mat'));
            % convert into full paths
            ntmp = length(tmplist);
            for ff=1:ntmp
                tmplist(ff).name = fullfile(ddpath,tmplist(ff).name);
            end
            % is this a mat-file?
        elseif exist(ddpath,'file')
            if ddpath(end-3:end) == '.mat'
                tmplist = dir(ddpath);
                tmplist.name = ddpath;
            end
        end
        % concatenate all directory listings. 
        flist = [flist;tmplist];
    end
end

% remove possible duplicates
[ii,ii]=unique({flist.name},'stable');
flist=flist(ii);

end