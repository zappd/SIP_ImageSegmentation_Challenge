function [varargout] = enviread(varargin)
% ENVIREAD Easily import ENVI raster files (BSQ,BIL,BIP) from header info.
%   Z = ENVIREAD(FILENAME); Reads an ENVI binary file (BSQ,BIL,BIP) into
%       an array using the information from the corresponding header
%       file FILENAME.hdr. The output array will be of dimensions
%       (m,n,b) where b is the number of bands, m is the number of samples
%        and n is the number of lines. 
%   Z = ENVIREAD(FILENAME,HEADERFILE); Uses the header information in
%       headerfile.
%   [Z,X,Y] = ENVIREAD(....); Returns the map coordinate vectors for geo-
%       registered data.
%   [Z,X,Y,info] = ENVIREAD(....); Returns the header information as a
%       structure.
%
%   NOTES:  -Requires READ_ENVIHDR to read header data.
%           -Geo-registration does not currently support rotated images.
%
% Version 1: 11-Jul-2007 15:11:13
% 17-Sep-2007 : modified line 37 (if isfield(info,'map_info')) to correct
%               error, 
%               modified fopen in line 94-96


%% check for header reader
if exist('read_envihdr.m','file') == 0
    error('This function requires READ_ENVIHDR.m')
end
%% READ HEADER INFO
% Read Filename
file = varargin{1};
hdrfile = strcat(file,'.hdr');
if nargin == 2
    hdrfile = varargin{2};
end
% Get image size and map data from header
info = read_envihdr(hdrfile);
%% Make geo-location vectors
if isfield(info,'map_info')
    if isfield(info.map_info,'mapx') && isfield(info.map_info,'mapy')
        xi = info.map_info.image_coords(1);
        yi = info.map_info.image_coords(2);
        xm = info.map_info.mapx;
        ym = info.map_info.mapy;
        %adjust points to corner (1.5,1.5)
        if yi > 1.5
            ym =  ym + ((yi*info.map_info.dy)-info.map_info.dy);
        end
        if xi > 1.5
            xm = xm - ((xi*info.map_info.dy)-info.map_info.dx);
        end

        varargout{3} = xm + ((0:info.samples-1).*info.map_info.dx); % originally varargout{2}
        varargout{4} = fliplr(ym - ((0:info.lines-1).*info.map_info.dy)); % originally varargout{3}
    end
end
%% Set binary format parameters
switch info.byte_order
    case {0}
        machine = 'ieee-le';
    case {1}
        machine = 'ieee-be';
    otherwise
        machine = 'n';
end
switch info.data_type
    case {1}
        format = 'int8';
    case {2}
        format= 'int16';
    case{3}
        format= 'int32';
    case {4}
        format= 'float';
    case {5}
        format= 'double';
    case {6}
        disp('>> Sorry, Complex (2x32 bits)data currently not supported');
        disp('>> Importing as double-precision instead');
        format= 'double';
case {9}
        error('Sorry, double-precision complex (2x64 bits) data currently not supported');
case {12}
        format= 'uint16';
case {13}
        format= 'uint32';
case {14}
         format= 'int64';
case {15}
        format= 'uint64';
    otherwise
        error(['File type number: ',num2str(dtype),' not supported']);
end
%% Read File

fid=fopen(file);
% read the whole file 
format=strcat(format,'=>',format);

Z = fread(fid,info.samples*info.lines*info.bands,format,0,machine); fclose all;

switch lower(info.interleave)
    case {'bsq'}
        % reshapes into samples x lines x bands
        Z = reshape(Z,[info.samples,info.lines,info.bands]);
        for k = 1:info.bands;
            %tmp(:,:,k) = rot90(Z(:,:,k)); %original
            tmp(:,:,k) = (Z(:,:,k))';  % with this transformation use imagesc(Z) to visualize

        end
    case {'bil'}
        Z = reshape(Z,[info.samples,info.lines*info.bands]);
        for k=1:info.bands
            %tmp(:,:,k) = rot90(Z(:,k:info.bands:end)); % original
            tmp(:,:,k) = Z(:,k:info.bands:end)';  % with this transformation use imagesc(Z) to visualize
        end
    case {'bip'}
        tmp = zeros(info.lines,info.samples);
        for k=1:info.bands
            tmp1 = Z(k:info.bands:end);
            tmp(:,:,k) = rot90(reshape(tmp1,[info.samples,info.lines]));
        end
end
Z = tmp;

varargout{1} = Z;
varargout{2} = info; % originally varargout{4}
