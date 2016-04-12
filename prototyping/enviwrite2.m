function i=enviwrite2(Z,fname,info)

% enviwrite          	- write ENVI image from MATLAB array (V. Guissard, Apr 29 2004)
%
% 				Write a MATLAB array to a file in ENVI standard format
%				from a [col x line x band] array
%
% SYNTAX
%
% image=freadenvi(fname)
% [image,p]=freadenvi(fname)
% [image,p,t]=freadenvi(fname)
%
% INPUT :
%
%
% image	c by l by b	name of the MATLAB variable containing the array to export
%				to an ENVI image, with c = cols, l the lines and b the bands
% fname	string	full pathname of the ENVI image to write.
%
% OUTPUT :
%
% i		integer	i = -1 if process fail
%
% NOTE :
%
%%%%%%%%%%%%%

% Parameters initialization

elements={'samples =' 'lines   =' 'bands   =' 'data type ='};
%d=[4 1 2 3 12 13];
% Check user input
if ~ischar(fname)
    error('fname should be a char string');
end




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
        format= 'float'; %maybe necessary float32
        format2='single';
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

% opening the file for writing

wfid = fopen(fname,'w');
if wfid == -1
    i=-1;
else i=40;
end
disp('Writing ENVI image ...');





switch lower(info.interleave)
    case {'bsq'}
        tmp=zeros([info.samples,info.lines,info.bands],format2);
        for k = 1:info.bands;

            tmp(:,:,k) = (Z(:,:,k))';
        end

        % reshapes into samples x lines x bands
        tmp = reshape(tmp,info.samples*info.lines*info.bands,1);

    case {'bil'}

        tmp=zeros([info.samples,info.lines*info.bands],format2);

        for k=1:info.bands
            tmp(:,k:info.bands:end)=squeeze(Z(:,:,k)');
        end

        tmp = reshape(tmp,info.samples*info.lines*info.bands,1);

    case {'bip'}
        disp('not yet implemented')

end
Z=tmp;
clear tmp
%format=strcat(format,'=>',format);
fwrite(wfid,Z,format,0,machine);


fclose(wfid);

% Write header file

fid = fopen(strcat(fname,'.hdr'),'w');
if fid == -1
    i=-1;
else i=40
end

fprintf(fid,'%s \n','ENVI');
fprintf(fid,'%s \n','description = {');
fprintf(fid,'%s \n','Exported from MATLAB}');
fprintf(fid,'%s %i \n',elements{1,1},info.samples);
fprintf(fid,'%s %i \n',elements{1,2},info.lines);
fprintf(fid,'%s %i \n',elements{1,3},info.bands);
fprintf(fid,'%s %i \n',elements{1,4},info.data_type);
fprintf(fid,'%s %s \n','interleave = ',info.interleave);
fclose(fid);