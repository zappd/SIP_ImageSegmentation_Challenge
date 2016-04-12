function hsi2gif( file_path, hdr_path )

dest = [file_path '.gif'];
raw = enviread(file_path, hdr_path);
depth = size(raw, 3);


%% Convert To GIF
min = double(Min3d(raw));
max = double(Max3d(raw));

for i = 1 : depth
    I = 1024*mat2gray(raw(:,:,i), [min max])'; 
    if i == 1;
        imwrite(I, dest,'gif', 'DelayTime', 0.1, 'Loopcount', inf);
    else
        imwrite(I, dest,'gif', 'DelayTime', 0.1, 'WriteMode', 'append');
    end
end


