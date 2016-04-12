function Y=cube2mat(X,output)
% if X is a p x q x m spectral cube then 

% if output is column then Y is a matrix m x (p*q) 
% where the j+(i-1) column is the spectrum for the (i,j) voxel
% i. e. 1st row, 2nd row , ith row ...

% if output is row then Y is a matrix m x (p*q) 
% where the i+(j-1) column is the spectrum for the (i,j) voxel
% i.e. 1st column, 2nd column, jth column ...

[ p, q, m ] = size(X);
if strcmp(output,'col')
    
Y=permute(X, [3 2 1]); % need to use ipermute(Y,[3 2 1]) to undo 
else 
    Y=permute(X, [3 1 2]);
end

Y=reshape(Y, [ m p*q]);


