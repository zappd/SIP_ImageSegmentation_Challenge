function [Z,M]=matMask2cube(Y,mask,rowN,colN,input)

% if Y is a m x (p*q) matrix then 

% if input is column then Z is a cube p x q x m


p=rowN;
q=colN;
m=size(Y,1);
if rem(size(Y,2),p*q)~=0
    error('matrix dimensions not congruend with input parameters')
end
if strcmp(input,'col')
Z=reshape(Y, [ m q p]);    
Z=ipermute(Z, [3 2 1]); % need to use ipermute(Y,[3 2 1]) to undo 
M=(reshape(mask,colN,rowN))';
else
    Z=reshape(Y,[m p q]);
    Z=ipermute(Z, [3 1 2]);
    M=reshape(mask',rowN,colN);
end

