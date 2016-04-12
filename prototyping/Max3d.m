 function [Max,index] = Max3d(x)
 N=size(x);
 if length(N)<3
    error(' Input 3D array');
 end
 [a,t]=max(x(:));
 Max=a;
 index1=ceil(t/(N(1)*N(2)));
 %Now we find the slice that contains the Max
 Temp=x(:,:,index1);
 [index2,index3]=find(Temp==max(Temp(:)));
 index=[index2;index3;index1]';