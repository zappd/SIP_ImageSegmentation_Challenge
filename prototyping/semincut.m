%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    CONFIG    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  clear; close all;

% Affinity function
  AFNT_TYPE='gauss';
  SIGMA=0.1; %constant between problem size

% Markov timestep
  %T=2400;  %good for 64x64
  T=600;   %good for 32x32
  %should scale with problem size

% Spectral approximation
  %SPEC_TYPE='exact'; 
  %SPEC_TYPE='both'; 
  SPEC_TYPE='approx'; 
  %EVMAXFACTOR=0.00; %all eigenvalues
  EVMAXFACTOR=0.05; 
  %EVMAXFACTOR=1/3; %suggested
  EVMAX=400; %use N eigenvalues
  %EVMAX=32; %use N eigenvalues
  %EVMAX=0;  %select eigenvalues based on EVMAXFACTOR

% Test case
  %DATA_TYPE='test'; %diffused 32z32 cube
  %DATA_TYPE='jesus'; %first mode of hyperspectra image cropped
  %DATA_TYPE='jesus_64'; %first mode of hyperspectra image cropped
  %DATA_TYPE='jesus_96'; %first mode of hyperspectra image cropped
  NORM_DATA=true;

% Report
  RES_FIG_INID=true; %Show initial data
  RES_FIG_EIGV=2;     %Show N eigenvalues
  RES_FIG_BCNT=3;    %Show blur kernels centered at N points
  RES_FIG_BTOP=false; %show blue on top of original
  RES_FIG_BTYP='linspace';    %Show blur kernels centered at N random points
  %RES_FIG_BTYP='rand';    %Show blur kernels centered at N random points
  %RES_FIG_BTYP='pick';    %pick points on figure
  RES_FIG_ZTRE=true;    %view clustering in first three z dimensions
  RES_FIG_STRE=true;    %view clustering in first three s dimensions (normalized Z) with cluster centers
  RES_FIG_CLUS=true;    %view clusters in original image
  RES_FIG_THRE=true;    %view thresholded clusters in original image
  %RES_FIG_THRW=true; %show edges per centroid connection
  RES_FIG_THRW=false; %show edges per centroid connection
  RES_FIG_TSEN=true;    %view threshold cluster sensitivity (TAU)

% Thresholds

% Kmeans config
  KMEANS_REPL=16; %number of times to replicates kmeans with kmeans++ to check for local minima
  KMEANS_ELBW=0.001; %elbow threshold
  %KMEANS_ELBW=0.05; %elbow threshold
  %KMEANS_ELBW=0.1; %elbow threshold

% Cluster Analysis
  CLUS_TAU=0.89; % Classification threshold
  CLUS_ELBW=0.5; %difference factor for elbow detection
  CLUS_NUMTEST=100; %number of points to check for elbow
  CLUS_OVRR=false; %override difference factor with CLUS_TAU

% Continuity Analysis
  CONT_ESTR=20; % Number of discretized points to test eqn 9
  CONT_FTHC=0.2; % Threshold for continuity
  CONT_FTHD=0.3; % threshold for discontinuity

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figcount=1;
rng('shuffle');

if(strcmp(DATA_TYPE,'jesus'))
  data=load('data2_wvn1');
  data=data(1:32,97:128);
elseif(strcmp(DATA_TYPE,'jesus_64'))
  data=load('data2_wvn1');
  data=data(1:64,65:128);
elseif(strcmp(DATA_TYPE,'jesus_96'))
  data=load('data2_wvn1');
  data=data(33:128,33:128);
elseif(strcmp(DATA_TYPE,'test'))
  data=ones(10,10);
  data=padarray(data,[1 1],0.9);
  data=padarray(data,[1 1],0.8);
  data=padarray(data,[1 1],0.8);
  data=padarray(data,[1 1],0.5);
  data=padarray(data,[1 1],0.2);
  data=padarray(data,[1 1],0.1);
  data=padarray(data,[1 1],0.1);
  data=padarray(data,[1 1],0.051);
  data=padarray(data,[1 1],0.111);
  data=padarray(data,[1 1],0.661);
  data=padarray(data,[1 1],0.661);
else
  fprintf('Config option invalid: %s=%s','DATA_TYPE',DATA_TYPE);
  return;
end

if(NORM_DATA==true)
  data=data-mean(reshape(data,[size(data,1)*size(data,2) 1]));
  stddev=sqrt(var(reshape(data,[size(data,1)*size(data,2) 1]),1));
  data=data/stddev;
end

if(RES_FIG_INID==true)
  ogdatafig=figure(figcount);
  imagesc(data');
  axis xy;
  title('Original data');
  colormap('gray');
  figcount=figcount+1;
end

% Add padding for boundary conditions
data=padarray(data,[1 1]);

N=(size(data,1)*size(data,2));
I=reshape(data,[N 1]);
A=zeros(N,N);
        
lx=size(data, 1);
ly=size(data, 2);
sz=size(data);

for i=2:lx-1
for j=2:ly-1
  idx_i=sub2ind(sz,i,j);

  ih=i;
  jh=j;
  A(idx_i,sub2ind(sz,i+1,j+1))=exp(-(data(ih,jh)-data(ih+1,jh+1))^2/(2*SIGMA^2));
  A(idx_i,sub2ind(sz,i  ,j+1))=exp(-(data(ih,jh)-data(ih,  jh+1))^2/(2*SIGMA^2));
  A(idx_i,sub2ind(sz,i-1,j+1))=exp(-(data(ih,jh)-data(ih-1,jh+1))^2/(2*SIGMA^2));
  
  A(idx_i,sub2ind(sz,i+1,j  ))=exp(-(data(ih,jh)-data(ih+1,jh  ))^2/(2*SIGMA^2));
  A(idx_i,sub2ind(sz,i  ,j  ))=exp(-(data(ih,jh)-data(ih  ,jh  ))^2/(2*SIGMA^2)); %unity
  A(idx_i,sub2ind(sz,i-1,j  ))=exp(-(data(ih,jh)-data(ih-1,jh  ))^2/(2*SIGMA^2));

  A(idx_i,sub2ind(sz,i+1,j-1))=exp(-(data(ih,jh)-data(ih+1,jh-1))^2/(2*SIGMA^2));
  A(idx_i,sub2ind(sz,i  ,j-1))=exp(-(data(ih,jh)-data(ih,  jh-1))^2/(2*SIGMA^2));
  A(idx_i,sub2ind(sz,i-1,j-1))=exp(-(data(ih,jh)-data(ih-1,jh-1))^2/(2*SIGMA^2));
end
end

% Toss out halo cells
A(1:lx             ,:)=[];
A((end-lx):(end),:)=[];

A(:            ,1:ly)=[];
A(:,(end-lx):(end))=[];

k=1;
for j=1:size(A,2)
  if(A(j,j)==0)
    idx(k)=j;
    k=k+1;
  end
end

A(:,[idx])=[];
A([idx],:)=[];

data(:,1)=[];
data(:,end)=[];
data(1,:)=[];
data(end,:)=[];

% Reset data stats without halos
N=(size(data,1)*size(data,2));
lx=size(data, 1);
ly=size(data, 2);
[sz]=size(data);

% Set normalizing factor so rows sum to unity
D=sum(A,2);

%Dont need to form n generally
M=zeros(N,N);
for j=1:N
  M(:,j)=A(:,j)./D(j);
end

if( strcmp(SPEC_TYPE,'exact') || strcmp(SPEC_TYPE,'both'))
end

%In production code we would get the eigenvalues through decomposition 
if( strcmp(SPEC_TYPE,'approx') || strcmp(SPEC_TYPE,'both'))
  
  if(EVMAXFACTOR==0.00)
    startcount=N-2;
    endcount=N-2;
  else
    startcount=16;
    endcount=N-2;
  end  

  if(EVMAX>0)
    startcount=EVMAX;
    endcount=EVMAX;
  end

opts.tol = 1e-4; 

    for lamcount=startcount:4:endcount
      [eVecs eVals]=eigs(M,lamcount,'lr',opts);
      etmp=diag(eVals);
      if( etmp(end)^(T) <= EVMAXFACTOR )
        break;
      end
    end
  
  Dmat=diag(D);
  
  for i=1:lamcount;
      tmp=(eVals^T)*eVecs';
      wti(:,i)=tmp(i,:)*D(i).^(-1/2);
  end
  
  qti=(Dmat.^(1/2)*eVecs)*wti';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%Visualize blur kernel
if(RES_FIG_BCNT>0)
  if(RES_FIG_BTOP==true)
    tfig=ogdatafig;
    figure(tfig)
  else
    tfig=figcount;
    figure(tfig)
    title('approximate markov matrix');
    figcount=figcount+1;
  end
  
  if(strcmp(RES_FIG_BTYP,'rand'))
    rpnts=squeeze(randi(N,RES_FIG_BCNT,1));
  elseif(strcmp(RES_FIG_BTYP,'pick'))
    figure(ogdatafig);
    for i=1:RES_FIG_BCNT
      while(1)
        [xx(i) yy(i)  ] = ginput(1);
        break;
      end
    end
    rpnts=sub2ind(sz,round(xx),round(yy));
    figure(tfig);
  else
    xx=round(linspace(1,lx,RES_FIG_BCNT));
    yy=round(linspace(1,ly,RES_FIG_BCNT));
    rpnts=sub2ind(sz,xx,yy);
  end
  
    color=zeros(3,lx,ly,3);
    pcolor=zeros(3,3);
    color(1,:,:,:) = cat(3, zeros(size(data)), ones(size(data)), zeros(size(data)));
    color(2,:,:,:) = cat(3, ones(size(data)), zeros(size(data)), zeros(size(data)));
    color(3,:,:,:) = cat(3, zeros(size(data)), zeros(size(data)), ones(size(data)));
    pcolor(1,:) = [ 0.3 1 0.3];
    pcolor(2,:) = [ 1 0.3 0.3];
    pcolor(3,:) = [ 0.3 0.3 1];
  
  
  if(strcmp(SPEC_TYPE,'approx') || strcmp(SPEC_TYPE,'both'))
  hold on
  for i=1:RES_FIG_BCNT
    itcol=mod(i,3)+1;
  
    datav=reshape(data,[N 1])'.*qti(rpnts(i),:);
  
    h(i) = imagesc(squeeze(color(itcol,:,:,:))); 
  
    [xx yy]=ind2sub(sz,rpnts(i));
    scatter(xx,yy,'MarkerEdgeColor','k','MarkerFaceColor',pcolor(itcol,:));
  
    datav=datav+min(min(datav));
    datav=datav/max(max(datav));
    set(h(i), 'AlphaData', reshape(datav,size(data))'/max(max(datav)));
  end
  hold off
  end
  
  if(strcmp(SPEC_TYPE,'exact') || strcmp(SPEC_TYPE,'both'))
    tfig=figcount;
      figure(tfig)
    hold on
    if(RES_FIG_BTOP==true)
      imagesc(data);
      title('Original Data, exact markov matrix');
      else
      title('exact markov matrix');
    end
    for i=1:RES_FIG_BCNT
      itcol=mod(i,3)+1;
    
      MM=M^T;
      datav=reshape(data,[N 1])'.*MM(rpnts(i),:);
    
      g(i) = imagesc(squeeze(color(itcol,:,:,:))); 
    
      [xx yy]=ind2sub(sz,rpnts(i));
      scatter(xx,yy,'MarkerEdgeColor','k','MarkerFaceColor',pcolor(itcol,:));
    
      datav=datav+min(min(datav));
      datav=datav/max(max(datav));
      set(g(i), 'AlphaData', reshape(datav,size(data))'/max(max(datav)));
    end
    hold off
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(RES_FIG_ZTRE==true)
  Q=eVecs'*Dmat*eVecs;
  z=Q^(1/2)*wti';
  h=figure();
  scatter3(z(1,:),z(2,:),z(3,:));
end


s=zeros(lamcount,N);
%S is normalized to be on the unit sphere
for i=1:size(z,2)
S(:,i)=z(:,i)/norm(z(:,i));
end


%find K mean clusters
Kmax=round(sqrt(N));


[junka,junkb,sum0]=kmeans(S',1,'Start','plus','Replicates',KMEANS_REPL);


%find elbow approx
oldsumd=sum(sum0);
sumlist(1)=oldsumd;
for K=2:2:Kmax
%idx gives labels by cluster index, C gives centroid locations, sumd gives sum of centroid distances
  [idxlabel, C, sumd, distm]=kmeans(S',K,'Start','plus','Replicates',KMEANS_REPL);;
  sumlist(K)=sum(sumd);
  dsumd=sum(sumd)-oldsumd;
  oldsumd=sum(sumd);
  if(abs(dsumd/sum0)<KMEANS_ELBW)
    break
  end
end

if (K==Kmax)
  display('Failed to find good elbow, change Kmax');
end

%plot clusters in S space
if(RES_FIG_STRE==true)
  h=figure();
  scatter(S(2,:),S(3,:));

  hold on
  scatter(C(:,2),C(:,3),'r');
  hold off
end

%plot clusters on image
if(RES_FIG_CLUS==true)
  figure();
  imagesc(data');
  axis xy;
  title('Original data');
  colormap('gray');

  hold on;

  CM=jet(K);
  
  for i=1:N
    [xx yy]=ind2sub(sz,i);
    scatter(xx,yy,25,CM(idxlabel(i),:));
  end

  figcount=figcount+1;
end

%%%%% WORKING SECTION

% find threshold knee
dotMet=S'*C';

if CLUS_OVRR==false
  testthresh=linspace(0.0,1,100);
  sumdefined=zeros(100,1);
  cluster_max=zeros(N,1);
  for i=1:N
    cluster_max(i)=dotMet(i,idxlabel(i));
  end
  
  %%%
  for i=1:100
    sumdefined(i)=1-sum((cluster_max<testthresh(i)))/sum((cluster_max<1.01));
  end
  
  dsumdefined=diff(sumdefined);
  for i=1:100
    if(sumdefined(i)<CLUS_ELBW)
      CLUS_TAU=testthresh(i);
      break;
    end
  end

  
  if RES_FIG_TSEN==true
    figure();
    plot(testthresh,sumdefined);
    xlabel('\tau threshold');
    ylabel('union cardinality');
    hold on
    plot([testthresh(i) testthresh(i)],[0 1],'r');
    hold off
  end

else
  %let CLUS__TAU be defined
  
end

%%%%%


% Threshold clusters
for i=1:N

  label=idxlabel(i);

  %dotMet=S'*C';

  if( dotMet(i,label) < CLUS_TAU)
    threshlabel(i)=0;
  else
    threshlabel(i)=label;
  end

  Kthresh=length(unique(threshlabel))-1;
end

[junk,activeK]=ismember(1:K,unique(threshlabel));



%calculate center of mass of sections in image space
  xxtmp=zeros(N,1);;
  yytmp=zeros(N,1);;

for i=1:K
k=0;
  for j=1:N
    if(threshlabel(j)==i)
      [xx yy]=ind2sub(sz,j);
      k=k+1;
      xxtmp(k)=xx;
      yytmp(k)=yy;
    end
  end
  XC(i)=mean(xxtmp(1:k));
  YC(i)=mean(yytmp(1:k));
end


%calculate continuity between nodes
r=linspace(0,1,CONT_ESTR+1);

for k=1:K
  for j=1:K
    for ri=1:CONT_ESTR+1
      n=r(ri)*C(k,:)+(1-r(ri))*C(j,:);
      m=n/norm(n);
      
      tmpDot=S'*m';
      tmpDot=tmpDot>CLUS_TAU;
      card=sum(tmpDot);
      
      ftmp(ri)=card/(sum(threshlabel==k)*r(ri)+sum(threshlabel==j)*(1-r(ri)));
    end
    [f(k,j),farg(k,j)]=min(ftmp);
  
  end
end

%plot thresholded clusters
if(RES_FIG_THRE==true)
  figure();
  imagesc(data');
  axis xy;
  title('Original data');
  colormap('gray');

  hold on;

  CM=jet(K);
  
  for i=1:N
    if(threshlabel(i)==0)
      continue;
    end
    [xx yy]=ind2sub(sz,i);
    scatter(xx,yy,25,CM(threshlabel(i),:));
  end

  scatter(XC,YC,55,'rx');

  if(RES_FIG_THRW==true)
  CM=jet(10);
  RM=linspace(1,0,10);
  GM=linspace(0,1,10);
  BM=linspace(0,0,10);
  CM=[ RM' GM' BM' ];


  for i=1:K
  for j=1:K

      indexc=floor(10*f(i,j))+1;
      if (isnan(indexc))
        continue;
      end
      if (indexc>10)
        indexc=10;
      end
        
      p=plot( [ XC(i) XC(j) ], [YC(i) YC(j)] );
  
      set(p,'Color', CM(indexc,:));

  end
  end
  end

  figcount=figcount+1;
end






%for i=1:idxlabel
%
%end


%S'*C'-tau0











