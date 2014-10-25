%% function [p,scale,dp] = wspect(y,dt,wname)

% [p,scale,dp] = wspect(y,dt,wname,displ);
%
% WSPECT estimates the spectral power density of a time
% series using the discrete wavelet decomposition. The power
% density is estimated at scales that are distributed on a dyadic grid
%   y       time series [array]
%   dt      sampling period, default is 1
%   wname   wavelet name, in Matlab syntax. Default is 'db4' for 4th
%           order Daubechies
%   displ   Display results if displ>0. Default is with display
%   p       power density : [n,1] array
%   scale   corresponding scale axis : [n,1] array
%   dp      standard deviation of P estimated from variability
%                                           ThDdW 6/07
%% load data

clear all; close all;

cd ~/khuromOLD/SOLARStuff/NewClusterWork/Data/BurstModeIntervals/20070130_0000
load BMerged.mat
y=BTot(:,4); 

dt=1/450; wname='db2'; displ=1; 

%% random argument crap
% if nargin<2, dt = 1; end
% if nargin<3, wname = 'db4'; end
% if nargin<4, displ = 1;     end


%% sets the maximum level for the swt
npts = length(y);
nlevel = floor(log(npts)/log(2))      % maximum nr of levels

L = wmaxlev(npts,wname) 

nlevel=L

%% Do the SWT decompositon
% The periodic extension of this needs to be sorted out. A temporary
% solution is to just chop 100 points off each end.

tic

WT = ndwt(y,nlevel,wname);

toc

%% reconstruct ndwt approximations

tic

A = cell(1,nlevel);
D = cell(1,nlevel);

for m = 1:nlevel
    A{m} = indwt(WT,'a',m);   % Approximations (low-pass components)
    D{m} = indwt(WT,'d',m);   % Details (high-pass components)
    m
end

toc

%% ignore this
% 
% [c,l] = wavedec(y(:),nlevel,wname);    % wavelet decomposition
% 
% % do some bookkeeping with the indices
% first = cumsum(l)+1;
% first = first(end-2:-1:1);
% ld   = l(end-1:-1:2);
% last = first+ld-1;
% cf = c;
% lf = l;

%% reconstruct all the approximations at all levels

tic

mzero = zeros(size(swd));
A = mzero;
A(nlevel,:) = iswt(swa,mzero,wname);

D = mzero;
for i = 1:nlevel
    swcfs = mzero;
    swcfs(i,:) = swd(i,:);
    D(i,:) = iswt(mzero,swcfs,wname);
end

for j=nlevel-1:-1:1
    A(j,:) = A(j+1,:) + D(j+1,:);
end

toc

%% write out to mat file

MASy=[BTot(100:numel(y)-100,1)';A(:,100:end-100);swd(:,100:end-100)];

%% Construct background field directions using external function
% This function will also output the angle between the magnetic field and
% the background solar wind plasma flow velocity direction which we can
% assume to be constant.

%% Decompose detail coefficients at all levels into perp and par components

nlevel=13; dt=1/450; wname='db2'; displ=1; 

load masx2.mat
load masy2.mat
load masz2.mat

coefspar=zeros(nlevel,numel(MASx(1,:)));
coefsperp=zeros(nlevel,numel(MASx(1,:)));

%%

for i=1:1:nlevel
    
    [parDir] =anisotropy1([MASx(i+1,:)',MASy(i+1,:)',MASz(i+1,:)']);
    [parperp] =anisotropy2(parDir,[MASx(i+nlevel+1,:)',MASy(i+nlevel+1,:)',MASz(i+nlevel+1,:)']);
    
    coefspar(i,:)=parperp(:,1);
    coefsperp(i,:)=parperp(:,2);
    
    clear('parDir','parperp');
    
end

clear('i','MASx','MASy','MASz');

%% Calculating power spectral density for par with errors

p_par = zeros(nlevel,1);
dp_par = zeros(nlevel,1);

for k = 1:nlevel
    
    cw = coefspar(k,:).*coefspar(k,:);
    p_par(k) = mean(cw);    % I think your one over N is being handled over here. So I think the
    dp_par(k) = std(cw);      % SWT should be ok
    
end

p_par=p_par.*(dt*2);

clear('k')

%% Calculating power spectral density for perp with errors

p_perp = zeros(nlevel,1);
dp_perp = zeros(nlevel,1);

for k = 1:nlevel
    
    cw = coefsperp(k,:).*coefsperp(k,:);
    p_perp(k) = mean(cw);    % I think your one over N is being handled over here. So I think the
    dp_perp(k) = std(cw);      % SWT should be ok
    
end

p_perp=p_perp.*(dt*2);

clear('k')

%% Translating scale to frequency and normalising to power spectral density

scale = 2.^[1:nlevel]';
frequency = scal2frq(scale,wname,dt); 

%% plotting 
if displ
    clf
    figure(1)
    plot(log10(frequency),log10(p_par),'k.-'); hold on;
    plot(log10(frequency),log10(p_perp),'r.-');
    grid on
    xlabel('scale')
    ylabel('spectral density')
    title(['nr of levels = ',int2str(nlevel),'   wavelet is ',wname])
    
    figure(2)
    plot(log10(frequency),(p_par./p_perp),'b.-');
end
    