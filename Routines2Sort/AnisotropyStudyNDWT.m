%% m-file for project: Solar wind turbulence anisotropy study
%**************************************************************************
% 
% File name: AnisotropyStudy
% Author: K. H. Kiyani
% File created: Sometime in March 2010
% Last updated: 13th May 2010
% 
% Description: (needs to be redone)
% ------------
% This m-file estimates the spectral power density of a time
% series using the discrete wavelet decomposition. The power
% density is estimated at scales that are distributed on a dyadic grid
%   y       time series [array]
%   dt      sampling period, default is 1
%   wname   wavelet name, in Matlab syntax. Default is 'db4' for 4th
%           order Daubechies
%   p       power density : [n,1] array
%   scale   corresponding scale axis : [n,1] array
%   dp      standard deviation of P estimated from variability
%
%**************************************************************************
% Notes:
% Remember that when you are calculating structure functions you will have
% to normalise the coefficients by square root of the scale you are looking
% at. This is because the wavelets coefficients are normalised w.r.t. this
% factor due to the L^2 Norm. This does not effect calculation of the power
% as the power is related to the L^2 norm.
% 
% Example code to plot structure functions:
% plot(log10(scale),log10((p./(scale))./(2*dt)),'.r')
% 
% Boundary effects on finite size samples on the interval:
% Check the Advanced concepts in MATLAB help or where to truncate at each
% level appropriately.
% Example code when you are checking approximations and boundary effects:
% plot(y); hold on; plot(A(16,:),'r','LineWidth',2);
% 
% We need to create seperate functions for structure functions, PSD, PDF
% based on different methods, and then plot them all together using
% subplot on matlab; and the output them to EPS files then going onto a
% LaTeX file. 
% 
%% load data

clear all; close all;

cd ~/khuromOLD/SOLARStuff/NewClusterWork/Data/BurstModeIntervals/20070130_0000

load BMerged.mat

samplelength=length(BTot(:,1));

% If length of data is odd, turn into even numbered sample by getting rid 
% of one point

if mod(samplelength,2)>0
    BTot=BTot(1:end-1,:);
end

samplelength=length(BTot(:,1));

dt=1/450; wname='db8';

nlevel = wmaxlev(samplelength,wname)   % maximum nr of levels

dwtmode('per'); % edge extension mode set to periodic extension

%% Do the SWT decompositon

% Gets the data size up to the next power of 2 due to SWT restrictions
% You might want to periodic extension here instead of a constant

pads=(2^(nextpow2(samplelength))) - samplelength;

SWT_parts = cell(1,4);

SWT_parts{1,1} = BTot(:,1);

%%

for m=2:1:2
    
    y=[BTot(1,m).*ones(pads/2,1); BTot(:,m); BTot(end,m).*ones(pads/2,1)];
    
    tic
    [swa,swd] = swt(y,nlevel,wname);
    toc
    
% reconstruct all the approximations and details at all levels

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
    
%     now get rid of the edge effects

    SWT_parts{1,m} = cell(1,3);
    
    SWT_parts{1,m}{1,1} = A(:,((pads/2)+1):(end-(pads/2)));
    SWT_parts{1,m}{1,2} = D(:,((pads/2)+1):(end-(pads/2)));
    SWT_parts{1,m}{1,3} = swd(:,((pads/2)+1):(end-(pads/2)));

    clear('y','swa','swd','A','D');

end

%%

clear('m','pads');

%% write out to mat file ignoring boundary effects (although I think I
%% should keep these)

MASy=[BTot(100:numel(y)-100,1)';A(:,100:end-100);swd(:,100:end-100)];

%% Construct background field directions using external function
% This function will also output the angle between the magnetic field and
% the background solar wind plasma flow velocity direction which we can
% assume to be constant.
%% Good

%% Decompose detail coefficients at all levels into perp and par components

nlevel=13; dt=1/450; wname='db2'; displ=1; 

load masx2.mat
load masy2.mat
load masz2.mat

% This si nto a biased towards the x-component, it's only to set the size
% of the array
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

p_perp=p_perp.*(dt*2);  % normalising to power spectral density

clear('k')

%% Translating scale to frequency

scale = 2.^[1:nlevel]';
frequency = scal2frq(scale,wname,dt); 

%% plotting 
    figure(1)
    
    plot(log10(frequency),log10(p_par),'k.-'); hold on;
    plot(log10(frequency),log10(p_perp),'r.-');
    
    grid on
    
    xlabel('frequency')
    ylabel('Power spectral density')
    title(['nr of levels = ',int2str(nlevel),'   wavelet is ',wname])
    
    figure(2)
    plot(log10(frequency),(p_par./p_perp),'b.-');
    

% Things to discuss in paper:
% ***************************
% Why use the wavelet strategy that you are doing: You will need to discuss
% orthogonal wavelets, perfect reconstructions, compact nature of wavelets,
% DWT strategy provides scaling funtions which cover the lower frequecies
% and nothing is thrown away. How you handled boundary effects at each
% scale. Why did you use Db6, mention tap filters. Might want to discuss
% conditional averages and Reynolds decomposition. Find out which one you
% are using, stationary wavelet transform, undecimated wavelet transform
% etc. and explain why you are using them i.e. to line up fetaures with
% the original signal. 