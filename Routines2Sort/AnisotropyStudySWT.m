%**************************************************************************
% m-file for project: Solar wind turbulence anisotropy study
%**************************************************************************
% 
% File name: AnisotropyStudySWT.m
% Author: K. H. Kiyani
% File created: March 2010
% Last updated: 17th May 2010
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
% to normalise the coefficients by 'square root of the scale' you are 
% looking at. This is because the wavelets coefficients are normalised 
% w.r.t. this factor due to the L^2 Norm. This does not effect calculation 
% of the power as the power is related to the L^2 norm.
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

c_load VCh1;

vel=mean(VCh1(:,2:4),1);

vel=vel./norm(vel);

clear('VCh1');

load BMerged.mat

%% Set parameters needed for SWT

samplelength=length(BTot(:,1));

% If length of data is odd, turn into even numbered sample by getting rid 
% of one point

if mod(samplelength,2)>0
    BTot=BTot(1:end-1,:);
end

samplelength=length(BTot(:,1));

dt=1/450; wname='db4';  

nlevel = wmaxlev(samplelength,wname)   % maximum nr of levels

nlevel=15; % Just a precaution for higher levels that take a long time. 

dwtmode('per'); % edge extension mode set to periodic extension

pads=(2^(nextpow2(samplelength))) - samplelength;   % for edge extension

%% Do the SWT decompositon and reconstruction

tempfiles = cell(1,3);

for m=2:1:4
    
% Gets the data size up to the next power of 2 due to SWT restrictions
% Although periodic extension is used for the wavelet edge handling we are
% getting the data up to the next power of 2 here by extending the data
% sample with a constant value
    
    y=[BTot(1,m).*ones(pads/2,1); BTot(:,m); BTot(end,m).*ones(pads/2,1)];
    
% Decompose the signal using the SWT
    
    tic
    [swa,swd] = swt(y,nlevel,wname);
    toc
    
    clear('y');
    
% Reconstruct all the approximations and details at all levels

    tic
    
    mzero = zeros(size(swd));
    A = mzero;
    A(nlevel,:) = iswt(swa,mzero,wname);
    
%     A = iswt(swa,mzero,wname);    % only use this if you run out of RAM
    
    clear('swa');
    
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
    
    tempfiles{m-1}=strcat('tempfile_',num2str(m),'_',wname,'.mat');
    save(tempfiles{m-1},'A','D','swd');
    clear('A','D','swd');
    
end
    
clear('m');

%% temporary for convenience if you had to start up due to restart

for m=2:1:4
    tempfiles{m-1}=strcat('tempfile_',num2str(m),'_',wname,'.mat');
end

clear('m');

%% Put all the files together into a cell structure

tic

SWT_parts = cell(1,4);

SWT_parts{1,1} = BTot(:,1);

for m=2:1:4
    
    load(tempfiles{m-1},'A','swd');
    
    SWT_parts{1,m} = cell(1,2);
    
%     only use this if you run out of RAM
%     SWT_parts{1,m}{1,1} = A(((pads/2)+1):(end-(pads/2)));

    SWT_parts{1,m}{1,1} = A(:,((pads/2)+1):(end-(pads/2)));
        clear('A');
        
%     SWT_parts{1,m}{1,2} = D(:,((pads/2)+1):(end-(pads/2)));
%         clear('D');

    SWT_parts{1,m}{1,2} = swd(:,((pads/2)+1):(end-(pads/2)));
        clear('swd');
        
end

toc

clear('m','pads','tempfiles','BTot');

%% Getting rid of the edge effects; to keep edges skip this section

[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(wname);

filterlength = length(Lo_D);

clear('Lo_D','Hi_D','Lo_R','Hi_R');

for j=1:1:nlevel
    
    extra = (2^(j-2))*filterlength; % give some reasoning for this eq
    
    for m=2:1:4
        
        % for approximations
        SWT_parts{1,m}{1,1}(j,1:extra)=NaN;
        SWT_parts{1,m}{1,1}(j,end-extra+1:end)=NaN;
        
        % for details
        SWT_parts{1,m}{1,2}(j,1:extra)=NaN;
        SWT_parts{1,m}{1,2}(j,end-extra+1:end)=NaN;
        
    end
    
    clear('m','extra');
    
end

clear('j','filterlength');

%% Decompose detail coefficients at all levels into perp and par components

coefspar=cell(1,nlevel);
coefsperpA=cell(1,nlevel);
coefsperpB=cell(1,nlevel);
coefsperpMag=cell(1,nlevel);

% coefsBVangle=cell(1,nlevel);

tic

for i=1:1:nlevel
    
% Using an external function, construct orthonormal basis comprising of 
% background field direction, and two perpendicular vectors using the solar
% wind plasma flow velocity direction.
% This function will also output the angle between the magnetic field and
% the background solar wind plasma flow velocity direction which we can
% assume to be constant.

    [parDir,perpA,perpB,~] = anisotropy1([SWT_parts{1,2}{1,1}(i,:)',...
        SWT_parts{1,3}{1,1}(i,:)',SWT_parts{1,4}{1,1}(i,:)'],vel);
    
% Project fluctuations given by wavelet details, parallel and perp to the
% scale dependent background field direction given by 'parDir'
    
    [parperp] = anisotropy2(parDir,perpA,perpB,...
        [SWT_parts{1,2}{1,2}(i,:)',...
        SWT_parts{1,3}{1,2}(i,:)',SWT_parts{1,4}{1,2}(i,:)']);
    
    clear('parDir','perpA','perpB');
    
    coefspar{1,i}=parperp(:,1);
    coefsperpA{1,i}=parperp(:,2);
    coefsperpB{1,i}=parperp(:,3);
    coefsperpMag{1,i}=parperp(:,4);
    
%     coefsBVangle{1,i}=BVangle;
    
    clear('parperp');
    
end

toc

clear('i','SWT_parts');

%% Calculating power spectral density for par with errors

p_par = zeros(nlevel,1);
dp_par = zeros(nlevel,1);

for k = 1:nlevel
    
    cw = coefspar{1,k}(:,1).*coefspar{1,k}(:,1);
    p_par(k) = nanmean(cw);    % I think your one over N is being handled over
    dp_par(k) = nanstd(cw);    %  here. So I think the SWT should be ok.
    
end

p_par=p_par.*(dt*2);

clear('k','cw');

%% Calculating power spectral density for perpA with errors

p_perpA = zeros(nlevel,1);
dp_perpA = zeros(nlevel,1);

for k = 1:nlevel
    
    cw = coefsperpA{1,k}(:,1).*coefsperpA{1,k}(:,1);
    p_perpA(k) = nanmean(cw);    % I think your one over N is being handled
    dp_perpA(k) = nanstd(cw);    % over here. So I think the SWT should be ok.
    
end

p_perpA=p_perpA.*(dt*2);  % normalising to power spectral density

clear('k','cw');

%% Calculating power spectral density for perpB with errors

p_perpB = zeros(nlevel,1);
dp_perpB = zeros(nlevel,1);

for k = 1:nlevel
    
    cw = coefsperpB{1,k}(:,1).*coefsperpB{1,k}(:,1);
    p_perpB(k) = nanmean(cw);    % I think your one over N is being handled
    dp_perpB(k) = nanstd(cw);    % over here. So I think the SWT should be ok.
    
end

p_perpB=p_perpB.*(dt*2);  % normalising to power spectral density

clear('k','cw');

%% Calculating power spectral density for perpMag with errors

p_perpMag = zeros(nlevel,1);
dp_perpMag = zeros(nlevel,1);

for k = 1:nlevel
    
    cw = coefsperpMag{1,k}(:,1).*coefsperpMag{1,k}(:,1);
    p_perpMag(k) = nanmean(cw);    % I think your one over N is being handled
    dp_perpMag(k) = nanstd(cw);    % over here. So I think the SWT should be ok.
    
end

p_perpMag=p_perpMag.*(dt*2);  % normalising to power spectral density

clear('k','cw');

%% Translating scale to frequency

scale = 2.^[1:nlevel]';
frequency = scal2frq(scale,wname,dt); 

%% plotting PSD
    figure(1)
    
    plot(log10(frequency),log10(p_par),'.-'); hold all;
    plot(log10(frequency),log10(p_perpA),'.-');
    plot(log10(frequency),log10(p_perpB),'.-');
    plot(log10(frequency),log10(p_perpMag),'k.-');
    
    grid on
    
    xlabel('frequency')
    ylabel('Power spectral density')
    title(['nr of levels = ',int2str(nlevel),'   wavelet is ',wname])
    
    figure(2)
    plot(log10(frequency),log10(p_par./p_perpMag),...
        'MarkerFaceColor',[1 1 1],...
        'MarkerEdgeColor',[0 0 0],...
        'MarkerSize',8,...
        'Marker','o',...
        'LineWidth',2,...
        'Color',[0 0 1]);
    
%% Structure functions
tic
[tausPar,AllMomentsPar]=sfunctions(coefspar,dt,scale,0);
toc
% tic
% [tausPerp,AllMomentsPerp]=sfunctions(coefsperpMag,dt,scale,1);
% toc

%% Zeta exponents

sfexponents(tausPar,AllMomentsPar)
% sfexponents(tausPerp,AllMomentsPerp)

%% PDF
    
    
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
% 
% for i=1:1:nlevel-8
%     subplot(nlevel-8,1,i)
%     plot(SWT_parts{1,1}, SWT_parts{1,4}{1,2}(i,:),'b');
% end