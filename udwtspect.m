function [p,frequency,scale,dp,swa,swd] = udwtspect(y,dt,nlevel,dbworder)

% [p,scale,dp] = udwtspect(y,dt,wname,displ);
%
% UDWTSPECT estimates the spectral power density of a time
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

if nargin<2, dt = 1; end
if nargin<3, nlevel=5; end
if nargin<4, dbworder = 4; end

%% Set parameters needed for UDWT

samplelength=length(y);

% If length of data is odd, turn into even numbered sample by getting rid 
% of one point

if mod(samplelength,2)>0
    y=y(1:end-1);
end

samplelength=length(y);

% edge extension mode set to periodic extension by default with this
% routine in the rice toolbox.

pads=(2^(nextpow2(samplelength))) - samplelength;   % for edge extension

[h_0,h_1] = daubcqf(2*dbworder);    % set up the filters
% [h_0,~,~,~] = wfilters('coif2');

%% Do the SWT decompositon

% Gets the data size up to the next power of 2 due to SWT restrictions
% Although periodic extension is used for the wavelet edge handling we are
% getting the data up to the next power of 2 here by extending the data
% sample with a constant value
    
    y=[y(1).*ones(pads/2,1); y; y(end).*ones(pads/2,1)];
    
% Decompose the signal using the SWT

tic

[swa,swd,L] = mrdwt(y,h_0,nlevel);

toc

swa=swa';
swd=swd';

%% Getting rid of the edge effects

% First get rid of the stuff which you added onto the ends

swd=swd(:,((pads/2)+1):(end-(pads/2)));
swa=swa(:,((pads/2)+1):(end-(pads/2)));

% Now get rid of the bonified edge effects

filterlength = length(h_0);

clear('h_0','h_1');

for j=1:1:nlevel
    
    extra = (2^(j-2))*filterlength; % give some reasoning for this eq
    
    swd(j,1:extra)=NaN;
    swd(j,end-extra+1:end)=NaN;
        
end

    swa(1,1:extra)=NaN;
    swa(1,end-extra+1:end)=NaN;

clear('j','filterlength','extra');

%% Calculating power spectral density

p = zeros(nlevel,1);
dp = zeros(nlevel,1);

pars=0.0;

for k = 1:nlevel

    cw = swd(k,:).*swd(k,:);
    p(k) = nanmean(cw);
    dp(k) = nanstd(cw);
    
end

scale = (2.^[1:nlevel]');
frequency = scal2frq(scale,strcat('db',num2str(dbworder)),dt); 

% Converting to power/energy spectral density
p=p.*(dt*2);

% Calculation of error
dp=dp.*(dt*2);  % turning standard deviation into PSD units
dp=dp./(sqrt(length(y)));  % turning into standard error of the mean
dp=dp./(p.*log(10));    % propagation of error after taking log operation
dp=dp.*2;   % 95% confidence interval within two standard deviations


    errorbar(log10(frequency),log10(p),dp,'r'); hold all;
    plot(log10(frequency),log10(p),'b.-');
    grid on
    xlabel('log_{10} frequency (Hz)')
    ylabel('log_{10} PSD (signal units)^{2}Hz^{-1}')
    title(['nr of levels = ',int2str(nlevel),...
        '   Daubechies wavelet is db',num2str(dbworder)])

%% Parsevals relation check
% This is not working properly yet as MATLAB's SWT routine does not spit
% out the approximation coefficients properly, although it does allow good
% reconstruction of approximations. Details, however, work perfectly fine.

answerreal=(sumsqr(y))/(length(y));  % Average energy in original signal.

% Energy from the detail coeficients

iu=0.0;

for i=1:1:nlevel
    iu=iu + (p(i)/(2^(i+1)));
end

iu=iu*(1/dt);   % Remember that sampling frequency = (1/dt)

% Energy from the approximation coeficients

ip=1/2;

for i=1:1:nlevel
    ip=ip-(1/(2^(i+1)));
end

ip=ip*(nanmean(swa(1,:).*swa(1,:)))*(dt*2);

ip=ip*(1/dt);

answerfreq=(ip+iu);

parsevalcheck=answerreal/answerfreq   