function [p,frequency,scale,dp,l,coefs] = wspect(y,dt,wname,fitline)

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

% notes: Need to get rid of edge effects -- This will make this routine
% perfect.
% Also note the effect from choices of wavelet bases with few zero moments.
% Notice that the Haar Wavelet overestimates the PSD especially when the
% steeper slope is encountered. This is because the Haar has only one zero
% moment. This means it can only detect PSD slopes of -3 or higher. To
% detect a power-law slope of f^-{\beta} one requires a wavelet with
% M>(\beta - 1)/2

if nargin<2, dt = 1; end
if nargin<3, wname = 'db4'; end
if nargin<4, fitline = 0; end

nlevel = wmaxlev(length(y),wname)   % maximum nr of levels

dwtmode('per'); % edge extension mode set to periodic extension

[c,l] = wavedec(y(:),nlevel,wname);    % wavelet decomposition

% do some bookkeeping with the indices
first = cumsum(l)+1;
first = first(end-2:-1:1);
ld   = l(end-1:-1:2);
last = first+ld-1;

p = zeros(nlevel,1);
dp = zeros(nlevel,1);

pars=0.0;   % For parseval relation check below

coefs=cell(1,nlevel);

for k = 1:nlevel
    w = first(k):last(k);
    coefs{1,k}=c(w);
    cw = c(w).*c(w);
    pars=pars+sum(cw);  % For parseval relation check below
    p(k) = mean(cw);
    dp(k) = std(cw);
end


scale = (2.^[1:nlevel]');
frequency = scal2frq(scale,wname,dt); 

% Converting to power/energy spectral density
p=p.*(dt*2);

% Calculation of error
dp=dp.*(dt*2);  % turning standard deviation into PSD units
dp=dp./(sqrt(l(end-1:-1:2)));  % turning into standard error of the mean
dp=dp./(p.*log(10));    % propagation of error after taking log operation
dp=dp.*2;   % 95% confidence interval within two standard deviations 
            % assuming Gaussian statistics. 


%     errorbar(log10(frequency),log10(p),dp,'r'); hold all;
    h=plot(log10(frequency),log10(p),'k.-');
%     loglog((frequency),(p),'k.-');
    grid on
    xlabel('log_{10} frequency (Hz)')
    ylabel('log_{10} PSD (signal units)^{2}Hz^{-1}')
    title(['nr of levels = ',int2str(nlevel),'   wavelet is ',wname])
    
%% Parsevals relation check

% Parsevals relation check based on summing the original signal and
% comparing with some of all coefficients at all levels all together. In
% fact this one is probably not needed as it gives no information about the
% correctness of the PSD normalisation. Anyway, we will include it, just to
% show that MATLAB is doing it discrete wavelet transform properly.

answerreal=(sumsqr(y));
answerfreq=(pars + sum(c(1:l(1)).*c(1:l(1))));

parsevalcheck1=answerreal/answerfreq

% Parsevals relation check based on energy contained under the PSD curve
% i.e. the integral of the PSD, which gives the average energy.

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

ip=ip*(mean(c(1:l(1)).*c(1:l(1))))*(dt*2);

ip=ip*(1/dt);

answerfreq=(ip+iu);

parsevalcheck2=answerreal/answerfreq

%% Line fitting

if fitline>0
    
    SlopeRegress(h);
    
end


end