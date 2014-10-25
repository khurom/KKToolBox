function [frequ,Pxx,X]=DataCubePSD(datacube,dx)

% ts; N(number of points for fft);

% datacube: grid of data
% dx: grid spacing
% N: Window size

fs = 1/dx;   
fN = 0.5*fs;        % Nyquist frequency

[N,~,~] = size(datacube);

nfft = pow2(nextpow2(N)); clear N;

frequ = (0:nfft/2)/(nfft/2); % build the normalized frequency vector, f1=[0,1]
frequ=frequ*fN; % normalized to the Nyquist frequency

% DON'T THINK WE NEED THIS
% x = shape_inp(datacube,N); % rearrange the data into an array of size (N,M)

X = eval_X(datacube); % calculate the Fourier transform of the array of samples 
X = conj(X);

clear datacube

size(X)

%
% Pxx=mean(X.*conj(X),2)./nfft; % The normalisation f.m.poli was using
% before; it does not take into account the windowing power, factor of 2 
% from it being a real signal and more
% importantly the sampling frequency. Note in the correct normalisation
% below, that there is only fs in the denominator; nfft has cancelled with
% the nfft involved in calculating the window power in eval_X (see below).

Pxx=( 2*mean(X.*conj(X),2) )./( fs );  

Pxx=((Pxx(2:end)));     % We do this because we will plot this on a loglog 
frequ=((frequ(2:end))); % plot and log of zero frequency will become a NaN

size(Pxx)
size(frequ)

% Plotting
figure(1);
plot(log10(frequ),log10(Pxx),'-b');
xlabel('log_{10} Frequency (Hz)','FontSize',16);
ylabel('log_{10} PSD ((signal units)^{2}Hz^{-1})','FontSize',16);

% % Phase stuff
% phase = unwrap(angle(mean(X,2)));
% phase=phase(2:end);
% 
% figure(2);
% plot(log10(frequ),(phase*180/pi),'-b');

end

% -------------------------------------------------------------------------
function X = eval_X(x)

%_________________________________________________________
%
%   function X = eval_X(x)
%   
%   calculate the FFT of array x, of size (N,M)
%   and gives the output X with size (nfft,M).
%   M is the number of subsamples extracted from a time sequence of length
%   nt. It is assumed that you have used the MATLAB routine setinp.m to
%   arrange your data in an array of size (N,M) before calling eval_X.m
%
%   A hanning window is applied to the input data before calculating the
%   Fourier Transform via the FFT algorithm. 
%   The number of points for the FFT transform (nfft) is taken as the 
%   nearest power of 2 of
%   N. Zeros are added to the samples in the case N is lower than the
%   number of nfft points.
%   Before transforming, the average value is subtracted from each
%   subsample. Similarly, any linear trend in the time evolution is
%   eliminated as well.
%
%   Francesca M. Poli
%_________________________________________________________

% The zero padding here might actually not be needed as the FFT routine in
% Matlab does this as standard.

[N,m,p] = size(x);
nfft = pow2(nextpow2(N));

F0 = tukeywin(N);      % Windowing function
% F0 = hanning(N);
% F0 = ones(N,1);
yn = zeros(nfft,m*p); % ready for zero padding

x=reshape(x,N,m*p);

xm = detrend(x,'constant');
xl = detrend(xm,'linear');

xlw = xl.*repmat(F0,[1,m*p]); % Multiple by the window defined earlier
yn(1:N,:,:)=xlw;       % zero padding

fft_y1 = fft(yn,nfft,1);     % FFT of signal along the first dimension

% X = fft_y1(1:nfft./2+1,:);  % this is how it was before
X = fft_y1(1:nfft./2+1,:)./sqrt(F0'*F0);    
% divide by square root of the window power would be equivalent to dividing
% by sqrt(F0'*F0)/sqrt(nfft), but we leave this nfft out as it will cancel
% later on above anyway. Also we use nfft here instead of N because we 
% have zero padded.  
%

% X=reshape(X,nfft,m*p);

end

