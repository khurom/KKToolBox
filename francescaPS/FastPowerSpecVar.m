function [frequ,Pxx,X,phase]=FastPowerSpecVar(data,ts)

% ts; N(number of points for fft);

% P=load(datafile);

% data=P(:,columNum); clear P;

fs = 1./ts;   
fN = 0.5*fs;        % Nyquist frequency
%
% N = 16384*4;   
N = 2^13;
N2 = nextpow2(N);   nfft = pow2(N2);
fx = eval_freq(nfft);       %   build the frequency vector, f1=[0,1] 
                            %   normalized to the Nyquist frequency
x = shape_inp(data,N);    % rearrange the data into an array of size (N,M)
X = eval_X(x);          % calculate the Fourier transform of the array of samples 
X = conj(X);        clear x

[nk,M] = size(X);
%   compute time vector over windows
tj=[];
for j=1:M
    tj = cat(1,tj,j*N/2./fs);
end
%
% Pxx=mean(X.*conj(X),2)./nfft; % The normalisation f.m.poli was using
% before; it does not take into account the windowing power, factor of 2 
% from it being a real signal and more
% importantly the sampling frequency. Note in the correct normalisation
% below, that there is only fs in the denominator; nfft has cancelled with
% the nfft involved in calculating the window power in eval_X (see below).
Pxx=( 2*mean(X.*conj(X),2) )./( fs );  

% Phase stuff
% ----------

phase = unwrap(angle(mean(X,2)));
phase=phase(2:end);

frequ=fx*fN;

Pxx=((Pxx(2:end)));

frequ=((frequ(2:end)));

size(Pxx)
size(frequ)

% % Create axes
% axes('FontSize',16);
% box('on');
% grid('on');
% hold('all');
% 
figure(1);
% loglog((frequ),(Pxx));
plot(log10(frequ),log10(Pxx),'-r');

xlabel('log_{10} Frequency (Hz)','FontSize',16);

ylabel('log_{10} PSD ((signal units)^{2}Hz^{-1})','FontSize',16);

figure(2);

plot(log10(frequ),(phase*180/pi),'-b');

% -------------------------------------------------------------------------
function f = eval_freq(nfft)

%______________________________________________________________________
%
%   function f=eval_freq(nfft) calculate the frequency vector given the
%   number of FFT points, nfft. The final frequency is in the range [0,1], 
%   where f=1 is the Nyquist frequency.
%
%   Francesca M. Poli
%______________________________________________________________________

f = [0:nfft/2]/(nfft/2); % normalized frequency


% -------------------------------------------------------------------------
function x = shape_inp(data,N)

%________________________________________________________________
%
%	function shape_inp(data,N)	is called by function setinp.m
%
%	reshapes arrays of size (nt,1) into array of size (N,M)
%	using overlapping of 50% between time windows	
%
%	Francesca M. Poli
%________________________________________________________________

nt = length(data);

m = fix(2*nt./N);

j = [1:m-1];
k = [0:N-1]';
jm = (j-1)*N/2+1;
ind = repmat(jm,N,1) + repmat(k,1,m-1);

x = data(ind);


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


[N,m] = size(x);
N0 = nextpow2(N);
nfft = pow2(N0);

F0 = hanning(N);
yn = zeros(nfft,m); % ready for zero padding

xm = detrend(x,'constant');
xl = detrend(xm,'linear');
xlw = xl.*repmat(F0,1,m);
yn(1:N,:)=xlw;       % zero padding

fft_y1 = fft(yn,nfft,1);     % FFT of signal #1

% X = fft_y1(1:nfft./2+1,:);  % this is how it was before
X = fft_y1(1:nfft./2+1,:)./sqrt(F0'*F0);    % divide by window power
%

% -------------------------------------------------------------------------
% wname='morl';
% A0 = [1 2 4 8 16 32:16:128 256 512 1024 2048 4096];
% F0 = scal2frq(A0,wname,ts);
% j0 = [F0 < fN];
% sj = A0(j0);   % scales
% fw = F0(j0);
% 
% Xw = cwt(dataSW,sj,wname);

% [EOF]

