function [frequ,Pxx]=FastPowerSpec(datafile,columNum,ts)

% ts; N(number of points for fft);

P=load(datafile);

data=P(:,columNum); clear P;

fs = 1./ts;   
fN = 0.5*fs;        % Nyquist frequency
%
N = 4000;    N2 = nextpow2(N);   nfft = pow2(N2);
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
Pxx=mean(X.*conj(X),2)./nfft;

frequ=fx*fN;

% % Create axes
% axes('FontSize',16);
% box('on');
% grid('on');
% hold('all');
% 
% figure(1);
% plot(log10(frequ(2:end)),log10(Pxx(2:end)));
% 
% xlabel('frequency (Hz)','FontSize',16);
% 
% ylabel('Power','FontSize',16);





% wname='morl';
% A0 = [1 2 4 8 16 32:16:128 256 512 1024 2048 4096];
% F0 = scal2frq(A0,wname,ts);
% j0 = [F0 < fN];
% sj = A0(j0);   % scales
% fw = F0(j0);
% 
% Xw = cwt(dataSW,sj,wname);



