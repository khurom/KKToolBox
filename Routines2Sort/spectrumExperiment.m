%
function [F,a]=spectrumExperiment(data,windowsize,ts)

% We will now proceed to construct spectrograms and spectrogram averages
% using short time Fourier transforms (STFT), continuous wavelet transforms
% (CWT) and lastly discrete wavelet transforms (DWT).

% Input: time=time vector; data=data vector; ts=sampling time interval;
% Output:

V=[];

for i=-2:0.1:1.5,
    V=[V 10^(i)];
end

%% STFT spectrogram
% loglog(F(4:end),P((4:end),(50)))
% [S,F,T,P]=spectrogram(data,windowsize,[],[],(1/ts));
[S,F,T,P]=spectrogram(data,windowsize,[],V,1/ts);
% plot(log10(F(4:end)),log10(P((4:end),(50))))
% surf(F,T,10*log10(abs(P')),'EdgeColor','none')
% view(0,90)

% plot(F,a);

% open FastPowerSpecVar.m
% [frequ,Pxx]=FastPowerSpecVar(data(:,11),0.125)

%%
% figure(1):subplot(222);
% plot(time,data);
% 
% figure(1):subplot(224);
% spectrogram(data,3000,[],[],(1/ts),'yaxis');t
% surf(F,T,10*log10(abs(P')),'EdgeColor','none')
% shading interp;

% 2nd argument sets the time resolution; it is the number of points in each
% window, so the higher this is the poorer the time resolution; conversely
% the lower it is the higher the timer resolution
% 5th argumnet is the sampling frequency

% a=(mean(P,2));

% figure(2)
% surf(F,T,10*log10(abs(P)))
           
% win=hann(windowsize);
% P=(S.*conj(S)).*(2/((1/ts)*(win'*win ) )) ;    % This is equivalent to
                                                % what is happening at the
                                                % moment automatically with
                                                % this routine.
            
a=(mean(P,2)); 

% If you are specifying a vector of frequencies then the MATLAB spectrogram routine
% uses a two-sided PSD estimate.
               
figure(1)%:subplot(223);
loglog((F(2:end)),(a(2:end)),'-r');
% plot(log10(F(2:end)),log10(a(2:end)),'*b');
% view(-90,90);t

%% Parsevals relations check

answerreal=(sumsqr(data))/(length(data));
answerfreq=((1/ts)/windowsize)*(sum(a(:)));

parsevalcheck=answerreal/answerfreq

%% pmtm -- psd using multitaper method
% next use pmtm and smoothing; weather movav or regression
% 
% [x1,y1]=pmtmsfr(data(:,11),0.125);
% yy1=smooth(y1,1000);
% xx1=smooth(x1,1000);
% loglog(xx1,yy1)
% hold
% loglog(x1,y1)
% xx1=smooth(x1,1000,'lowess');
% yy1=smooth(y1,1000,'lowess');
% loglog(xx1,yy1)
% yy1=smooth(log10(y1),1000);
% xx1=smooth(log10(x1),1000);
% plot(xx1,yy1)
% yy1=smooth(log10(y1(12:end)),300,'lowess');
% xx1=smooth(log10(x1(12:end)),300,'lowess');
% plot(xx1,yy1)