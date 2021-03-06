function [F,a,era]=spectro(data,windowsize,ts,plotoption,fitline)
% We will now proceed to construct spectrograms and spectrogram averages
% using a short time Fourier transform (STFT).
% Input: time=time vector; data=data vector; ts=sampling time interval;
% Output:
% Setting up font sizes for plots
set(0,'defaultaxesfontsize',16);
set(0,'defaulttextfontsize',16);

if nargin<4, plotoption = 'b'; end
if nargin<5, fitline = 0; end

%% STFT spectrogram

% %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% % Comment out this section if you are not providing frequencies. This uses
% % the Goertzel algorithm and the PSD becomes a two sided one (hence
% % factor of 2 with P. note:Parsevals relation doesn't work very well here.
% %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% 
% V=[];
% 
% for i=-2:0.1:2,
%     V=[V 10^(i)];
% end
% 
% [S,F,T,P]=spectrogram(data,windowsize,[],V,1/ts);
% P=P.*2;
% %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Comment out this section if you ARE providing a frequency vector.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[S,F,T,P]=spectrogram(data,windowsize,[],[],(1/ts));

          
% win=hann(windowsize);
% P=(S.*conj(S)).*(2/((1/ts)*(win'*win ) )) ;    % This is equivalent to
                                                % what is happening at the
                                                % moment automatically with
                                                % this routine.
                                                
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            
a=(nanmean(P,2));

% disp('Number of windows being averaged: ' );
% disp(size(P,2));

fprintf('%s %d ','Number of windows being averaged: ', size(P,2));

% Calculate standard error of the mean assuming the statistics are Gaussian
era=(nanstd(P,1,2));
era=(2.*log10(exp(1))).*(era./((sqrt(numel(a))).*a));
% Factor 2 gives 95% confidence interval

% subplot(211);
% errorbar(log10(F(2:end)),log10(a(2:end)),era(2:end),'-k'); hold on;
h=loglog(F(2:end),a(2:end),plotoption);
% loglog((F(2:end)),(a(2:end)),plotoption);
% plot((F(2:end)),(a(2:end)),plotoption);

ylabel('PSD ((signal units)^{2}Hz^{-1})','FontSize',16);

% subplot(212);
% errorbar(log10(F(2:end)),zeros((numel(a)-1),1),era(2:end));

xlabel('Frequency (Hz)','FontSize',16);
% title('95% confidence interval error ((signal units)^{2} Hz^{-1})','FontWeight','bold',...
%     'FontSize',14);

%% Parsevals relation check

answerreal=(nansum(data.*data))/(length(data));
answerfreq=((1/ts)/windowsize)*(nansum(a(:)));

parsevalcheck=answerreal/answerfreq

%% Line fitting

if fitline>0
    
    SlopeRegress(h);
    
end

end