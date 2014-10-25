clear all; close all;

% cd ~/khuromOLD/SOLARStuff/NewClusterWork/Data/BurstModeIntervals/20070120_1200
cd ~/khuromOLD/SOLARStuff/NewClusterWork/Data/BurstModeIntervals/20070130_0000

BFGM = c_load('B4','var');

BFGM=BFGM(40355:end,:);

% BSC = c_load('BSC4','var');
% 
% BFGM=BSC(256463:1876571,:);

[S,F,T,P]=spectrogram(BFGM(:,4),4096,[],[],67);

% [S,F,T,P]=spectrogram(BFGM(:,4),4096*8,[],[],450);

imagesc(T,log10(F(5:100)),log10(P(5:100,:))); % ,[-6 1])
%%
for i=1:117, P(:,i)=P(:,i).*(F(:).^(5/3)); end
%%
imagesc(T,log10(F(5:50)),log10(P(5:50,:)))

%%
figure(2);
a=mean(P,2);

loglog(F(3:end),a(3:end))
%%
contourf(T,log10(F(2:400)),log10(P(2:400,:)))

%%

XI=[15:5:3591];

YI=[0:0.01:33.5];
YI=YI';

ZI = interp2(T,F,P,XI,YI,'linear');

imagesc(XI,log10(YI(2:200)),log10(ZI(2:200,:)))

%%
for i=1:716, ZI(:,i)=ZI(:,i).*(YI(:).^(5/3)); end
%%
imagesc(XI,log10(YI(5:200)),log10(ZI(5:200,:)))


