% %% Old Load data
% 
% clear all; close all;
% 
% % cd ~/khuromOLD/SOLARStuff/NewClusterWork/Data/BurstModeIntervals/20070120_1200
% cd ~/khuromOLD/SOLARStuff/NewClusterWork/Data/BurstModeIntervals/20070130_0000
% % 
% BSC = c_load('BSC4','var');
% BFGM = c_load('B4','var');
% 
% % You should only do this for the 20070130 date
% BSC=BSC(256463:1876571,:);
% BFGM=BFGM(40355:end,:);
%  
% figure(1); hold on;
% spectro(BFGM(:,4),2^14,1/67,'b');
% spectro(BSC(:,4),2^16,1/450,'r');
% % 
% % hold off;
% 
% % 20070120 interval A
% % BSC=BSC(1:2054276,:);
% % BFGM=BFGM(1:307226,:);
% 
% % 20070120 interval B
% % BSC=BSC(2258533:end,:);
% % BFGM=BFGM(337759:end,:);
% 
% cd /Users/khurom/Dropbox/work.plasma/ProjectSurveyDissapationRange/data/Interval_01
% 
% load C4_FGM_20070130_0000_0120.mat;
% load C4_HBR_20070130_0000_0120.mat;
% 
% figure(2); hold on;
% spectro(FGM(:,4),2^14,1/67,'b');
% spectro(STAFF(:,4),2^16,1/450,'r');

%% Put your interval here

clear all; close all;
cd ~/Dropbox/work.plasma/ProjectSurveyDissapationRange/data;

load Interval_ReadFile.mat

q=28;
sc=4;

cd(strcat('Interval_',Interval{q}));

load(strcat('C',num2str(sc),'_FGM_',DateTime{q},'.mat'));
load(strcat('C',num2str(sc),'_HBR_',DateTime{q},'.mat'));

BFGM=FGM;
BSC=STAFF;

clear('FGM','STAFF','sc','q');

% Set parameters

wname='db20';   % wavelet for DWT

% Normal mode 

% leve=4;                     % wavelet frequency cut level
% desiredCadence=1/25;        % final time resolution required
% fourierwindowsize=2^14;     % Size of overlapping windows in PSD

% Burst mode

leve=[8,8,8];               % wavelet frequency cut level for each 
                            % component x,y,z
desiredCadence=1/450;       % final time resolution required
fourierwindowsize=2^16;     % Size of overlapping windows in PSD


%% check for gappy data. This routine will not work well with data gaps due
% to the extensive filtering used.
subplot(2,2,1); plot(BSC(:,1),BSC(:,2)); axis tight;
subplot(2,2,2); plot(BFGM(:,1),BFGM(:,2)); axis tight;
subplot(2,2,3); plot(BSC(1:end-1,1),1./(24*60*60.*(diff(BSC(:,1))))); axis tight;
1/(24*60*60*nanmean(diff(BSC(:,1))))
subplot(2,2,4); plot(BFGM(1:end-1,1),1./(24*60*60.*(diff(BFGM(:,1))))); axis tight;
1/(24*60*60*nanmean(diff(BFGM(:,1))))

%% PSD checker

close all;

for m=2:1:4
    
component=m;

figure(m-1);
% Burstmode
spectro(BFGM(:,component),fourierwindowsize/4,1/67,'b'); hold on;
spectro(BSC(:,component),fourierwindowsize,1/450,'r');

% Normalmode
% spectro(BFGM(:,component),fourierwindowsize,1/22.5,'b'); hold on;
% spectro(BSC(:,component),fourierwindowsize,1/25,'r');

[p,pfreq,~,~,~,~]=wspect(BSC(:,component),desiredCadence,wname);
plot(log10(pfreq(leve(m-1))),log10(p(leve(m-1))),...
    'MarkerFaceColor',[0.0431372560560703 0.517647087574005 0.780392169952393],...
    'MarkerEdgeColor',[0.847058832645416 0.160784319043159 0],...
    'MarkerSize',8,...
    'Marker','o',...
    'LineWidth',2,...
    'LineStyle','none',...
    'Color',[1 1 0]);

clear('component','p','pfreq');

end

%% Close figures

close all;

%% interpolate FGM onto the same timebase as STAFF-SC

BFGMx = interp1(BFGM(:,1),BFGM(:,2),BSC(:,1));
BFGMy = interp1(BFGM(:,1),BFGM(:,3),BSC(:,1));
BFGMz = interp1(BFGM(:,1),BFGM(:,4),BSC(:,1));

BFGM_HF=[BSC(:,1), BFGMx, BFGMy, BFGMz];

clear('BFGMx','BFGMy','BFGMz');

%% Get rid of the NaN's in the timeseries
% These are really at the edges due to the interpolation

a=isnan(BFGM_HF(:,2));
[row,col]=find(a==0);
BFGM_HFnew=BFGM_HF(row,:);
BSC_new=BSC(row,:);

clear('a','row','col','BFGM','BSC','BFGM_HF');

%% low-pass filter FGM

BFGM_LpF=zeros(size(BFGM_HFnew));
BFGM_LpF(:,1)=BFGM_HFnew(:,1);

for i=2:1:4
    
    [dwtwc,l] = wavedec(BFGM_HFnew(:,i),leve(i-1),wname);
    BFGM_LpF(:,i) = wrcoef('a',dwtwc,l,wname,leve(i-1));
    clear('dwtwc','l');
    
end

clear ('i');

%% check spectra

spectro(BFGM_HFnew(:,2),fourierwindowsize,desiredCadence,'-r'); hold on;
spectro(BFGM_LpF(:,2),fourierwindowsize,desiredCadence,'-b');

%% high-pass filter STAFF-SC

BSC_HpF=zeros(size(BSC_new));
BSC_HpF(:,1)=BSC_new(:,1);

for i=2:1:4
    
    [dwtwc,l] = wavedec(BSC_new(:,i),leve(i-1),wname);
    
    for j=1:1:leve(i-1)
        
        BSC_HpF(:,i) =  BSC_HpF(:,i) + wrcoef('d',dwtwc,l,wname,j);
        
    end
    
    clear('dwtwc','l','j');
    
end

clear ('i');

%% check spectra

spectro(BSC_new(:,2),fourierwindowsize,desiredCadence,'-r'); hold on;
spectro(BSC_HpF(:,2),fourierwindowsize,desiredCadence,'-b');

%% glue data together

BTot=zeros(size(BSC_HpF));
BTot(:,1)=BSC_HpF(:,1);

for i=2:1:4
    
    BTot(:,i) = BFGM_LpF(:,i) + BSC_HpF(:,i);
    
end

clear('i');

%% check spectra

for m=2:1:4;

component=m;
    
figure(m-1);
    
spectro(BSC_new(:,component),fourierwindowsize,desiredCadence,'-r'); hold on;
spectro(BFGM_HFnew(:,component),fourierwindowsize,desiredCadence,'-b');
spectro(BTot(:,component),fourierwindowsize,desiredCadence,'-k');

[p,pfreq,~,~,~,~]=wspect(BSC_new(:,component),desiredCadence,wname);
plot(log10(pfreq(leve(m-1))),log10(p(leve(m-1))),...
    'MarkerFaceColor',[0.0431372560560703 0.517647087574005 0.780392169952393],...
    'MarkerEdgeColor',[0.847058832645416 0.160784319043159 0],...
    'MarkerSize',8,...
    'Marker','o',...
    'LineWidth',2,...
    'LineStyle','none',...
    'Color',[1 1 0]);

end

clear component;

%% save and then delete crap

save('C4_BTot_20070328_0400_0440.mat','BTot');

clear all;
