
clear all; close all;

cd ~/Dropbox/ProjectSurveyDissapationRange/data;

load('Interval_ReadFile.mat','Interval','DateTime');

q=03;
sc=4;

cd(strcat('Interval_',Interval{q}));

load(strcat('C',num2str(sc),'_FGM_',DateTime{q},'.mat'));
load(strcat('C',num2str(sc),'_NBR_',DateTime{q},'.mat'));

BFGM=FGM;
BSC=STAFF;

clear('FGM','STAFF');

%%

% Set parameters

wname='db20';   % wavelet for DWT

% Normal mode 

leve=[4,4,4];               % wavelet frequency cut level for each component x,y,z
desiredCadence=1/25;        % final time resolution required
fourierwindowsize=2^12;     % Size of overlapping windows in PSD

% Burst mode
% 
% leve=[7,7,7];               % wavelet frequency cut level for each component x,y,z
% desiredCadence=1/450;       % final time resolution required
% fourierwindowsize=2^13;     % Size of overlapping windows in PSD

% PSD checker

close all;

for m=2:1:4
    
component=m;

figure(m-1);
% Burstmode
% spectro(BFGM(:,component),fourierwindowsize/4,1/67,'b'); hold on;
% spectro(BSC(:,component),fourierwindowsize,1/450,'r');

% Normalmode
spectro(BFGM(:,component),fourierwindowsize,1/22.5,'b'); hold on;
spectro(BSC(:,component),fourierwindowsize,1/25,'r');

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

%% Close figures

close all;

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

ind=2;

spectro(BSC_new(:,ind),fourierwindowsize,desiredCadence,'-r'); hold on;
spectro(BSC_HpF(:,ind),fourierwindowsize,desiredCadence,'-b');

clear ind;

%% Close figures

close all;

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

save(strcat('C',num2str(sc),'_BTot_',DateTime{q},'.mat'),'BTot');

