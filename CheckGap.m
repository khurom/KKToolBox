
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

% check for gappy data. This routine will not work well with data gaps due
% to the extensive filtering used.
subplot(2,2,1); plot(BSC(:,1),BSC(:,2)); axis tight;
subplot(2,2,2); plot(BFGM(:,1),BFGM(:,2)); axis tight;
subplot(2,2,3); plot(BSC(1:end-1,1),1./(24*60*60.*(diff(BSC(:,1))))); axis tight;
1/(24*60*60*nanmean(diff(BSC(:,1))))
subplot(2,2,4); plot(BFGM(1:end-1,1),1./(24*60*60.*(diff(BFGM(:,1))))); axis tight;
1/(24*60*60*nanmean(diff(BFGM(:,1))))