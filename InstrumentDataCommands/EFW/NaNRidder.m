% load mEDSIf.mat;

clear all; close all;

% cd ~/khuromOLD/SOLARStuff/NewClusterWork/Data/BurstModeIntervals/20070120_1200
cd ~/khuromOLD/SOLARStuff/NewClusterWork/Data/BurstModeIntervals/20070130_0000
% cd ~/khuromOLD/SOLARStuff/NewClusterWork/Data/BurstModeIntervals/noise/20070630_1420

c_load diE4p1234
%% Bit which gets rid of the NaNs
a=isnan(diE4p1234);
[row,col]=find(a(:,2:3)==0);
b=unique(row);
Efield=diE4p1234(b,:);

comp=3; % x component=2, y component=3, z component=4

% irf_plot a4
% irf_psd(a4(:,1:3),16384*4);hold on
% FastPowerSpecVar(a4(:,3),1/450); hold on;
% spectro(Efield(:,comp),16384*4,1/450);
% y=diff(a4(:,1));
% plot(y);

number_of_NaNs_removed_is=numel(diE4p1234(:,1))-numel(Efield(:,1))

clear('a','b','row','col','ans','diE4p1234','number_of_NaNs_removed_is');

%%

ts=1/450;

time=Efield(:,1);

%% Bit to handle gaps                            
                            
time_size=numel(time);

gap=10*ts;

a=diff(time);

b=find(a>gap); % array indices of TIME where a gap BEGINS

b=[0;b;time_size];  % pad on front and end with (beginning-1)=0 and end 
                    % indices of 'time' array.
                    
% Ignore the first one and the last one as they consist of small samples

clear('time');

%%

time=[];
timesum=0.0;

intervals=cell(numel(b)-1,1);

for i=1:1:(numel(b)-1)
    
    intervals{i}=[Efield(b(i)+1:b(i+1),1) Efield(b(i)+1:b(i+1),comp)];
    timesum=timesum + ( Efield(b(i+1),1) - Efield(b(i)+1,1) );
    time=[time timesum];
    
end

clear('i');

%% PSD calculation
PSD=[];

for i=2:1:numel(intervals)-1
    
    [f,a]=spectro(intervals{i}(:,2),32768,1/450,'-b'); hold on;
    PSD=[PSD a];
    
end

a=mean(PSD,2);

%%

c_load diBSC4;
c_load diB4; 
%% 
figure(1);
spectro(diB4(:,4),16384*2,1/67,'.b'); hold on;

spectro(diBSC4(:,4),16384*2,1/450,'.r'); hold on;

% spectro(Efield(:,comp),16384*4,1/450,'.r'); hold on; % with gaps

plot(log10(f(2:end)),log10(a(2:end)),'.k');     % without gaps
% loglog((f(2:end)),(a(2:end)),'.r');



