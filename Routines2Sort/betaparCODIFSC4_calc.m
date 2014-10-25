%%
% m-file to calculate the parallel and perpendicular plasma beta as well as
% proton temperature anisotropies.
% Author: Khurom Kiyani
% Date: 31/03/09
% note: to do ion calculations change Cp4 to Ch1
%%

%% Load particle data

load mCISR.mat;

%%

ts_Tr=timeseries((TperpCp4(:,2)./TparCp4(:,2)),TperpCp4(:,1));

ts1 = timeseries(B4(:,5),B4(:,1));
ts2 = timeseries(NCp4(:,2),NCp4(:,1));
% [ts1 ts2] = synchronize(ts1,ts2,'Intersection','tolerance',1e-11);

ts_Tpar=timeseries(TparCp4(:,2),TparCp4(:,1));
ts_Bmag=resample(ts1,TparCp4(:,1),'linear');
ts_n=resample(ts2,TparCp4(:,1),'zoh');

% This const is a constant which combines all the necessary unit
% conversions with the Boltzman constant and 8*pi
const=8*pi*1.3807;

beta_par=(const.*ts_n.data.*ts_Tpar.data)./(ts_Bmag.data.*ts_Bmag.data);
ts_betapar=timeseries(beta_par,ts_Tpar.time);

clear('const','beta_par','ts1','ts2','ts_Tpar','ts_Bmag','ts_n')

figure(1);
plot(ts_Tr,'.b'); hold all; plot(ts_betapar,'.r');

figure(2);
loglog(ts_betapar.data,ts_Tr.data,'.b')