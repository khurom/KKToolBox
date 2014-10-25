%%
% m-file to calculate the parallel and perpendicular plasma beta as well as
% proton temperature anisotropies.
% Author: Khurom Kiyani
% Date: 31/03/09
% note: to do ion calculations change Cp4 to Ch1
%%

% Note: you have to truncate the electron data time-series if you want to
% have it on the same time interval as the magnetic field for the below
% calculations.

ts_Tr=timeseries((tperp./tpar),time);

ts1 = timeseries(B4(:,5),B4(:,1));
ts2 = timeseries(edens,time);
% [ts1 ts2] = synchronize(ts1,ts2,'Intersection','tolerance',1e-11);

ts_Tpar=timeseries(tpar,time);
ts_Bmag=resample(ts1,time,'linear');
ts_n=resample(ts2,time,'zoh');

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