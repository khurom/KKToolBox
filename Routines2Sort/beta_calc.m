%%
% m-file to calculate the plasma beta
% Author: Khurom Kiyani
% Date: 31/03/09
% note: to do ion calculations change Cp4 to Ch1
%%

ts1 = timeseries(Bmag(:,2),Bmag(:,1));
ts2 = timeseries(NCh1(:,2),NCh1(:,1));
% [ts1 ts2] = synchronize(ts1,ts2,'Intersection','tolerance',1e-11);

Tmag=sqrt( (TperpCh1(:,2).*TperpCh1(:,2)) + (TparCh1(:,2).*TparCh1(:,2)) );

ts_Tmag=timeseries(Tmag(:),TparCh1(:,1));
ts_Bmag=resample(ts1,TparCh1(:,1),'linear');
ts_n=resample(ts2,TparCh1(:,1),'zoh');

% This const is a constant which combines all the necessary unit
% conversions with the Boltzman constant and 8*pi
const=8*pi*1.3807;

beta=(const.*ts_n.data.*ts_Tmag.data)./(ts_Bmag.data.*ts_Bmag.data);
ts_beta=timeseries(beta,ts_Tmag.time);

clear('const','ts1','ts2','ts_Bmag','ts_n')

figure(1);
plot(ts_beta,'.r');
