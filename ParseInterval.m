%**************************************************************************
function [signalOut]=ParseInterval(signalIn,startTime,endTime)
%**************************************************************************
% 
% File name: ParseInterval.m
% Author: Khurom H. Kiyani
% File created: 22nd October 2013
% Last updated: 22nd October 2013
% Updated by: Khurom H. Kiyani
% 
% Description:
% ------------
% This m-file takes in a signal which has NaNs and replaces them with a constant value if at the 
% edges of the signal, or replaces them with a linearly interpolated segment if somewhere in the 
% middle of the signal. You basically want to use something like this when a./ you want to do 
% Fourier filtering and b./ it doesn't matter that you have padded NaNs out in this way.
% 
% Input variables/signals:
% 
% 'time':           signal time vector
% 'signal':         signal vector
% 
%**************************************************************************

startTime=datenum(startTime,'dd-mmm-yyyy HH:MM:SS');
endTime=datenum(endTime,'dd-mmm-yyyy HH:MM:SS');

signalOut=signalIn((signalIn(:,1)>startTime),:);
signalOut=signalOut((signalOut(:,1)<endTime),:);

end