%**************************************************************************
function [out]=ConstNan(time,signal)
%**************************************************************************
% 
% File name: ConstNan.m
% Author: Khurom H. Kiyani
% File created: 21st October 2013
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

% linearly interpolating all the middle NaNs
timetemp=time(~isnan(signal));
signal(isnan(signal))=[];
out=interp1(timetemp,signal,time);

clear timetemp;

% interpolating the edges with constants
ind=find(~isnan(out));
out(1:ind(1)-1)=out(ind(1));
out(ind(end)+1:end)=out(ind(end));

clear ind;

end