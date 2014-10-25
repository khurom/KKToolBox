%%  MOVING AVERAGE
%   This function implements a simple moving average using a low-pass
%   eliptic filter.
%   Input: 'time' variable (array column)
%          'variable' to average (array column)
%          'parts' in variable to average over
%   Output:'varavr' averaged variable column
%   Author: Khurom Kiyani
%   Date: 18 October 2013
%**************************************************************************

function [varavr]=movav(time,variable,parts,samplefreq)

varavr=irf_filt(variable,0,1/parts,samplefreq,3);

plot(time,variable); hold on;
plot(time,varavr,'r','Linewidth',2.0);

end