%**************************************************************************
function SlopeRegress(plothandle)
%**************************************************************************
% 
% File name: SlopeRegress.m
% Author: Khurom H. Kiyani
% File created: 14th March 2014
% Last updated: 16th March 2014
% Updated by: Khurom H. Kiyani
% 
% Description:
% ------------
% This function plots a linear fit on a curve that you select from an 
% existing line in a figure. The linear fit is done using ordinary least 
% squares. The errors provided are the average of the upper and lower 95% 
% quantiles. Make sure that a line/curve is selected already from the 
% figure on the screen OR that a handle for such a plot for a line is 
% passed as an argument.
% 
% Input variable:
% ---------------
% 'plothandle': handle to a plot for a line/curve.
% 
% Notes:
% ------
% 1./ Need to recognise and transform accordingly log coordinate axes.
% 2./ Need to use your own regression algorith -- it's fairly easy. This
% current one is a more complciated linear regression and not a OLS
% 
%**************************************************************************

if nargin<1,
    
    x=get(gco,'XData'); x=x';
    y=get(gco,'YData'); y=y';

else
    
    x=get(plothandle,'XData'); x=x';
    y=get(plothandle,'YData'); y=y';
    
end

shift=0.0; %0.3

[xin,~]=ginput(2);

idx=x>xin(1) & x<xin(2);

xnew=[ones(size(x(idx))) x(idx)];
[b,berr]=regress(y(idx),xnew);

x1=xin(1):abs(xin(2)-xin(1))/10:xin(2);
y1=b(2).*x1 + b(1);
hold on
plot(x1,y1+shift,'r')

% intercept=b(1);
slope=b(2);
% errintercept=(abs(berr(1,1)-b(1))+abs(berr(1,2)-b(1)))/2;
errslope=(abs(berr(2,1)-b(2))+abs(berr(2,2)-b(2)))/2;

disp(['Slope=',num2str(slope),' ','+/-',num2str(errslope)]);

annotation('textbox',...
    [0.44 0.71 0.22 0.13],'String',{['Slope= ',num2str(slope),' ',...
    '\pm',num2str(errslope)]},...
    'FitBoxToText','off',...
    'LineStyle','none');

end

