function SlopeRegress(arg1,arg2,arg3)
%
% SLOPEREGRESS Simple and flexible linear regression to estimate slope and
% provide associate errors
% 
% File name: SlopeRegress.m
% Author: Khurom H. Kiyani
% File created: 14th March 2014
% Institute: University of Warwick
% E-mail: k.kiyani@warwick.ac.uk
% 
% Last updated: 25th May 2015
% Updated by: Khurom H. Kiyani
% 
% Description:
% ------------
% This function plots a linear fit on a curve that you select from an 
% existing line in a figure. The linear fit is done using ordinary least 
% squares. The errors provided are the average of the upper and lower 95% 
% quantiles and assume that the underlying distribution of residuals is a 
% Gaussian with zero mean and constant variance. Make sure that a line is 
% selected already from the figure on the screen OR that a handle for such 
% a plot for a line is passed as an argument.
% 
% Input variable:
% ---------------
% If only one argument is provided: 
%   'arg1': handle to a plot for a line/curve.
% 
% If two arguments are provided:
%   'arg1': beginning x-coordinate of the fit.
%   'arg2': end x-coordinate of the fit.
% 
% If three arguments are provided:
%   'arg1': handle to a plot for a line/curve.
%   'arg2': beginning x-coordinate of the fit.
%   'arg3': end x-coordinate of the fit.
% 
% Notes:
% ------
% Needs the Statistics toolbox for the Quantile function 'tinv'. This is
% needed to calculate the t-score at which one has 95% confidence intervals
% for the Student's t-distribution; this is the underlying distribution for
% the slope and intercept parameters calculated. 'tinv' is a complicated 
% function which is computed with different sets of approximations.  
% 
% See also TINV

% -------------------------------------------------------------------------
% ARGUMENT CHECKING
% -------------------------------------------------------------------------

if nargin<2,
    
    if nargin<1,
    
        arg1=gco;
    
    elseif ~ishandle(arg1)
        
        error('myApp:argChk', 'Wrong type of input argument');
        
    end
    
    x=get(arg1,'XData'); x=x';
    y=get(arg1,'YData'); y=y';
    
    clear arg1
    
    [xin,~]=ginput(2);  % Take input from the GUI cross-hairs
    
    if xin(1)>xin(2)
       
        error('myApp:argChk', '1st x-coord has to be less than 2nd x-coord');
        
    end
    
    idx=x>xin(1) & x<xin(2);
    
elseif nargin<3

    if arg1>arg2
       
        error('myApp:argChk', '1st x-coord has to be less than 2nd x-coord');
        
    end

    x=get(gco,'XData'); x=x';
    y=get(gco,'YData'); y=y';
    
    xin=[arg1,arg2];
    
    clear arg1 arg2
    
    idx=x>xin(1) & x<xin(2);
    
else
    
    if ~ishandle(arg1)
       
        error('myApp:argChk', 'Wrong type of input argument');
        
    end
    
    if arg2>arg3
       
        error('myApp:argChk', '1st x-coord has to be less than 2nd x-coord');
        
    end
    
    x=get(arg1,'XData'); x=x';
    y=get(arg1,'YData'); y=y';
    
    xin=[arg2,arg3];
    
    clear arg1 arg2 arg3
    
    idx=x>xin(1) & x<xin(2);
    
end

    % Find if the axes are linear or logarithmic
    xAxisScale=get(gca,'Xscale');
    yAxisScale=get(gca,'Yscale');
    

% -------------------------------------------------------------------------
% FIT LINE PLOT ASTHETICS
% -------------------------------------------------------------------------

shift=0.0;  % If you want to move the fit line above the curve that you are 
            % fitting change this to something like 0.3

% -------------------------------------------------------------------------
% THE REGRESSION
% -------------------------------------------------------------------------

if (strcmp(xAxisScale,'log'))>0 && (strcmp(yAxisScale,'log'))>0

    xnew=[ones(size(x(idx))) log10(x(idx))];
    [b,bstderr]=lscov(xnew,log10(y(idx)));

    x1=logspace(log10(xin(1)),log10(xin(2)),10);
    y1=(10^(b(1) + shift)).*(x1.^b(2));
    hold on
    loglog(x1,y1,'r')

elseif (strcmp(xAxisScale,'linear'))>0 && (strcmp(yAxisScale,'linear'))>0

    xnew=[ones(size(x(idx))) x(idx)];
    [b,bstderr]=lscov(xnew,y(idx));
    
    x1=linspace(xin(1),xin(2),10);
    y1=b(2).*x1 + b(1) + shift;
    hold on
    plot(x1,y1,'r');
    
elseif (strcmp(xAxisScale,'linear'))>0 && (strcmp(yAxisScale,'log'))>0
    
    xnew=[ones(size(x(idx))) x(idx)];
    [b,bstderr]=lscov(xnew,log10(y(idx)));
    
    x1=linspace(xin(1),xin(2),10);
    y1=10.^(x1.*b(2)+b(1)+10.^(shift));
    hold on
    semilogy(x1,y1,'r');
    
elseif (strcmp(xAxisScale,'log'))>0 && (strcmp(yAxisScale,'linear'))>0
    
    xnew=[ones(size(x(idx))) log10(x(idx))];
    [b,bstderr]=lscov(xnew,y(idx));
    
    x1=logspace(log10(xin(1)),log10(xin(2)),10);
    y1=b(2).*log10(x1) + b(1) + shift;
    hold on
    semilogx(x1,y1,'r');
    
end
% -------------------------------------------------------------------------

confInt=95;     % confidence interval in percent

alpha=1-(confInt/100); clear confInt  

% for two-sided Student's t-distribution with 'length(y)-2' degrees of
% freedom

tval=tinv((1-(alpha/2)),length(y)-2);

intercept=b(1);
errintercept=tval*bstderr(1);

slope=b(2);
errslope=tval*bstderr(2);

% Display result on the command
disp(['Slope=',num2str(slope),' ','+/-',num2str(errslope)]);
disp(['Intercept=',num2str(intercept),' ','+/-',num2str(errintercept)]);

% Display result on the figure
annotation('textbox',...
    [0.44 0.71 0.22 0.13],'String',{['Slope= ',num2str(slope),' ',...
    '\pm',num2str(errslope)]},...
    'FitBoxToText','off',...
    'LineStyle','none');

end

