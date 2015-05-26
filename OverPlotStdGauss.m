function OverPlotStdGauss
%
% SLOPEREGRESS Plots a standardised Gaussian PDF
% File name: OverPlotStdGauss.m
% Author: Khurom H. Kiyani
% File created: 22nd March 2015
% Institute: University of Warwick
% E-mail: k.kiyani@warwick.ac.uk
% 
% Last updated: 22nd May 2015
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
% requires Statistics Toolbox

x = [-10:.1:10];
norm = normpdf(x,0,1);
semilogy(x,norm,'--')

clear x norm