%**************************************************************************
% 
% Function:     OutlierThresh.m
% 
% Description:  This function implements a thresholding or conditioning
%               scheme, such that a percentage of extreme values or
%               outliers of the input data set are excluded.
% 
% Inputs:       'ThreshIN' is the data (column vector) to be thresholded.
%               'perc' is the percentage value (real number) of points to
%               be excluded.
% 
% Outputs:      'ThreshOUT' is the data (column vector) after 'perc'
%               percentage of points have been excluded.
% 
% Author:       Khurom Kiyani
% Date:         6th September 2011
% Last updated: 28th February 2013
% 
% Copyright 2011 Imperial College London. All rights reserved.
% 
%**************************************************************************


function [ThreshOUT]=OutlierThresh(ThreshIN,perc)

perc=perc/100; 

if iscell(ThreshIN)>0

for m=1:1:length(ThreshIN)
    cut_off=round(perc*(length(RemoveNaN(ThreshIN{m}))));  
    [a,b]=sort(abs(ThreshIN{m}));
    g=find(isnan(a));
    ThreshIN{m}(b(g(1)-cut_off:g(1)-1))=NaN;
    clear('cut_off','a','b','g');
end

clear('m');

else

    cut_off=round(perc*(length(RemoveNaN(ThreshIN))));  
    [a,b]=sort(abs(ThreshIN));
    g=find(isnan(a));
    ThreshIN(b(g(1)-cut_off:g(1)-1))=NaN;
    clear('cut_off','a','b','g');
    
end
    
ThreshOUT=ThreshIN;


end