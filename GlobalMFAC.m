%**************************************************************************
function [MeanField,Fluctuations] = GlobalMFAC(Data,LPFiltFreq,Fs,plotopt)
%**************************************************************************
% 
% File name: GlobalMFAC.m
% Author: Khurom H. Kiyani
% File created: 10th September 2013
% Last updated: 18th March 2014
% 
% Description:
% ------------
% Global Mean Field Alligned Coordinates (GlobalMFAC).
% This m-file essentially performs a low-pass filter at a chosen frequency
% to generate the background field. The low-pass filter used is an Elliptic
% filter. The remaining field is the Fluctuations. 
% 
% Input variables:
% 
% 'Data':           Field *vector* time-series. Assumes that the first
%                   column is the time vector. 
% 'LPFiltFreq':     Low-pass filter cut-off frequency; if this is exactly
%                   zero the low pass filter will not be used and a mean
%                   average will be performed over the entire interval.
% 'Fs':             Sampling frequency of the signal.
% 'plotopt':        If this is '0' no plot will appear. If '1' (actually 
%                   just greater than '0', a plot will appear. 
% 
% Output variables:
% 
% 'MeanField':      Resultant mean background field. 
% 'Fluctuations':   Resultant fluctuations.
% 
%**************************************************************************
% 
% Additional Notes:
% -----------------
% 1./ Need to make a better criteria for when it goes to a linear trend.
%     LPFiltFreq<0.0005 will not do.
% 2./ Use better filters? Or will elliptic do? [Consider changing this to a
%     Chebychev Type II].
% 
%**************************************************************************

if nargin<4, plotopt = 1; end

cols=size(data); cols=cols(2);

% Do filtering

if LPFiltFreq==0    % i.e. Global average field
    
    MeanField = repmat(mean(Data(:,2:end)),length(Data),1);
    MeanField = [Data(:,1) , MeanField];
    
elseif LPFiltFreq<0.0005
    
    MeanField = [Data(:,1) , Data(:,2:end) - detrend(Data(:,2:end))];
    
else                % the low-pass filter. Order 3 elliptic filter, assumes 
                    % first column is time. 
                    
    MeanField = irf_filt(Data,0,LPFiltFreq,Fs,3);
    
end

% Construct fluctuations

Fluctuations = [Data(:,1) , Data(:,2:end) - MeanField(:,2:end)];

if plotopt>0
    
    for m=2:cols
        
        subplot(3,2,m+(m-3));
        plot(Data(:,1),Data(:,m)); hold on; plot(MeanField(:,1),...
            MeanField(:,m),'r','LineWidth',3);
        axis tight;
        subplot(3,2,m+(m-3)+1);
        plot(Fluctuations(:,1),Fluctuations(:,m));
        axis tight;
        
    end
    
end

end

