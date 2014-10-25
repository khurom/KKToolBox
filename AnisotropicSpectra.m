%**************************************************************************
function [p,dp]=AnisotropicSpectra(coefs,nlevel,lengthofsamples,dt)
%**************************************************************************
% 
% File name: AnisotropicSpectra.m
% Author: Khurom H. Kiyani
% File created: 2nd December 2013
% Last updated: 3rd December 2013
% Updated by: Khurom H. Kiyani
% 
% Description:
% ------------
% Given the associate real wavelet coefficients, this function calculates 
% the Power Spectral Density along with its errors. It assumes that the
% errors on the mean are Gaussian distributed. The scales are assumed to be
% calculated elsewhere.
% 
% Input variables:
% 
% 'coefs':              A cell array of the real wavelet coefficients at
%                       each scale.
% 'nlevel':             Maximum level of the wavelet decomposition.
% 'lengthofsamples':    A vector array containing the length of the samples
%                       at each scale.
% 
% Output variables:
% 
% 'p':      Power Spectral Density at a particular scale.
% 'dp':     Error associated with the above.
% 
%**************************************************************************

p = zeros(nlevel,1);
dp = zeros(nlevel,1);

for k = 1:nlevel
    
    cw = coefs{1,k}(:,1).*coefs{1,k}(:,1);
    p(k) = nanmean(cw);
    dp(k) = nanstd(cw);
    
end

% errors are calculated assuming that we will be plotting on a loglog plot.

p=p.*(dt*2);

dp=dp.*(dt*2);
dp=dp./(sqrt(lengthofsamples)); 
dp=dp./(p.*log(10));   
dp=dp.*2;  

clear('k','cw');

end