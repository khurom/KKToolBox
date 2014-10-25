%**************************************************************************
function [taus,differ]=DifferOlap(data,dt,base,minExp,maxExp)
%**************************************************************************
% 
% File name: DifferOLap.m
% Author: Khurom H. Kiyani
% File created: 21st April 2014
% Last updated: 21st April 2014
% Updated by: Khurom H. Kiyani
% 
% Description:
% ------------
% This function calculates the increments or 'fluctuations' of a data
% signal for a range of scales. The fluctuations are calculated with
% overlapping increments. One should be wary of this regarding
% correlations, but we will go ahead and do it here because everyone else
% does. 
% 
% Input variables:
% ----------------
% 'data':   The data vector assuming a one-dimensional column vector.
% 'dt':     The time sampling of the data vector.
% 'base':   The base of the geometrically progressing scales.
% 'minExp': The exponent defining the minimum scale of interest.
% 'maxExp': The exponent defining the maximum scale of interest.
% 
% Output variables:
% -----------------
% 'taus':   The scale of interest.
% 'differ': A cell array containing the fluctuations (increments) at all
%           the scales of interest.
% 
% Notes:
% ------
% This routine can easily be vectorised using repmat - do it at some point
% and stop being lazy.
% 
%**************************************************************************

if nargin<2, dt=1; end
if nargin<3, base=1.2; end
if nargin<4, minExp=0; end
if nargin<5, maxExp=10; end

dataSize=numel(data);

taus=(base).^(minExp:maxExp); 
clear base minExp maxExp;
taus=floor(taus); taus=unique(taus);

differ=cell(1,numel(taus)); 

index=1:dataSize;
index=index';

for n=1:numel(taus),
  
    indexn=index+taus(n);
       
    diffdata=data( indexn(1:end-taus(n)) ) - data( index(1:end-taus(n)) );
    
    differ{n}=diffdata;

end

taus=taus.*dt; taus=taus';

end


