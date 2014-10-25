function gaps=GapLabelling(in,threshold,tolerance)
%**************************************************************************
% 
% File name: GapLabelling.m
% Author: K. H. Kiyani
% File created: 21st December 2012
% Last updated: 29th December 2012
% 
% Description:
% ------------
% function which outputs an array with the beginning and end indices of the
% data gaps which are larger than the time variable 'threshold'
% 
% Input variables:
% ----------------
% in: vector of times; assumed to be monotonically increasing
% threshold: threshold period on which to accept gaps, default=1
% tolerance: tolerance value of acceptance within threshold, default=0.2
% 
% Output variables:
% -----------------
% gaps: an array with three columns; 1st row containing the beginning 
%       indices of the gaps in the timeseries; 2nd row containing the end
%       indices; and the last row containing the size of the gap.
% 
%**************************************************************************
% Notes:
% ------
% This adjustment of 'tolerance' is the tolerance we are willing to accept 
% on the labelling i.e. if we had a threshold=1 specified, then a sequence 
% such as 3.9 and 5.1 would produce a gap of 1.2. It is clear in the case 
% of the sequence 3.9 and 5.1 that these are tolerable representations of 4
% and 5 so should be within the threshold of 1; hence we pick a 
% tolerance=0.2. This is what will be used when dealing with ULYSES data 
% for example. 
% 
%**************************************************************************

threshold=threshold+tolerance

a=diff(in);
b=numel(find(a>threshold));
gaps=zeros(b,3);
clear('b');
gaps(:,1)=find(a>threshold);
gaps(:,3)=a(gaps(:,1));
clear('a');
gaps(:,2)=gaps(:,1)+1;

end
