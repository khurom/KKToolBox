function data=RegulariseGrid2(in,threshold,tolerance)
%**************************************************************************
% 
% File name: RegulariseGrid2.m
% Author: K. H. Kiyani
% File created: 21st December 2012
% Last updated: 2nd January 2013
% 
% Description:
% ------------
% function which outputs an array with the beginning and end indices of the
% data gaps which are larger than the time variable 'threshold'
% 
% Input variables:
% ----------------
% in: data array with first column assumed to be time, monotonically 
%     increasing.
% threshold: threshold period on which to accept gaps
% tolerance: tolerance value of acceptance within threshold
% 
% Output variables:
% -----------------
% data: original data put onto a regularised grid using very conservative 
%       interpolation; with NaNs labelling the gaps.
% 
% Auxillary files:
% ----------------
% GapLabelling.m
% 
%**************************************************************************
% Notes:
% ------
% Make sure to choose 'threshold' appropriately, as some data (e.g. ULYSES)
% are a bit shit in time labelling i.e time stamps are not exact integers. 
%
%**************************************************************************

newtime=[in(1,1):1:in(end,1)]';
% New time vector (interpolation onto 1 second); for a more coarse-grained 
% timeeseries with 'threshold' time seperation use: 
% newtime=[in(1,1):threshold:in(end,1)]';

gaps=GapLabelling(in(:,1),threshold,tolerance);

% Put a half second time stamp and the beginning and end of the gap and 
% give it the value NaN. Doesn't matter if things overlap here, as long as 
% the NaNs stay  within the gap.

gaplabelsStart=[in(gaps(:,1),1)+0.5,nan(numel(gaps(:,1)),size(in,2)-1)];
gaplabelsEnd=[in(gaps(:,2),1)-0.5,nan(numel(gaps(:,2)),size(in,2)-1)];

clear('gaps','threshold','tolerance');

gapslabeltot=[gaplabelsStart;gaplabelsEnd];

clear('gaplabelsStart','gaplabelsEnd');

in=[in;gapslabeltot];

clear('gapslabeltot');

in=sortrows(in,1);

data=interp1(in(:,1),in(:,2:size(in,2)),newtime);
% Currently using linear interpolation. However, if wanting to use nearest
% neighbour interpolation, then use the following line instead:
% data=interp1(in(:,1),in(:,2:5),newtime,'nearest');

data=[newtime,data];

clear('newtime');

end