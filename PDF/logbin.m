
% This script takes in data, and then produces a logarithmically binned
% histogram:
% [midpoints,histsum,edges] =logbin(data,logbase,binNum)
% where the inputs are 
% data: the input data from some experiment which will be binned
% logbase: the log base used in the binning
% binNum: the number of bins used in the data reduction process.
% 
% The outputs are:
% 
% midpoints := the centre of each log bin
% histsum := the number of counts in each bin
% edges := a 2D vector containing the left and right edge values of each
% bin.
% 
% One should note that this script is a slightly modified version of
%  'logbinfromlinbin.m' which can be found here:
% http://www-personal.umich.edu/~ladamic/courses/si614w06/matlab/index.html
% 
% Modified by LPS April, 2010.
%  The rounding in the original file ('logbinfromlinbin.m') created some
% extra unnecessary bins, this was fixed by the author.
% 
function [midpoints,histsum,edges] =logbin(data,logbase,binNum)

a=data;
b=ones(length(data),1);
if (logbase <=1 )
    disp('warning: logbase should be > 1. exiting...');
    return;
end

[a1,ia] = sort(a);
a2 = b(ia);

max2=log(ceil(max(a1)+1))/log(logbase);
step=max2/binNum;
tmp = 0:step:max2;

x = (logbase).^(tmp);

widths=(x(2:length(x))-x(1:(length(x)-1)));
midpoints=(widths./2)+x(1:(length(x)-1));
edges=[x(1:(length(x)-1)) ; x(2:length(x))];

histsum = zeros(size(midpoints));

i = 1;
j = 1;

while (i <= length(a1)) 
   if (a1(i) < x(j+1)) 
      histsum(j) = histsum(j) + a2(i);
   else
      j = j + 1;
      histsum(j) = histsum(j) + a2(i);
   end
   i = i+1;
end
   
histsum = histsum./widths;
