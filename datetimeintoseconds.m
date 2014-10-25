function time=datetimeintoseconds(in)
%**************************************************************************
% 
% File name: datetimeintoseconds.m
% Author: K. H. Kiyani
% File created: 2nd January 2013
% Last updated: 15th October 2014
% 
% Description:
% ------------
% function which changes a time vector from the format 'dd-mm-yyyy
% HH:MM:SS.FFF', into serial date number with a format 'SS'.
% 
% Input variables:
% ----------------
% in: vector of times; assumed to be monotonically increasing
% 
% Output variables:
% -----------------
% time: 
% 
% Auxillary files:
% ----------------
% datenum.m (from MATLAB)
% 
% Notes:
% ------
% 1. To convert back to a normal date-time, here is an example for the
% first ten entries:
% datestr(time(1:10)./(24*60*60),'dd-mmm-yyyy HH:MM:SS:FFF')
% 
% 2. You will need to comment your original text file with '%' so that this
% routine ignores the stuff that you don;t need to read. 
% 
% 3. You might need to put more '%*f' in below if you have lots of
% variables. Then again you might not. Try it out and see if it works.
% 
%**************************************************************************

tic

fid=fopen(in,'r');
C=textscan(fid,'%s %s %*f %*f %*f %*f %*f %*f %*f','CommentStyle','%');
D=strcat(C{1,1},'SPACE',C{1,2});
D=regexprep(D,'\SPACE',' ');
fclose(fid);
E=datenum(D,'dd-mm-yyyy HH:MM:SS.FFF');
time=E.*24*60*60;

toc

end