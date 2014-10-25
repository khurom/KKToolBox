%**************************************************************************
% 
% Function:     RemoveNaN.m
% 
% Description:  This function removes all NaNs from a single column/row
%               data vector.
% 
% Inputs:       'in' is the input data vector.
% 
% Outputs:      'out' is the data vector after all NaNs have been removed.
% 
% Author:       Khurom Kiyani
% Date:         28th February 2013
% Last updated: 28th February 2013
% 
% Copyright 2011 Imperial College London. All rights reserved.
% 
%**************************************************************************

function [out]=RemoveNaN(in)

in(isnan(in))=[];
out=in;

end