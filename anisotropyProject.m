%**************************************************************************
function [Bpar,BperpA,BperpB,BperpMag]=...
    anisotropyProject(parDir,perpA,perpB,details)
%**************************************************************************
% 
%   Filename:       Anisotropy2.m
%   Author:         K. H. Kiyani
%   Last updated:   26/05/2010
% 
% Description:
% ------------
% Projection function which takes the details and the unit vector from the 
% approximations and projects the details onto the unit vector. This will 
% churn out the parallel and perpendicular details.
% 
% Inputs:-
% 
% parDir:
% perpA:
% perpB:
% details:
% 
% Outputs:-
% 
% Bpar:
% BperpA:
% BperpB:
% BperpMag:
% 
%**************************************************************************

Bpar=details(:,1).*parDir(:,1) + details(:,2).*parDir(:,2)...
    + details(:,3).*parDir(:,3);

BperpA=details(:,1).*perpA(:,1) + details(:,2).*perpA(:,2)...
    + details(:,3).*perpA(:,3);

BperpB=details(:,1).*perpB(:,1) + details(:,2).*perpB(:,2)...
    + details(:,3).*perpB(:,3);

BperpMag=sqrt(BperpA.*BperpA + BperpB.*BperpB);


% Another way of getting Bperp magnitude
% **************************************
% Bperp2=[details(:,1) - Bpar(:).*parDir(:,1),details(:,2)...
%     - Bpar(:).*parDir(:,2),details(:,3) - Bpar(:).*parDir(:,3)];
% 
% Bperp2=sqrt(Bperp(:,1).*Bperp(:,1) + Bperp(:,2).*Bperp(:,2)...
%     + Bperp(:,3).*Bperp(:,3));
% 
% max(BperpMag - Bperp2)   % check for difference
