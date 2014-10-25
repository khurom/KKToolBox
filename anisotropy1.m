%**************************************************************************
function [parDir,perpA,perpB,BVangle] =anisotropy1(Bmean,vel)
%**************************************************************************
% 
% Filename:     Anisotropy1.m
% Author:       K. H. Kiyani
% File created: March 2010
% Last updated: 26/05/2010
% 
% Description:
% ------------
% Calculates the mean magnetic field directional unit vector 'parDir'; the 
% two perpendicular unit vectors 'perpA' and 'perpB' to this using the
% solar wind velocity unit vector 'vel'; and finally the angle 'BVangle' of
% 'parDir' to 'vel'.
% 
% Inputs:-
% 
% Bmean: Mean background magnetic field
% vel: Solar wind velocity unit vector
% 
% Outputs:-
% 
% parDir: unit vector parallel to magnetic field direction
% perpA: unit vector = crossproduct of 'parDir' and 'vel'
% perpB: unit vector = crossproduct of 'parDir' and 'perpA'
% BVangle: angle in degrees between 'parDir' and 'vel'
% 
%**************************************************************************

Bmag=sqrt(Bmean(:,1).*Bmean(:,1) + Bmean(:,2).*Bmean(:,2)...
    + Bmean(:,3).*Bmean(:,3));

parDir1=Bmean(:,1)./Bmag;
parDir2=Bmean(:,2)./Bmag;
parDir3=Bmean(:,3)./Bmag;

parDir=[parDir1,parDir2,parDir3];

clear('parDir1','parDir2','parDir3','Bmean','Bmag');

% angle between mean solar wind velcity unit vector and mean magnetic field
% unit vector

cosBVangle = parDir(:,1).*vel(:,1) + parDir(:,2).*vel(:,2)...
    + parDir(:,3).*vel(:,3);

BVangle = acos(cosBVangle);
BVangle = BVangle.*(360/(2*pi)); % conversion from radians to degrees

clear('cosBVangle');

% create new orthonormal basis from parDir, parDir X Vel
% and parDir X (parDir X Vel). Here are the additional perp vectors. 

perpA1 = parDir(:,2).*vel(:,3) - parDir(:,3).*vel(:,2);
perpA2 = parDir(:,3).*vel(:,1) - parDir(:,1).*vel(:,3);
perpA3 = parDir(:,1).*vel(:,2) - parDir(:,2).*vel(:,1);

% We need to normalise by magnitude here as the resultant perpA vector will
% not be a unit vector unless we do so. This is needed because 'vel' and
% 'pardir', even though they are unit vectors, are not necessarily
% orthogonal to each other and hence a cross product will produce something
% which is not exactly norm=1 (which is expected for a unit vector).

perpAmag=sqrt(perpA1.*perpA1 + perpA2.*perpA2 + perpA3.*perpA3);

perpA1 = perpA1./perpAmag;
perpA2 = perpA2./perpAmag;
perpA3 = perpA3./perpAmag;

perpA=[perpA1,perpA2,perpA3];

clear('perpA1','perpA2','perpA3','perpAmag');

% no need to normalise by magnitude here as the two vectors being crossed
% are exactly orthogonal; hence this will produce an exact unit vector.

perpB1 = parDir(:,2).*perpA(:,3) - parDir(:,3).*perpA(:,2);
perpB2 = parDir(:,3).*perpA(:,1) - parDir(:,1).*perpA(:,3);
perpB3 = parDir(:,1).*perpA(:,2) - parDir(:,2).*perpA(:,1);

perpB=[perpB1,perpB2,perpB3];

clear('perpB1','perpB2','perpB3');
