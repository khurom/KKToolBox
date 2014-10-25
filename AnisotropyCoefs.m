%**************************************************************************
function [coefsBVangle,coefspar,coefsperpA,coefsperpB,coefsperpMag,...
          coefsmag,lengthofsamples]=AnisotropyCoefs(SWT_parts,nlevel,vel)
%**************************************************************************
% 
% File name: AnisotropyCoefs.m
% Author: Khurom H. Kiyani
% File created: 1st December 2013
% Last updated: 3rd December 2013
% Updated by: Khurom H. Kiyani
% 
% Description:
% ------------
% This function decompose the fluctuations at all levels into perpendicular 
% and parallel components to the background magnetic field. It uses a local
% scale-dependent orthonormal basis constructed from the the local background
% magnetic field and the direction of average solar wind velocity over the 
% interval being studies. It also calculates the angles between the local
% background magnetic field and the average solar wind velocity direction, 
% which is taken as a proxy for the r or k vector spatial measurement 
% direction. 
% 
% Input variables:
% 
% 'SWT_parts':  A cell containing, for each component of the signal, the
%               decomposition of the signal into fluctuations and
%               background fields at each scale -- decomposed using the
%               UDWT. Edge effects are also removed.
% 'nlevel':     Maximum level of the wavelet decomposition.
% 'vel':        Average unit direction vector of solar wind velocity over
%               the interval being analysed.
% 
% Output variables:
% 
% 'coefsBVangle':   At each scale the local angle between the local magnetic
%                   field and the average solar wind velocity direction unit
%                   vector.
% 'coefspar':       At each scale the fluctuations projected in the parallel
%                   direction to the local background magnetic field.
% 'coefsperpA':     At each scale the fluctuations projected in one of the
%                   perpendicular directions to the local background 
%                   magnetic field.
% 'coefsperpB':     At each scale the fluctuations projected in the other
%                   perpendicular direction to the local background magnetic
%                   field.
% 'coefsperpMag':   At each scale the fluctuations projected in the total
%                   perpendicular plane to the local background magnetic field.
% 'coefsmag':     At each scale the total magnitude of the fluctuations.
% 'lengthofsamples':At each scale the length of the vector of fluctuations
%                   excluding NaNs.
% 
% Notes:
% ------
% We have to program in this particular way with the multiple writing to 
% the hardisk, inorder to overcome memory problems due to many variables of
% very large data sets. You pribably don;t have to do this with newer
% machines.
% 
% Dependencies/Auxillary files:
% -----------------------------
% anisotropy1.m
% anisotropyProject.m
% 
%**************************************************************************


%**************************************************************************
% Angle between local background magnetic field and the average direction 
% of the solar wind velocity over the interval 
%**************************************************************************

coefsBVangle=cell(1,nlevel);

for i=1:1:nlevel

    [~,~,~,BVangle] = anisotropy1([SWT_parts{1,2}{1,1}(i,:)',...
        SWT_parts{1,3}{1,1}(i,:)',SWT_parts{1,4}{1,1}(i,:)'],vel);
    
    coefsBVangle{1,i}=BVangle;
    
    clear('BVangle');
    
end

save('BVangle.mat','coefsBVangle');

clear('i','coefsBVangle');

%**************************************************************************
% Parallel fluctuations
%**************************************************************************

coefspar=cell(1,nlevel);

for i=1:1:nlevel
    
% Using an external function, construct orthonormal basis comprising of 
% background field direction, and two perpendicular vectors using the solar
% wind plasma flow velocity direction.
% This function will also output the angle between the magnetic field and
% the background solar wind plasma flow velocity direction which we can
% assume to be constant.

    [parDir,perpA,perpB,~] = anisotropy1([SWT_parts{1,2}{1,1}(i,:)',...
        SWT_parts{1,3}{1,1}(i,:)',SWT_parts{1,4}{1,1}(i,:)'],vel);
    
% Project fluctuations given by wavelet details, parallel and perp to the
% scale dependent background field direction given by 'parDir'
    
    [Bpar,~,~,~] = anisotropyProject(parDir,perpA,perpB,...
        [SWT_parts{1,2}{1,2}(i,:)',...
        SWT_parts{1,3}{1,2}(i,:)',SWT_parts{1,4}{1,2}(i,:)']);
    
    clear('parDir','perpA','perpB');
    
    coefspar{1,i}=Bpar;
    
    clear('Bpar');
    
end

save('parcoefs.mat','coefspar');

clear('i','coefspar');

%**************************************************************************
% Perpendicular fluctuations in one of the perpendicular directions
%**************************************************************************

coefsperpA=cell(1,nlevel);

for i=1:1:nlevel
    
    [parDir,perpA,perpB,~] = anisotropy1([SWT_parts{1,2}{1,1}(i,:)',...
        SWT_parts{1,3}{1,1}(i,:)',SWT_parts{1,4}{1,1}(i,:)'],vel);
    
    [~,BperpA,~,~] = anisotropyProject(parDir,perpA,perpB,...
        [SWT_parts{1,2}{1,2}(i,:)',...
        SWT_parts{1,3}{1,2}(i,:)',SWT_parts{1,4}{1,2}(i,:)']);
    
    clear('parDir','perpA','perpB');
    
    coefsperpA{1,i}=BperpA;
    
    clear('BperpA');
    
end

save('perpAcoefs.mat','coefsperpA');

clear('i','coefsperpA');

%**************************************************************************
% Perpendicular fluctuations in the other perpendicular direction
%**************************************************************************

coefsperpB=cell(1,nlevel);

for i=1:1:nlevel
    
    [parDir,perpA,perpB,~] = anisotropy1([SWT_parts{1,2}{1,1}(i,:)',...
        SWT_parts{1,3}{1,1}(i,:)',SWT_parts{1,4}{1,1}(i,:)'],vel);
    
    [~,~,BperpB,~] = anisotropyProject(parDir,perpA,perpB,...
        [SWT_parts{1,2}{1,2}(i,:)',...
        SWT_parts{1,3}{1,2}(i,:)',SWT_parts{1,4}{1,2}(i,:)']);
    
    clear('parDir','perpA','perpB');
    
    coefsperpB{1,i}=BperpB;
    
    clear('BperpB');
    
end

save('perpBcoefs.mat','coefsperpB');

clear('i','coefsperpB');

%**************************************************************************
% Magnitude of perpendicular fluctuations
%**************************************************************************

coefsperpMag=cell(1,nlevel);

for i=1:1:nlevel
    
    [parDir,perpA,perpB,~] = anisotropy1([SWT_parts{1,2}{1,1}(i,:)',...
        SWT_parts{1,3}{1,1}(i,:)',SWT_parts{1,4}{1,1}(i,:)'],vel);
    
    [~,~,~,Bperpmag] = anisotropyProject(parDir,perpA,perpB,...
        [SWT_parts{1,2}{1,2}(i,:)',...
        SWT_parts{1,3}{1,2}(i,:)',SWT_parts{1,4}{1,2}(i,:)']);
    
    clear('parDir','perpA','perpB');
    
    coefsperpMag{1,i}=Bperpmag;
    
    clear('Bperpmag');
    
end

save('perpMagcoefs.mat','coefsperpMag');

clear('i','coefsperpMag');

%***************************

% Load temporary files from hardisk.

load('BVangle.mat'); load('parcoefs.mat'); load('perpAcoefs.mat');
load('perpBcoefs.mat'); load('perpMagcoefs.mat');

% Delete temporary files from hardisk.

delete('BVangle.mat'); delete('parcoefs.mat'); delete('perpAcoefs.mat');
delete('perpBcoefs.mat'); delete('perpMagcoefs.mat');


%**************************************************************************
% Total magnitude of fluctuations
%**************************************************************************

coefsmag=cell(1,nlevel);

for m=1:1:nlevel
   
    coefsmag{m}=sqrt(coefspar{m}.*coefspar{m}+coefsperpMag{m}.*coefsperpMag{m});
    
end

clear('m');

%**************************************************************************
% Calculating length of vector at each wavelet stage, excluding NaNs
%**************************************************************************

lengthofsamples=zeros(nlevel,1);

for j=1:1:nlevel
    
    a=isnan(coefspar{1,j});
    [row,~]=find(a==0);
    b=coefspar{1,j}(row,:);
    lengthofsamples(j)=length(b);
    
    clear('a','b','row');
    
end

end