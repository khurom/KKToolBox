%   Filename:           Anisotropy.m
%   Author:               K. H. Kiyani
%   Last updated:     25/02/2010
%Projection function which takes the details and the unit vector from the approximations
% and projects the details onto the unit vector. This will need to churn out the parallel and 
% perpendicular details.


%%  Decompose to parallel and perp to mean field in interval

Bmean=mean(BFGM,1);

Bmag=sqrt(Bmean(:,2).*Bmean(:,2)+Bmean(:,3).*Bmean(:,3)+Bmean(:,4).*Bmean(:,4));

parDir=Bmean(2:4)./Bmag;

Bpar=[BFGM(:,1) , BFGM(:,2:4)*parDir'];
Bperp=[BFGM(:,1) , BFGM(:,2:4) - Bpar(:,2)*parDir];
Bperp=[Bperp , sqrt(Bperp(:,2).*Bperp(:,2) + Bperp(:,3).*Bperp(:,3) + Bperp(:,4).*Bperp(:,4))];

BSCpar=[BSC(:,1) , BSC(:,2:4)*parDir'];
BSCperp=[BSC(:,1) , BSC(:,2:4) - BSCpar(:,2)*parDir];
BSCperp=[BSCperp , sqrt(BSCperp(:,2).*BSCperp(:,2) + BSCperp(:,3).*BSCperp(:,3) + BSCperp(:,4).*BSCperp(:,4))];

B_ANIS=[Bpar, Bperp(:,5)];
BSC_ANIS=[BSCpar, BSCperp(:,5)];

%%  Power spectral density
figure(6);
[fs,ps,ers]=spectro(B_ANIS(:,2),16384,1/67,'-k'); hold on;
%
% figure(7);
[fsc,psc,ersc]=spectro(BSC_ANIS(:,2),16384*4,1/450,'-k'); hold on;

% STAFF-SC noise floor
% plot([0 1 2],[-5 -7.1 -9]);

% irf_psd(BSC(:,1:4),16384*4); hold on;
% irf_psd(BFGM(:,1:4),16384*20,67,[]); hold on;


%% Output the BSC and BFGM data to an ascii file


fid=fopen('BSTAFFSC4_20070120_ANIS_A.dat','w');
fprintf(fid,'%f %f %f\n',BSC_ANIS');
fclose(fid);

fid=fopen('BFGMSC4_20070120_ANIS_A.dat','w');
fprintf(fid,'%f %f %f\n',B_ANIS');
fclose(fid);




