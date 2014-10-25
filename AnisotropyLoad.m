%**************************************************************************
function [BTot,Bmag,TperpAv,Vdir,velmag,dt,samplelength]...
                                =AnisotropyLoad(Bfield,Particles,fs)
%**************************************************************************
% 
% File name: AnisotropyLoad.m
% Author: Khurom H. Kiyani
% File created: 14th November 2013
% Last updated: 1st December 2013
% Updated by: Khurom H. Kiyani
% 
% Description:
% ------------
% This m-file takes in the Cluster Spacecraft file addresses and loads the 
% relvant variables into the workspace; ready for the other functions which 
% will perfomr the wavelet transforms on the data.
% 
% Input variables:
% 
% 'Bfield':     String with the file address for the magnetic field
% 'Particle':   String with the file address for the particle fields
% 'fs':         Sampling frequency of the magnetic field data
% 
% Output variables:
% 
% 'BTot':           B field vector
% 'Bmag':           Average value of the magnitude of B over the interval        
% 'TperpAv':        Average value of the perp temperature over the interval
% 'Vdir':           Averaged unit vector of the direction of the SW velocity 
% 'velmag':         Average value of the SW speed
% 'dt':             Sampling time interval between measurements
% 'samplelength':   Number of measurements in signal
% 
%**************************************************************************

load(Particles);

Vdir=nanmean(CIS(:,3:5),1); velmag=norm(Vdir); Vdir=Vdir./velmag;

TperpAv=nanmean(CIS(:,8)); clear CIS;

load(Bfield);

Bmag=sqrt(BTot(:,2).*BTot(:,2) + BTot(:,3).*BTot(:,3) + BTot(:,4).*BTot(:,4));
Bmag=nanmean(Bmag);

dt=1/fs;

samplelength=length(BTot(:,1));

% If length of data is odd, turn into even numbered sample by getting rid 
% of one point

if mod(samplelength,2)>0
    BTot=BTot(1:end-1,:);
end

samplelength=length(BTot(:,1));


end