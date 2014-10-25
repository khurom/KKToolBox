%**************************************************************************
% m-file for project: Solar wind turbulence anisotropy study
%**************************************************************************
% 
% File name: AnisotropyStudyUDWT.m
% Author: K. H. Kiyani
% File created: 7th June 2010
% Last updated: 12th March 2012
% 
% Description:
% ------------
% This m-file decomposes the vector magnetic field fluctuations onto the
% vector background magnetic field using the undecimated discrete wavelet
% decomposition. The scale-by-scale 'details' and 'approximations' provide
% the scale-by-scale fluctuations and background field, respectively. This
% type of anisotropy study is historically called 'variance anisotropy'.
% This is basically the anisotropy study of magnetic field quantities
% rather than the wavevector 'k_i' values that it depends upon.
% 
% The following quantities are estimated:

% * the power spectral density parallel and perpendicular to the background
% magnetic field estimated at scales that are distributed on a dyadic grid
% 2^n:{1 ... m};
% 
% * the ratio of the parallel and perpendicular power spectral densities;
% 
% * the structure functions and resulting exponents of the fluctuations
% parallel and perpendicular to the background magnetic field estimated at 
% scales that are distributed on a dyadic grid 2^n:{1 ... m};
% 
% * the probablity density functions for the fluctuations parallel and 
% perpendicular to the background magnetic field estimated at scales that 
% are distributed on a dyadic grid 2^n:{1 ... m}.
% 
% Auxillary files: (need to be renamed)
% ----------------
% anisotropy1.m
% anisotropyProject.m
% 
% Related files:
% --------------
% MagCompressibility.m
% MagCompressibilityGlobalField.m
%
%**************************************************************************
% 
% Additional Notes:
% -----------------
% Remember that when you are calculating structure functions you will have
% to normalise the coefficients by 'square root of the scale' you are 
% looking at. This is because the wavelets coefficients are normalised 
% w.r.t. this factor due to the L^2 Norm. This does not effect calculation 
% of the power as the power is related to the L^2 norm.
% 
% Example code to plot structure functions:
% plot(log10(scale),log10((p./(scale))./(2*dt)),'.r')
% 
% Boundary effects on finite size samples on the interval:
% Check the Advanced concepts in MATLAB help or where to truncate at each
% level appropriately.
% Example code when you are checking approximations and boundary effects:
% plot(y); hold on; plot(A(16,:),'r','LineWidth',2);
% 
%**************************************************************************
% 
% To do list:
% -----------
% 1./   We need to create seperate functions for structure functions, PSD, 
%       PDF based on different methods, and then plot them all together 
%       using subplot on matlab; and the output them to EPS files then 
%       going onto a LaTeX file. 
% 2./   This whole m-file needs to be made into a function which outputs
%       the above into a PDF file ready for viewing. 
%       
%       It must take the following data inputs:
%       
%       a./ Velocity field vector
%       b./ Perpendicular ion temperature
%       c./ Magnetic field vector
% 
%       and provide the following outputs:
% 
%       a./ Scaling exponents with errors
%       b./ PDF's fitted with t-location scale distributions ??
% 
%% load 20070130 interval data

clear all; close all;

cd ~/khuromOLD/SOLARStuff/NewClusterWork/Data/BurstModeIntervals/20070130_0000

c_load VCh1;

vel=nanmean(VCh1(:,2:4),1);

velmag=norm(vel);

vel=vel./norm(vel);

clear('VCh1');

c_load TperpCh1;

ionTperp=mean(TperpCh1(:,2));

clear('TperpCh1');

cd ~/khuromOLD/SOLARStuff/NewClusterWork/project2Anisotropy/data/interval20070130

load BMerged.mat

Bmag=sqrt(BTot(:,2).*BTot(:,2) + BTot(:,3).*BTot(:,3) + BTot(:,4).*BTot(:,4));
Bmag=mean(Bmag);

%% load 20070120ALL interval data

clear all; close all;

cd ~/khuromOLD/SOLARStuff/NewClusterWork/Data/BurstModeIntervals/20070120_1200

c_load VCh1;

vel=nanmean(VCh1(:,2:4),1);

vel=vel./norm(vel);

clear('VCh1');

cd ~/khuromOLD/SOLARStuff/NewClusterWork/project2Anisotropy/data/interval20070120ALL

load BMerged.mat

%% load 20070120A interval data

clear all; close all;

cd ~/khuromOLD/SOLARStuff/NewClusterWork/Data/BurstModeIntervals/20070120_1200

c_load VCh1;

vel=nanmean(VCh1(1:1061,2:4),1);

vel=vel./norm(vel);

clear('VCh1');

cd ~/khuromOLD/SOLARStuff/NewClusterWork/project2Anisotropy/data/interval20070120A

load BMerged.mat

%% load 20070120B interval data

clear all; close all;

cd ~/khuromOLD/SOLARStuff/NewClusterWork/Data/BurstModeIntervals/20070120_1200

c_load VCh1;

vel=nanmean(VCh1(1239:end,2:4),1);

vel=vel./norm(vel);

clear('VCh1');

cd ~/khuromOLD/SOLARStuff/NewClusterWork/project2Anisotropy/data/interval20070120B

load BMerged.mat

%% Set parameters needed for UDWT

samplelength=length(BTot(:,1));

% If length of data is odd, turn into even numbered sample by getting rid 
% of one point

if mod(samplelength,2)>0
    BTot=BTot(1:end-1,:);
end

samplelength=length(BTot(:,1));

dt=1/450;

wname='coif2';

Lps=7;          %   Low pass filter phase shift for level 1 Coiflet2
Hps=4;          %   High pass filter phase shift for level 1 Coiflet2

[h_0,h_1,~,~] = wfilters(wname);   % set up the filters

nlevel = wmaxlev(samplelength,wname)   % maximum nr of levels

% edge extension mode set to periodic extension by default with this
% routine in the rice toolbox.

pads=(2^(nextpow2(samplelength))) - samplelength;   % for edge extension

%% Do the UDWT decompositon and reconstruction

tempfiles = cell(1,3);

for m=2:1:4
    
% Gets the data size up to the next power of 2 due to UDWT restrictions
% Although periodic extension is used for the wavelet edge handling we are
% getting the data up to the next power of 2 here by extending the data
% sample with a constant value
    
    y=[BTot(1,m).*ones(pads/2,1); BTot(:,m); BTot(end,m).*ones(pads/2,1)];
    
% Decompose the signal using the UDWT
    
    tic
    [swa,swd,L] = mrdwt(y,h_0,nlevel);
    toc
    
    clear('y','L');
    
% Reconstruct all the approximations and details at all levels

    tic
    
    mzero = zeros(size(swd));
    A = mzero;
    A(:,nlevel) = mirdwt(swa,mzero,h_0,nlevel);
    
    D = mzero;
    
    for i = 1:nlevel
        
        swcfs = mzero;
        swcfs(:,i) = swd(:,i);
        D(:,i) = mirdwt(zeros(length(swa),1),swcfs,h_0,nlevel);
        
    end
    
    clear('swa','i');
    
    for j=nlevel-1:-1:1
        
        A(:,j) = A(:,j+1) + D(:,j+1);
    
    end
    
    clear('j');
    
    toc

% *************************************************************************
% VERY IMPORTANT: LINEAR PHASE SHIFT CORRECTION
% *************************************************************************
% Correct for linear phase shift in wavelet coefficients at each level. No
% need to do this for the low-pass filters approximations as they will be
% reconstructed and the shift will automatically be reversed. The formula
% for the shift has been taken from Walden's paper, or has been made up by
% me (can't exactly remember) -- but it is verified and correct.
% *************************************************************************

    for j=1:nlevel

    shiftfac=Hps*(2^(j-1));
    
    for l=1:1:j
    
        shiftfac=shiftfac + Lps*(2^(l-2))*((l-2)>=0) ;
    
    end
    
    swd(:,j)=circshift(swd(:,j),shiftfac);
    clear('shiftfac','l');
    
    end
    
% *************************************************************************
    
% Save variables to external file on hardisk
    
    A=A';
    D=D';
    swd=swd';
    
    tempfiles{m-1}=strcat('UDWTtempfile_',num2str(m),'_',wname,'.mat');
    save(tempfiles{m-1},'A','D','swd');
    clear('A','D','swd');
    
end
    
clear('m');

%% temporary for convenience if you had to start up due to restart

for m=2:1:4
    tempfiles{m-1}=strcat('UDWTtempfile_',num2str(m),'_',wname,'.mat');
end

clear('m');

%% Put all the files together into a cell structure

tic

SWT_parts = cell(1,4);

SWT_parts{1,1} = BTot(:,1);

for m=2:1:4
    
    load(tempfiles{m-1},'A','D');
    
    SWT_parts{1,m} = cell(1,2);
    
%     only use this if you run out of RAM
%     SWT_parts{1,m}{1,1} = A(((pads/2)+1):(end-(pads/2)));

    SWT_parts{1,m}{1,1} = A(:,((pads/2)+1):(end-(pads/2)));
        clear('A');
        
    SWT_parts{1,m}{1,2} = D(:,((pads/2)+1):(end-(pads/2)));
        clear('D');

%     SWT_parts{1,m}{1,2} = swd(:,((pads/2)+1):(end-(pads/2)));
%         clear('swd');
        
end

toc

clear('m','pads','tempfiles','BTot');

%% Getting rid of the edge effects; to keep edges skip this section

filterlength = length(h_0);

clear('h_0','h_1');

for j=1:1:nlevel
    
    extra = (2^(j-2))*filterlength; % give some reasoning for this eq
    
    for m=2:1:4
        
        % for approximations
        SWT_parts{1,m}{1,1}(j,1:extra)=NaN;
        SWT_parts{1,m}{1,1}(j,end-extra+1:end)=NaN;
        
        % for details
        SWT_parts{1,m}{1,2}(j,1:extra)=NaN;
        SWT_parts{1,m}{1,2}(j,end-extra+1:end)=NaN;
        
    end
    
    clear('m','extra');
    
end

clear('j','filterlength');

%% Decompose detail coefficients at all levels into perp and par components
% Have to program in this particular way to overcome memory problems due to
% many variables of very large data sets.

%***************************

coefsBVangle=cell(1,nlevel);

tic

for i=1:1:nlevel

    [~,~,~,BVangle] = anisotropy1([SWT_parts{1,2}{1,1}(i,:)',...
        SWT_parts{1,3}{1,1}(i,:)',SWT_parts{1,4}{1,1}(i,:)'],vel);
    
    coefsBVangle{1,i}=BVangle;
    
    clear('BVangle');
    
end

toc

save('BVangle.mat','coefsBVangle');

clear('i','coefsBVangle');

%***************************

%***************************

coefspar=cell(1,nlevel);

tic

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

toc

save('parcoefs.mat','coefspar');

clear('i','coefspar');

%***************************

%***************************

coefsperpA=cell(1,nlevel);

tic

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

toc

save('perpAcoefs.mat','coefsperpA');

clear('i','coefsperpA');


%***************************

%***************************

coefsperpB=cell(1,nlevel);

tic

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

toc

save('perpBcoefs.mat','coefsperpB');

clear('i','coefsperpB');

%***************************

%***************************

coefsperpMag=cell(1,nlevel);

tic

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

toc

save('perpMagcoefs.mat','coefsperpMag');

clear('i','coefsperpMag');

%***************************

clear('SWT_parts');

% Load temporary files from hardisk.

load('BVangle.mat'); load('parcoefs.mat'); load('perpAcoefs.mat');
load('perpBcoefs.mat'); load('perpMagcoefs.mat');

% Delete temporary files from hardisk.

delete('BVangle.mat'); delete('parcoefs.mat'); delete('perpAcoefs.mat');
delete('perpBcoefs.mat'); delete('perpMagcoefs.mat');

%% Calculating length of vector at each wavelet stage, excluding NaNs

lengthofsamples=zeros(nlevel,1);

for j=1:1:nlevel
    
    a=isnan(coefspar{1,j});
    [row,col]=find(a==0);
    b=coefspar{1,j}(row,:);
    lengthofsamples(j)=length(b);
    
    clear('a','b','row','col');
    
end

%% Calculating power spectral density for par with errors

p_par = zeros(nlevel,1);
dp_par = zeros(nlevel,1);

for k = 1:nlevel
    
    cw = coefspar{1,k}(:,1).*coefspar{1,k}(:,1);
    p_par(k) = nanmean(cw);
    dp_par(k) = nanstd(cw);
    
end

p_par=p_par.*(dt*2);

dp_par=dp_par.*(dt*2);
dp_par=dp_par./(sqrt(lengthofsamples)); 
dp_par=dp_par./(p_par.*log(10));   
dp_par=dp_par.*2;  

clear('k','cw');

%% Calculating power spectral density for perpA with errors

p_perpA = zeros(nlevel,1);
dp_perpA = zeros(nlevel,1);

for k = 1:nlevel
    
    cw = coefsperpA{1,k}(:,1).*coefsperpA{1,k}(:,1);
    p_perpA(k) = nanmean(cw);
    dp_perpA(k) = nanstd(cw);
    
end

p_perpA=p_perpA.*(dt*2);

dp_perpA=dp_perpA.*(dt*2);
dp_perpA=dp_perpA./(sqrt(lengthofsamples)); 
dp_perpA=dp_perpA./(p_perpA.*log(10));   
dp_perpA=dp_perpA.*2;  

clear('k','cw');

%% Calculating power spectral density for perpB with errors

p_perpB = zeros(nlevel,1);
dp_perpB = zeros(nlevel,1);

for k = 1:nlevel
    
    cw = coefsperpB{1,k}(:,1).*coefsperpB{1,k}(:,1);
    p_perpB(k) = nanmean(cw);
    dp_perpB(k) = nanstd(cw);
    
end

p_perpB=p_perpB.*(dt*2);

dp_perpB=dp_perpB.*(dt*2);
dp_perpB=dp_perpB./(sqrt(lengthofsamples)); 
dp_perpB=dp_perpB./(p_perpB.*log(10));   
dp_perpB=dp_perpB.*2;  

clear('k','cw');

%% Calculating power spectral density for perpMag with errors

p_perpMag = zeros(nlevel,1);
dp_perpMag = zeros(nlevel,1);

for k = 1:nlevel
    
    cw = coefsperpMag{1,k}(:,1).*coefsperpMag{1,k}(:,1);
    p_perpMag(k) = nanmean(cw);
    dp_perpMag(k) = nanstd(cw);
    
end

p_perpMag=p_perpMag.*(dt*2);

dp_perpMag=dp_perpMag.*(dt*2);
dp_perpMag=dp_perpMag./(sqrt(lengthofsamples)); 
dp_perpMag=dp_perpMag./(p_perpMag.*log(10));   
dp_perpMag=dp_perpMag.*2;  

clear('k','cw');


%% Translating scale to frequency

scale = 2.^[1:nlevel]';
frequency = scal2frq(scale,wname,dt);
% wavenumber = frequency.*(2*pi/(velmag*(1e3)));
% 
% aconstant = 9.4848e-7;    % sqrt(k_B * m_p / e^2)
% iongyrorad = aconstant*sqrt(ionTperp*(1e6)/((Bmag*(1e-9))^2));
% 
% krhoi = wavenumber.*iongyrorad;


%% plotting PSD
    
    set(0,'defaultaxesfontsize',20);
    set(0,'defaulttextfontsize',20);

    figure(1)
    subplot(2,1,1)
    
    plot(log10(frequency(2:end)),log10(p_par(2:end)),'bo-','LineWidth',2,...
        'MarkerFaceColor','w','MarkerEdgeColor','b'); hold all;
%     errorbar(log10(frequency),log10(p_par),dp_par,'k');

    plot(log10(frequency(2:end)),log10(p_perpA(2:end)),'g.-'); hold all;
%     errorbar(log10(frequency),log10(p_perpA),dp_perpA,'k');
    plot(log10(frequency(2:end)),log10(p_perpB(2:end)),'k.-');
%     errorbar(log10(frequency),log10(p_perpB),dp_perpB,'k');

    plot(log10(frequency(2:end)),log10(p_perpMag(2:end)),'ro-','LineWidth',2,...
        'MarkerFaceColor','w','MarkerEdgeColor','r');
%     errorbar(log10(frequency),log10(p_perpMag),dp_perpMag,'k');
    
    axis tight;
    
    ylabel('log_{10} PSD [nT^{2} Hz^{-1}]')    
    
    
%     axis square

%--------------------------------------------------------------------------
% Uncomment this bit if you want to put an axis on top with krhoi
% 
%     ax1 = gca;
%     
%     ax2 = axes('Position',get(ax1,'Position'),...
%            'XAxisLocation','top',...
%            'YAxisLocation','right',...
%            'Color','none',...
%            'XColor','k','YColor','k');
%     
%     
%     hold all;
%        
%     plot(log10(krhoi(2:end)),log10(p_perpMag(2:end)),'o','LineWidth',2,...
%         'MarkerFaceColor','none','MarkerEdgeColor','none','Parent',ax2);
%     
%     axis tight;     % absolutely necessary - do not leave out, otherwise
%                     % second axis goes a bit awry, as everything else is
%                     % axis tight also.
%--------------------------------------------------------------------------
    
    grid on
    
%     axis square

    subplot(2,1,2)
    
  
  plot(log10(frequency(2:end)),log10(p_par(2:end)./p_perpMag(2:end)),...
        'MarkerFaceColor',[1 1 1],...
        'MarkerEdgeColor',[0 0 0],...
        'MarkerSize',8,...
        'Marker','o',...
        'LineWidth',2,...
        'Color',[0 0 1]);
    
    grid on
    
    xlabel('log_{10} spacecraft frame frequency [Hz]')
    ylabel('log_{10} PSD_{||/\perp}')
      
    axis tight
    
%     axis square
    
%% Structure functions
tic
[tausPar,AllMomentsPar]=swfunctions(coefspar,dt,scale,0);
toc
tic
[tausPerp,AllMomentsPerp]=swfunctions(coefsperpMag,dt,scale,0);
toc

%% Zeta exponents

% sfexponents(tausPar,AllMomentsPar,0)
sfexponents(tausPerp,AllMomentsPerp,15)

%% PDF
