%**************************************************************************
% m-file for project: Solar wind turbulence magnetic compressibility study
%**************************************************************************
% 
% File name: MagCompressibility.m
% Author: K. H. Kiyani
% File created: 6th August 2010
% Last updated: 11th August 2010
% 
% Description: (needs to be redone)
% ------------
% This m-file estimates the spectral power density of a time
% series using the discrete wavelet decomposition. The power
% density is estimated at scales that are distributed on a dyadic grid
%   y       time series [array]
%   dt      sampling period, default is 1
%   wname   wavelet name, in Matlab syntax. Default is 'db4' for 4th
%           order Daubechies
%   p       power density : [n,1] array
%   scale   corresponding scale axis : [n,1] array
%   dp      standard deviation of P estimated from variability
%
%**************************************************************************
% Notes:
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
% We need to create seperate functions for structure functions, PSD, PDF
% based on different methods, and then plot them all together using
% subplot on matlab; and the output them to EPS files then going onto a
% LaTeX file. 
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
Bmag=mean(Bmag)  

%% load the new data survey

clear all; close all;

cd /Users/khurom/Dropbox/work.plasma/ProjectSurveyDissapationRange/data/Interval_01
load C3_CIS_20060307_1400_1540.mat;

% cd /Users/khurom/Dropbox/work.plasma/ProjectSurveyDissapationRange/data/Interval_02
% load C1_CIS_20070127_1800_1950.mat

% cd /Users/khurom/Dropbox/work.plasma/ProjectSurveyDissapationRange/data/Interval_03
% load C3_CIS_20070127_2120_2320.mat

% cd /Users/khurom/Dropbox/work.plasma/ProjectSurveyDissapationRange/data/Interval_04
% load C3_CIS_20040222_0310_0450.mat

% cd /Users/khurom/Dropbox/work.plasma/ProjectSurveyDissapationRange/data/Interval_13
% load C1_CIS_20010405_2235_2332.mat

% cd /Users/khurom/Dropbox/work.plasma/ProjectSurveyDissapationRange/data/Interval_14
% load C1_CIS_20020211_1919_2029.mat

% cd /Users/khurom/Dropbox/work.plasma/ProjectSurveyDissapationRange/data/Interval_17
% load C1_CIS_20070328_0400_0440.mat

vel=nanmean(CIS(:,3:5),1);

velmag=norm(vel);

vel=vel./norm(vel);

ionTperp=nanmean(CIS(:,8));

clear('CIS');

load C3_BTot_20060307_1400_1540.mat   % Interval_01
% load C1_BTot_20070127_1800_1950.mat   % Interval_02
% load C3_BTot_20070127_2120_2320.mat   % Interval_03
% load C3_BTot_20040222_0310_0450.mat   % Interval_04
% load C1_BTot_20010405_2235_2332.mat   % Interval_13
% load C1_BTot_20020211_1919_2029.mat   % Interval_14
% load C1_BTot_20070328_0400_0440.mat   % Interval_17

Bmag=sqrt(BTot(:,2).*BTot(:,2) + BTot(:,3).*BTot(:,3) + BTot(:,4).*BTot(:,4));
Bmag=mean(Bmag)

cd /Users/khurom/Dropbox/work.plasma/ToolBox/

%% Set parameters needed for UDWT

samplelength=length(BTot(:,1));

% If length of data is odd, turn into even numbered sample by getting rid 
% of one point

if mod(samplelength,2)>0
    BTot=BTot(1:end-1,:);
end

samplelength=length(BTot(:,1));

% dt=1/450; 
dt=1/25; 

wname='coif2';

Lps=7;          %   Low pass filter phase shift for level 1 Coiflet2
Hps=4;          %   High pass filter phase shift for level 1 Coiflet2

[h_0,h_1,~,~] = wfilters(wname);   % set up the filters

% dbworder=10; % The number of zero moments of the Daubechies wavelet
% 
% wname=strcat('db',num2str(dbworder));  
% 
% [h_0,h_1] = daubcqf(2*dbworder);    % set up the filters

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
    
% Correct for linear phase shift in wavelet coefficients at each level

    for j=1:nlevel

    shiftfac=Hps*(2^(j-1));
    
    for l=1:1:j
    
        shiftfac=shiftfac + Lps*(2^(l-2))*((l-2)>=0) ;
    
    end
    
    swd(:,j)=circshift(swd(:,j),shiftfac);
    clear('shiftfac','l');
    
    end
    
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
    
    load(tempfiles{m-1},'A','swd');
    
    SWT_parts{1,m} = cell(1,2);
    
%     only use this if you run out of RAM
%     SWT_parts{1,m}{1,1} = A(((pads/2)+1):(end-(pads/2)));

    SWT_parts{1,m}{1,1} = A(:,((pads/2)+1):(end-(pads/2)));
        clear('A');
        
%     SWT_parts{1,m}{1,2} = D(:,((pads/2)+1):(end-(pads/2)));
%         clear('D');

    SWT_parts{1,m}{1,2} = swd(:,((pads/2)+1):(end-(pads/2)));
        clear('swd');
        
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

%% Getting the indices of the arrays at different angles

tic

anglebins=[60.0 70.0 80.0 90.0];

absize=length(anglebins);

angleindices=cell(1,absize);

for m=1:absize

angleindices{1,m}=cell(1,nlevel);

startangle=anglebins(m)-5.0;
endangle=anglebins(m)+5.0;

for k = 1:nlevel
    
    angleindices{1,m}{1,k}=find((coefsBVangle{1,k}(:,1))>startangle & ...
        (coefsBVangle{1,k}(:,1))<endangle);
    
end

clear('k');

end

clear('m');

toc

%% Calculating power spectral density for par with errors

tic

p_par = zeros(nlevel,absize);
dp_par = zeros(nlevel,absize);

for m=1:absize

for k = 1:nlevel
    
    cw = coefspar{1,k}((angleindices{1,m}{1,k}),1)...
        .*coefspar{1,k}((angleindices{1,m}{1,k}),1);
    p_par(k,m) = (nanmean(cw)).*(dt*2);
    dp_par(k,m) = (nanstd(cw)).*(dt*4)...
        ./(sqrt(length(angleindices{1,m}{1,k})));
    dp_par(k,m)=dp_par(k,m)./(p_par(k,m).*log(10));   
    
end

clear('k','cw');

end

clear('m');

toc

%% Calculating power spectral density for perpA with errors

tic

p_perpA = zeros(nlevel,absize);
dp_perpA = zeros(nlevel,absize);

for m=1:absize

for k = 1:nlevel
    
    cw = coefsperpA{1,k}((angleindices{1,m}{1,k}),1)...
        .*coefsperpA{1,k}((angleindices{1,m}{1,k}),1);
    p_perpA(k,m) = (nanmean(cw)).*(dt*2);
    dp_perpA(k,m) = (nanstd(cw)).*(dt*4)...
        ./(sqrt(length(angleindices{1,m}{1,k})));
    dp_perpA(k,m)=dp_perpA(k,m)./(p_perpA(k,m).*log(10));   
    
end

clear('k','cw');

end

clear('m');

toc

%% Calculating power spectral density for perpB with errors

tic

p_perpB = zeros(nlevel,absize);
dp_perpB = zeros(nlevel,absize);

for m=1:absize

for k = 1:nlevel
    
    cw = coefsperpB{1,k}((angleindices{1,m}{1,k}),1)...
        .*coefsperpB{1,k}((angleindices{1,m}{1,k}),1);
    p_perpB(k,m) = (nanmean(cw)).*(dt*2);
    dp_perpB(k,m) = (nanstd(cw)).*(dt*4)...
        ./(sqrt(length(angleindices{1,m}{1,k})));
    dp_perpB(k,m)=dp_perpB(k,m)./(p_perpB(k,m).*log(10));   
    
end

clear('k','cw');

end

clear('m');

toc

%% Calculating power spectral density for perpMag with errors

tic

p_perpMag = zeros(nlevel,absize);
dp_perpMag = zeros(nlevel,absize);

for m=1:absize

for k = 1:nlevel
    
    cw = coefsperpMag{1,k}((angleindices{1,m}{1,k}),1)...
        .*coefsperpMag{1,k}((angleindices{1,m}{1,k}),1);
    p_perpMag(k,m) = (nanmean(cw)).*(dt*2);
    dp_perpMag(k,m) = (nanstd(cw)).*(dt*4)...
        ./(sqrt(length(angleindices{1,m}{1,k})));
    dp_perpMag(k,m)=dp_perpMag(k,m)./(p_perpMag(k,m).*log(10));   
    
end

clear('k','cw');

end

clear('m');

toc

%% Calculate magnetic compressibility

tic

magcom = zeros(nlevel,absize);
dmagcom = zeros(nlevel,absize);

for m=1:absize

for k = 1:nlevel

magcom(k,m) = (p_par(k,m)./(p_perpMag(k,m)+p_par(k,m)));
dmagcom(k,m) = ((dp_par(k,m)).*(dp_par(k,m))).*((p_par(k,m)).*(p_par(k,m)));
dmagcom(k,m) = dmagcom(k,m) + ...
    ((dp_perpMag(k,m)).*(dp_perpMag(k,m))).*((p_perpMag(k,m)).*(p_perpMag(k,m)));
dmagcom(k,m) = dmagcom(k,m)./...
    ((p_par(k,m)+p_perpMag(k,m)).*(p_par(k,m)+p_perpMag(k,m)));
dmagcom(k,m) = sqrt(((dp_par(k,m)).*(dp_par(k,m))) + dmagcom(k,m)); 
   

end

clear('k');

end

clear('m');

toc

%% Alternative magnetic compressibility

tic

magcom = zeros(nlevel,absize);
dmagcom = zeros(nlevel,absize);

for m=1:absize

for k = 1:nlevel

        cwperp = coefsperpMag{1,k}((angleindices{1,m}{1,k}),1)...
        .*coefsperpMag{1,k}((angleindices{1,m}{1,k}),1);
    
        cwpar = coefspar{1,k}((angleindices{1,m}{1,k}),1)...
        .*coefspar{1,k}((angleindices{1,m}{1,k}),1);
    
    cw = cwpar./(cwpar+cwperp);
    
    clear('cwpar','cwperp');
    
    magcom(k,m) = (nanmean(cw));
    dmagcom(k,m) = (nanstd(cw)).*2./(sqrt(length(angleindices{1,m}{1,k})));
    dmagcom(k,m) = dmagcom(k,m)./((magcom(k,m)).*log(10));

end

clear('k','cw');

end

clear('m');

toc

%% Fouad magnetic compressibility for linear theory

tic

magcom = zeros(nlevel,absize);
dmagcom = zeros(nlevel,absize);

for m=1:absize

for k = 1:nlevel

        cwperp = coefsperpA{1,k}((angleindices{1,m}{1,k}),1)...
        .*coefsperpA{1,k}((angleindices{1,m}{1,k}),1);
    
        cwpar = coefspar{1,k}((angleindices{1,m}{1,k}),1)...
        .*coefspar{1,k}((angleindices{1,m}{1,k}),1);
    
    cw = cwpar./(cwpar+cwperp);
    
    clear('cwpar','cwperp');
    
    magcom(k,m) = (nanmean(cw));
    dmagcom(k,m) = (nanstd(cw)).*2./(sqrt(length(angleindices{1,m}{1,k})));
    dmagcom(k,m) = dmagcom(k,m)./((magcom(k,m)).*log(10));

end

clear('k','cw');

end

clear('m');

toc


%% Translating scale to frequency

scale = 2.^[1:nlevel]';
frequency = scal2frq(scale,wname,dt);
wavenumber = frequency.*(2*pi/(velmag*(1e3)));

aconstant = 9.4848e-7;    % sqrt(k_B * m_p / e^2)
iongyrorad = aconstant*sqrt(ionTperp*(1e6)/((Bmag*(1e-9))^2));

krhoi = wavenumber.*iongyrorad;


%% plotting PSD
    
    set(0,'defaultaxesfontsize',20);
    set(0,'defaulttextfontsize',20);

    figure(1)
    
    for k=1:absize
    
    plot(log10(frequency(:)),log10(p_par(:,k)),'bo-','LineWidth',2,...
        'MarkerFaceColor','w','MarkerEdgeColor','b'); hold all;
    errorbar(log10(frequency),log10(p_par(:,k)),dp_par(:,k),'k');
%     plot(log10(frequency),log10(p_perpA),'g.-');
%     errorbar(log10(frequency),log10(p_perpA),dp_perpA,'k');
%     plot(log10(frequency),log10(p_perpB),'k.-');
%     errorbar(log10(frequency),log10(p_perpB),dp_perpB,'k');

    plot(log10(frequency(:)),log10(p_perpMag(:,k)),'ro-','LineWidth',2,...
        'MarkerFaceColor','w','MarkerEdgeColor','r');
%     plot(log10(krhop),log10(p_par),'ro-','LineWidth',2,...
    errorbar(log10(frequency),log10(p_perpMag(:,k)),dp_perpMag(:,k),'k');
    
    end

    axis tight;
    
    ylabel('log_{10} PSD [nT^{2} Hz^{-1}]')    
    
    clear('k');

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
%     plot(log10(krhoi(2:end-1)),log10(p_perpMag(2:end-1)),'o','LineWidth',2,...
%         'MarkerFaceColor','none','MarkerEdgeColor','none','Parent',ax2);
%     
%     axis tight;
%--------------------------------------------------------------------------
    
    grid on

    figure(2)
    
    for k=1:absize
    
%     errorbar(log10(krhoi(:)),log10(magcom(:,k)),dmagcom(:,k),'k');
    
    hold all;
    
    plot(log10(krhoi(:)),log10(magcom(:,k)),...
        'MarkerFaceColor',[1 1 1],...
        'MarkerEdgeColor',[0 0 0],...
        'MarkerSize',8,...
        'Marker','o',...
        'LineWidth',2,...
        'Color',[1 0 0]);
    
    end
    
    clear('k')    
    
    grid on
    
    xlabel('k\rho_{i}')
    ylabel('C_{||}')
    
    axis tight
    
%% Delete crappy files

delete('UDWTtempfile_2_coif2.mat'); delete('UDWTtempfile_3_coif2.mat'); delete('UDWTtempfile_4_coif2.mat');
