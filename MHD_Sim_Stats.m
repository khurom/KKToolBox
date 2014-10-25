
%% Load data

% cd /media/DATA/Work.SpacePlasma/DataCubes/;
% cd /media/Extra/DataCubes/;
cd /home/khurom/Desktop/PROJECT_OF_THE_MILENIUM/Strong/

%% load perp data

% clear all; close all;

% load bPerpStrong.mat; 

mwt_perp=test00_b00; clear test00_b00;

% load mwt_perp0.mat; mwt_perp=mwt_perp0; clear mwt_perp0;
% load mwt_perp1.mat; mwt_perp=mwt_perp1; clear mwt_perp1;

%% load parallel data

load mwt_para0.mat; mwt_para=mwt_para0; clear mwt_para0;
% load mwt_para1.mat; mwt_para=mwt_para1; clear mwt_para1;

%% Prepare perp data TWO PERP

component=4;

% First perp component
bx=mwt_perp(1:end/2,component);

% Second perp component
% bx=mwt_perp((end/2)+1:end,1);

sizebx=size(bx);

clear mwt_perp component;

% Bx_perpA=reshape(bx,512,64,32); % For EMHD 512 cube

Bx_perpA=reshape(bx,1536,192,8);

clear bx;

% Calculate different spatial scales 'l'

base=1.2; maxexp=35;

l=unique(floor(base.^(0:1:maxexp))); 

clear base maxexp;

%% Prepare perp data ONE PERP

component=4;

% First perp component
bx=mwt_perp(1:end,component);

sizebx=size(bx);

clear mwt_perp component;

% Bx_perpA=reshape(bx,2048,192,8); %Without
% Bx_perpA=reshape(bx,2048,192,96); %With
% Bx_perpA=reshape(bx,1024,128,64); %Strong EMHD

clear bx;

% Calculate different spatial scales 'l'

base=1.2; maxexp=35;

l=unique(floor(base.^(0:1:maxexp))); 

clear base maxexp;

%% prepare parallel data

component=4;

bx=mwt_para(1:end,component);

sizebx=size(bx);

clear mwt_para component;

% Bx_perpA=reshape(bx,512,64,16); % For EMHD 512 cube

Bx_perpA=reshape(bx,128,192,192);

clear bx;

% Calculate different spatial scales 'l'

base=1.2; maxexp=23;

l=unique(floor(base.^(0:1:maxexp)));

clear base maxexp;

%% Fluctuations

Fluctuations=cell(size(l));

for comp=1:1:numel(l)
    
    temp=circshift(Bx_perpA,-l(comp));
    temp(end-l(comp)+1:end,:,:)=NaN;
    Fluctuations{comp}=reshape(temp-Bx_perpA,sizebx);
    clear temp;
    
end

clear m;

%% Use overlapping for calculating fluctuations and populate cell.

Fluctuationsx=cell(size(l));

for comp=1:1:numel(l)
    
    temp=circshift(Bx_perpA,-l(comp));
    temp(end-l(comp)+1:end,:,:)=NaN;
    Fluctuationsx{comp}=reshape(temp-Bx_perpA,sizebx);
    clear temp;
    
end

clear m;

Fluctuationsy=cell(size(l));

for comp=1:1:numel(l)
    
    temp=circshift(By_perpA,-l(comp));
    temp(end-l(comp)+1:end,:,:)=NaN;
    Fluctuationsy{comp}=reshape(temp-By_perpA,sizebx);
    clear temp;
    
end

clear m;

Fluctuationsz=cell(size(l));

for comp=1:1:numel(l)
    
    temp=circshift(Bz_perpA,-l(comp));
    temp(end-l(comp)+1:end,:,:)=NaN;
    Fluctuationsz{comp}=reshape(temp-Bz_perpA,sizebx);
    clear temp;
    
end

clear m;

%% Magnitude of fluctuations

Fluctuations=cell(size(l));

for comp=1:1:numel(l)
    
Fluctuations{comp}=sqrt(Fluctuationsx{comp}.*Fluctuationsx{comp}+Fluctuationsy{comp}.*Fluctuationsy{comp}+Fluctuationsz{comp}.*Fluctuationsz{comp});
    
end

clear m;

%% Plot the rescaled PDF for the fluctuations

for comp=1:1:10   %1:1:numel(l)
    
    NormHisto(Fluctuations{1,comp},300,0,1,'r');
    
end

%% Plot the rescaled PDF for the fluctuations

for comp=8:1:20   %1:1:numel(l)
    
    NormHisto(differ{1,comp},300,0,1,'r');
    
end

%% PDFs without scaling 

NormHisto(Fluctuationsx{1,4},400,0,0,'b');
NormHisto(Fluctuationsy{1,4},400,0,0,'g');
NormHisto(Fluctuationsz{1,4},400,0,0,'r');

%% These conditions are given in percentages -- change as you will.
% conditions=[0, 0.01, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0];
conditions=[0];

ConditionedMoments=cell(1,length(conditions));

tic

for comp=1:1:length(conditions)
    [momentorder,ConditionedMoments{comp}]=sfunctions(OutlierThresh(Fluctuations,...
        conditions(comp)),l,0);
end

toc

clear('m');

%% Structure functions

[momentorder,Moments]=sfunctions(Fluctuations,l,5,1);

clear('m');

%% Exponents zeta(p)

% sfexponents((1/1536).*l,Moments,1,0);
sfexponents((1/1024).*l,Moments,1,0);

%% Kurtosis

kurt(Fluctuations,(1/512).*l,1);

%% ****************************************************************
% WAVELET STUFF
% *************************************************************************

%% EMHD Strong

clear all; close all;

load test00_b00.147_slyces_perp;

data=test00_b00; clear test00_b00;

dt=1/1024;
nx=1024; ny=1024; nz=1024;


data=data(:,4);

data=reshape(data,[],nz/32);

% data=data(:,1);

% spect=cell(32,1);
% 
% for k=1:1:32
%     
%     [f,P]=spectro(data(:,k),2^11,dt,'b',0);
%     spect{k,1}=P'; clear P;
%     
% end

%% With/Without

clear all; close all;

% % With
load test00_b00.366_slyces_perp;
load test00_u00.366_slyces_perp;

% Without
% load test00_b00.305_slyces_perp;
% load test00_u00.305_slyces_perp;

data=test00_u00-test00_b00; clear test00_b00 test00_u00;
% data=test00_b00; clear test00_b00

dt=1/1536; 
nx=1536; ny=1536; nz=128;

[N,D]=rat(logspace(log10(0.5),0,10));
N=fliplr(N(2:end-1)); D=fliplr(D(2:end-1));

% [nx*(ny/32)]*(nz/32)

% Calculate magnitude of Z_plus
data=sqrt(data(:,1).*data(:,1)+data(:,2).*data(:,2)+data(:,3).*data(:,3));

% data=data(:,4);

data=reshape(data,[],nz/32);

% data=data(:,1);

%% Resample

m=4;

% dt=dt/(N(m)/D(m));    % Undersample
dt=dt/(D(m)/N(m));  % Oversample

data=resample(data(11:end-10,:),N(m),D(m));

%% resample using spline interpolation [much better than above]

resamplefac=logspace(log10(0.5),0,10);

m=3;

dt=dt*resamplefac(m)

xnew=1:resamplefac(m):73728;

ynew=interp1([1:1:73728],data,xnew,'spline');
%%
plot([1:1:73728],data(:,1),'.-b');
hold on
plot(xnew,ynew(:,1),'.r')


%% Set parameters needed for UDWT

samplelength=length(data(:,1));

% If length of data is odd, turn into even numbered sample by getting rid 
% of one point

if mod(samplelength,2)>0
    data=data(1:end-1,:);
end

samplelength=length(data(:,1));

wname='coif2';

Lps=7;          %   Low pass filter phase shift for level 1 Coiflet2
Hps=4;          %   High pass filter phase shift for level 1 Coiflet2

[h_0,h_1,~,~] = wfilters(wname);   % set up the filters

nlevel = wmaxlev(samplelength,wname)   % maximum nr of levels

% edge extension mode set to periodic extension by default with this
% routine in the rice toolbox.

pads=(2^(nextpow2(samplelength))) - samplelength;   % for edge extension

%% Do the UDWT decompositon and reconstruction
    
% Gets the data size up to the next power of 2 due to UDWT restrictions
% Although periodic extension is used for the wavelet edge handling we are
% getting the data up to the next power of 2 here by extending the data
% sample with a constant value
    
coefsTot=cell(1,nlevel);

for m=1:1:(nz/32)

y=[data(1,m).*ones(pads/2,1); data(:,m); data(end,m).*ones(pads/2,1)];

% Decompose the signal using the UDWT
    
tic
[swa,swd,L] = mrdwt(y,h_0,nlevel);
toc
    
clear('y','L','swa');
    
% *************************************************************************
% VERY IMPORTANT: LINEAR PHASE SHIFT CORRECTION
% *************************************************************************
% Correct for linear phase shift in wavelet coefficients at each level. The
% formula for the shift has been taken from Walden's paper, or was 
% constructed by myself (can't exactly remember) -- but it is verified and 
% correct.
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

swd=swd';

swd = swd(:,((pads/2)+1):(end-(pads/2)));

% Getting rid of the edge effects

filterlength = length(h_0);

coefs=cell(1,nlevel);

for j=1:1:nlevel
    
    extra = (2^(j-2))*filterlength;

    swd(j,1:extra)=NaN;
    swd(j,end-extra+1:end)=NaN;
    
    coefs{1,j}=swd(j,:)';
    
    clear('extra');
    
end

clear('swd','j','filterlength');

for j=1:1:nlevel
    
    coefsTot{1,j}=[coefsTot{1,j};coefs{1,j}];
    
end

clear('coefs','j');

end

clear('Hps','Lps'); clear('pads'); clear('h_0','h_1');

%% Structure functions

scale = 2.^(1:nlevel)';
frequency = scal2frq(scale,wname,dt);

[taus,AllMoments]=swfunctions(coefsTot,dt,scale,0);

%% Add SFs

load SFTot.mat

tausTot=[taus;tausTot]; AllMomentsTot=[AllMoments; AllMomentsTot];

[tausTot,ind]=sort(tausTot);
AllMomentsTot=AllMomentsTot(ind,:);

save('SFTot.mat','tausTot','AllMomentsTot');

%% Zeta exponents

sfexponents(taus,AllMoments,0,0);
% sfexponents(AllMoments(:,21),AllMoments,0,5);

%% PDF

%% Calculate PDFs at dfferent scales for the Price Coastline

figure(2); hold on;

for k=3:1:6
    
    NormHisto(coefsMinus{1,k},200,0,1,'.-b'); hold on;

end

xlabel('x/\sigma(\tau)');
ylabel('log_{10} \sigma(\tau) PDF(x,\tau)');

%% Wavelet Power spectra

for k = 1:nlevel

    cw = coefsTot{k}.*coefsTot{k};
    p(k) = nanmean(cw);
    
end

% Converting to power/energy spectral density
p=p.*(dt*2);

    loglog(frequency,p,'b');
    grid on
    xlabel('log_{10} frequency (Hz)')
    ylabel('log_{10} PSD (signal units)^{2}Hz^{-1}')
    title(['nr of levels = ',int2str(nlevel),...
        '   Wavelet is Coiflet2'])
    
%% Spectral simulation power spectrum

load test00_spp_split.366

sppSpec=reshape(test00_spp_split(:,3),769,[]);
totalspec=nansum(sppSpec(:,2:50),2);
loglog([0:1:768],128.*totalspec,'r');
loglog([0:1:768],128.*sppSpec(:,1),'k');

%% Spectra of all slices

[F,a1,era]=spectro(data(:,1),2^12,dt,'b'); hold on
[F,a2,era]=spectro(data(:,2),2^12,dt,'r');
[F,a3,era]=spectro(data(:,3),2^12,dt,'b');
[F,a4,era]=spectro(data(:,4),2^12,dt,'r');
a=(a1+a2+a3+a4)./4;
loglog((F(2:end)),(a(2:end)),'b');

%% fit nonlinear regression model for the co-dimension parameter estimation

modelfun = @(b,x)x(:,1)/8 + b(1) - b(1)*((1 - (3/(b(1)*4))).^(x(:,1)/2));
beta0 = [1.2];
mdl = fitnlm(zeta(:,1),zeta(:,2),modelfun,beta0)

%% multifractals
hold on;
p=0:0.1:5;
czero=1.0419;
zeta_p = p./8 + czero - czero*((1 - (3/(czero*4))).^(p/2));
% zeta_p = p./8 +1 -(1/4).^(p/2);
% zeta_p = (p./16) +1 -((1/4).^(p/4));
% zeta_p = p./4 +1 -(1/2).^(p/2);
% zeta_p = p.*(0.9/2) - p.*(3/8) + 1 - (1/4).^(p./2);
plot(p,zeta_p,'r');

%% fractals

hold on;
p=0:0.1:5;
% zeta_p = p.*(0.9/2) - p.*(3/8) + 1 
zeta_p=p.*0.666;
plot(p,zeta_p,'r');
