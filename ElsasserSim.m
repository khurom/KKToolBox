%% Load data

cd ~/Desktop/PROJECT_OF_THE_MILENIUM/Without/;
% cd ~/Desktop/PROJECT_OF_THE_MILENIUM/With/;

%% load perp data

clear all; close all;

load zperpWout.mat;

% load zperpWith.mat;

%% Prepare perp data ONE PERP

component=4;

% Zplus
zp=zp(1:end,component);

sizezp=size(zp);

zp=reshape(zp,1536,192,8);

% Zminus
zm=zm(1:end,component);

zm=reshape(zm,1536,192,8);

clear component;

% Calculate different spatial scales 'l'

base=1.2; maxexp=35;

l=unique(floor(base.^(0:1:maxexp))); 

clear base maxexp;

%% Elssaser fluctuations

FluctuationsZp=cell(size(l));

for m=1:1:numel(l)
    
    temp=circshift(zp,-l(m));
    temp(end-l(m)+1:end,:,:)=NaN;
    FluctuationsZp{m}=reshape(temp-zp,sizezp);
    clear temp;
    
end

clear m;

FluctuationsZm=cell(size(l));

for m=1:1:numel(l)
    
    temp=circshift(zm,-l(m));
    temp(end-l(m)+1:end,:,:)=NaN;
    FluctuationsZm{m}=reshape(temp-zm,sizezp);
    clear temp;
    
end

clear m;

FluctuationsZ=cell2mat(FluctuationsZm).*cell2mat(FluctuationsZp);
Epp=FluctuationsZ.*cell2mat(FluctuationsZp);
Epm=FluctuationsZ.*cell2mat(FluctuationsZm);

% clear FluctuationsZm FluctuationsZp
clear FluctuationsZ

Epp=mat2cell(Epp,[2359296],ones(1,31));
Epm=mat2cell(Epm,[2359296],ones(1,31));

%% Plot the rescaled PDF for the fluctuations

figure(1)
hold on;

for m=1:1:numel(l)
    
    NormHisto(Epp{1,m},300,0,1,'r');
    
end

figure(2)
hold on;

for m=1:1:numel(l)
    
    NormHisto(Epm{1,m},300,0,1,'b');
    
end


%% Structure Functions for Epp

[momentorderEpp,MomentsEpp]=sfunctions(Epp,l,5,1);

clear('m');

%% Structure Functions for Epm

[momentorderEpm,MomentsEpm]=sfunctions(Epm,l,5,1);

clear('m');

%% Exponents zeta(p)

sfexponents((1/1536).*l,MomentsEpp,0,0);

%% Exponents zeta(p)

sfexponents((1/1536).*l,MomentsEpm,0,0);

%% PDFs for Elssaser variables seperately

figure(1)
hold on;

for m=1:1:numel(l)
    
    NormHisto(FluctuationsZm{1,m},300,0,1,'r');
    
end

figure(2)
hold on;

for m=1:1:numel(l)
    
    NormHisto(FluctuationsZp{1,m},300,0,1,'b');
    
end

%% Structure Functions for Zm

[momentorderZm,MomentsZm]=sfunctions(FluctuationsZm,l,6,1);

clear('m');

%% Exponents zeta(p)

sfexponents((1/1536).*l,MomentsZm,0,0);

%% Structure Functions for Zp

[momentorderZp,MomentsZp]=sfunctions(FluctuationsZp,l,6,1);

clear('m');


%% Exponents zeta(p)

sfexponents((1/1536).*l,MomentsZp,0,0);

%%
% \langle |(\deltaz^{-})^2 \deltaz^{+}|^{m}  \rangle
