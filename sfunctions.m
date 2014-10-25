% Function to calculate structure functions

function [scales,AllMoments]=sfunctions(data,scales,maxmom,displ)

% displ takes the values '1' for plots, or '0' for no plots.

% This code needs to be vectorised. 

% Setting up font sizes for plots
set(0,'defaultaxesfontsize',16);
set(0,'defaulttextfontsize',16);

% Structure functions

maxmom=(maxmom*10) + 1;

AllMoments=ones(length(scales),maxmom);

for n=1:1:length(scales),
     
    tempdiffinc=(abs(data{1,n})).^(0.1);
    
    tempdiff=tempdiffinc;
    
    AllMoments(n,1)=1.0;
    AllMoments(n,2)=nanmean(tempdiff);
    
    for j=3:1:maxmom,
        
        tempdiff=tempdiff.*tempdiffinc;
        AllMoments(n,j)=nanmean(tempdiff);
                
    end

end

% moment plots

if displ
    
    clf
    
    figure(1);
    plot(log10(scales(:)),log10(AllMoments(:,11)),'.r'); hold on;
    plot(log10(scales(:)),log10(AllMoments(:,21))+2,'.b'); hold on;
    plot(log10(scales(:)),log10(AllMoments(:,31))+4,'.b'); hold on;
    plot(log10(scales(:)),log10(AllMoments(:,41))+6,'.b'); hold on;
    plot(log10(scales(:)),log10(AllMoments(:,51))+8,'.k');
    xlabel('log_{10}(\tau [secs])'),ylabel('log_{10}(S_m)');

end


