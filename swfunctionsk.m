% Function to calculate structure functions from wavelet coefficients

function [taus,AllMoments]=swfunctionsk(data,dt,scales,angles,dir,displ)

%% Setting up font sizes for plots
set(0,'defaultaxesfontsize',16);
set(0,'defaulttextfontsize',16);

%% Structure functions

taus=scales.*dt;

AllMoments=ones(length(taus),51);

for n=1:1:length(taus),
     
    tempdiffinc=((1/sqrt(taus(n))).*(abs(data{1,n}((angles{1,dir}{1,n}),1)))).^(0.1);
    
    tempdiff=tempdiffinc;
    
    AllMoments(n,1)=1.0;
    AllMoments(n,2)=nanmean(tempdiff);
    
    for j=3:1:51,
        
        tempdiff=tempdiff.*tempdiffinc;
        AllMoments(n,j)=nanmean(tempdiff);
                
    end

end

%% moment plots

if displ
    
    clf
    
    figure(1);
    plot(log10(taus(:)),log10(AllMoments(:,11)),'.r'); hold on;
    plot(log10(taus(:)),log10(AllMoments(:,21))+2,'.b'); hold on;
    plot(log10(taus(:)),log10(AllMoments(:,31))+4,'.b'); hold on;
    plot(log10(taus(:)),log10(AllMoments(:,41))+6,'.b'); hold on;
    plot(log10(taus(:)),log10(AllMoments(:,51))+8,'.k');
    xlabel('log_{10}(\tau [secs])'),ylabel('log_{10}(S_m)');

end


