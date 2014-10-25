% I should put some comments here

function [taus,variance,differ,AllMoments]=NEWdiffer(data,dt,minExp,maxExp)

dataSize=numel(data);

differ=[];
taus=[]; 
sigmas=[];
AllMoments=[];

%% Main loop

for n=minExp:maxExp,

tau=(1.2).^n;
tau=floor(tau);
taus=[taus tau];

end

taus=unique(taus);

index=1:dataSize;
index=index';

for n=1:numel(taus),
  
    indexn=index+taus(n);
       
    diffdata=data( indexn(1:end-taus(n)) ) - data( index(1:end-taus(n)) );
    
    sigmas=[sigmas kurtosis(diffdata)];
    
    tempdiffinc=abs(diffdata).^(0.1);
    
    tempdiff=tempdiffinc;
    
    moment=[1.0, nanmean(tempdiff)];
    
    for j=0.2:0.1:5,
        
        tempdiff=tempdiff.*tempdiffinc;
        moment=[moment, nanmean(tempdiff)];
        
        
    end
    
    AllMoments=[AllMoments;  moment];
    
    differ=[differ,{diffdata}];

% [nc,bins,sigma]=NormHistoNoPlot(diffdata,300);
% 
% % produces a histogram
% figure1=semilogy(bins,nc,'Color',[0 0 0.5625]);
% 
% sigmas=[sigmas sigma];

end

taus=taus.*dt; taus=taus';

sigmas=sigmas';
variance=sigmas;%.*sigmas;


figure(1);
semilogx((taus),(variance),'.b');

% figure(2);
% loglog((1./taus),(variance),'.r')

%% moment plots
figure(3);
plot(log10(taus(:,1)),log10(AllMoments(:,11)),'.r'); hold on;
plot(log10(taus(:,1)),log10(AllMoments(:,21)),'.b'); hold on;
plot(log10(taus(:,1)),log10(AllMoments(:,31)),'.b'); hold on;
plot(log10(taus(:,1)),log10(AllMoments(:,41)),'.b'); hold on;
plot(log10(taus(:,1)),log10(AllMoments(:,51)),'.k');

