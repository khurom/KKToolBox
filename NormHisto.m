% ------------------------------------------------------------------
% produces a normalised histogram i.e. PDF of data in a vector array
% ------------------------------------------------------------------

function [bins,nc,errnc] = NormHisto(Xin,numbBins,errorflag,scalingflag,plotoption,pushup)

a=isnan(Xin);
[row,col]=find(a==0);
Xin=Xin(row,:);

clear('a','row','col');

% Only takes column vectors. Has a spas attack if you pass a row vector.

% %---------------------------------------------------
% % Summarises some statitical parameters for the data
% %---------------------------------------------------
% % calculate standard deviation
% mean_Xin=mean(Xin)
% var_Xin=var(Xin)
std_Xin=std(Xin);
% kurt_Xin=kurtosis(Xin)
% % -----------------------------


% calculate bin size for uniform bins
binSize=(max(Xin)-min(Xin))./numbBins;
n=hist(Xin,numbBins);
%[n,bins]=hist(X,numbBins);
bins = (min(Xin)+(binSize*0.5)):binSize:(max(Xin)-(binSize*0.5));

% These are uniformaly spaced bins at the moment -- you can
% generalise and change this file later if you wish

area=sum(binSize.*n);

nc=n./area;

% This bit just confirms that the total area is 1 i.e. the
% PDF is correctly normalised
Norm=sum(binSize*nc);

% %% Create figure
% figure1 = figure('PaperPosition',[0.6345 6.345 20.3 15.23],'PaperSize',[20.98 29.68]);
%  
% %% Create axes
% axes1 = axes(...
%   'FontSize',20,...
%   'YMinorTick','on',...
%   'Parent',figure1);
% xlabel(axes1,'x');
% ylabel(axes1,'P(x)');
% box(axes1,'on');
% hold(axes1,'all');

% Zeros will cause problem when taking logs and calculating errors

nc(find(nc==0))=NaN;
n(find(n==0))=NaN;

errnc=1./(n.*(sqrt(n)));    % Error if plotting log density

if scalingflag
    
% % Uncomment the next two lines if you want to get a pdf collapse on \tau
%     nc=nc*(tau.^H);
%     bins=bins./(tau.^H);
    nc=nc*(std_Xin);
    bins=bins/(std_Xin);
    
end

%figure
plot(bins,log10(nc)+pushup,plotoption); hold on;

if errorflag
    
    errorbar((bins),log10(nc)+pushup,errnc,'k'); hold all;
    
end

