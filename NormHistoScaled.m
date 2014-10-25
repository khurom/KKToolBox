% ------------------------------------------------------------------
% produces a normalised histogram i.e. PDF of data in a vector array
% ------------------------------------------------------------------

function NormHistoScaled(Xin,numbBins,tau,H)

a=isnan(Xin);
[row,col]=find(a==0);
Xin=Xin(row,:);

clear('a');

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
Norm=sum(binSize*nc)


%% Create figure
figure1 = figure('PaperPosition',[0.6345 6.345 20.3 15.23],'PaperSize',[20.98 29.68]);
 
%% Create axes
axes1 = axes(...
  'FontSize',20,...
  'YMinorTick','on',...
  'Parent',figure1);
xlabel(axes1,'x/\sigma(\tau)');
ylabel(axes1,'log_{10} {\sigma}P(x,\tau)');
box(axes1,'on');
hold(axes1,'all');

% Uncomment the next two lines if you want to get a pdf collapse
% according to the scaling exponent tau
%
% nc=nc*(tau.^H);
% bins=bins./(tau.^H);

% Uncomment the next two lines if you want to get a pdf collapse
% according to the standard deviation
%
nc=nc*(std_Xin);
bins=bins/(std_Xin);

% produces a histogram
figure(1)
bar1=plot(bins,log10(nc));

hold('all');





