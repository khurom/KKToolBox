function [i,zeta,er]=sfexponents(taus,moments,fitline,pushup)

% Needs commenting!!! This is currently using standard errors of the slope
% for calculating confidence intervals.

if nargin<3, fitline = 0; end
if nargin<4, pushup = 0; end

if isrow(taus)>0
    taus=taus';
end

sizemom=size(moments);
moms=sizemom(2)-1;

taus=log10(taus);
moments=log10(moments);

scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)-(scrsz(3)/3) scrsz(4)/5 scrsz(3)/7 scrsz(4)/1.5])
% rect = [left, bottom, width, height]

% moment plots
subplot(2,1,1)

plot(taus,moments(:,11)+pushup,'ks');hold on;
plot(taus,moments(:,21)+2+pushup,'ks');hold on;
plot(taus,moments(:,31)+4+pushup,'ks');hold on;
plot(taus,moments(:,41)+6+pushup,'ks');hold on;
plot(taus,moments(:,51)+8+pushup,'ks');hold on;

title('Structure functions');

xlabel('log_{10}(\tau [secs])'),ylabel('log_{10}(S_m)');

% Ask which points (start and end) the user wants to fit lines to.

st=input('START point of linear fit: ');
en=input('END point of linear fit: ');

carryon=1;

while (carryon==1)

% Ordinary linear regression estimates based on inner points

x=cell(moms,1); y=cell(moms,1); a=cell(moms,1); b=cell(moms,1);

n=en-st+1;

for i=1:1:moms
    
    x{i,1}=[ones(n,1) taus(st:en)]; y{i,1}=moments(st:en,i+1);
    [a{i,1},b{i,1}]=regress(y{i,1},x{i,1},0.32);
%   using 68% confidence intervals i.e. 1 sigma if it was normal, slightly
%   more than 1 if it is a t-statistic as we are using, depending on the
%   fit points.
    
end

clear('x','y','i');

% points for plotting line of best fit

subplot(2,1,1,'replace')

plot(taus,moments(:,11)+pushup,'ks',taus(st:en),...
    taus(st:en).*a{10,1}(2)+a{10,1}(1)+pushup,'r-');hold on;
plot(taus,moments(:,21)+2+pushup,'ks',taus(st:en),...
    taus(st:en).*a{20,1}(2)+a{20,1}(1)+2+pushup,'r-');hold on;
plot(taus,moments(:,31)+4+pushup,'ks',taus(st:en),...
    taus(st:en).*a{30,1}(2)+a{30,1}(1)+4+pushup,'r-');hold on;
plot(taus,moments(:,41)+6+pushup,'ks',taus(st:en),...
    taus(st:en).*a{40,1}(2)+a{40,1}(1)+6+pushup,'r-');hold on;
plot(taus,moments(:,51)+8+pushup,'ks',taus(st:en),...
    taus(st:en).*a{50,1}(2)+a{50,1}(1)+8+pushup,'r-');hold on;

title('Structure functions');

xlabel('log_{10}(\tau [secs])'),ylabel('log_{10}(S_m)');

% Exponents \zeta of structure functions

zeta=ones(moms+1,1);

zeta(1,1)=0.0;

for i=2:1:moms+1
    
    zeta(i,1)=a{i-1,1}(2);
    
end

clear('i');

% Errors on the exponents \zeta using the confidence bounds from the
% regression analysis above

er=ones(moms+1,1);

er(1,1)=0.0;

for i=2:1:moms+1
    
    er(i,1)=max(abs(zeta(i,1)-b{i-1,1}(2,1)),abs(zeta(i,1)-b{i-1,1}(2,2)));
    
end

% Both the upper and lower bound intervals are the same anyway, so this
% above line is a bit redundant

clear('i');

% Plot exponents \zeta against moment order

% Counter to represent moment order in the x-axis of the plots below
i=0:moms;
i=i./10;

% Plotting exponents \zeta against moment order 'i' -- with errors
subplot(2,1,2,'replace')
errorbar(i(1:moms+1),zeta(1:moms+1,1),er(1:moms+1,1),'ks'); hold all;
h=plot(i(1:moms+1),zeta(1:moms+1,1),'ks'); %hold on;

title('Scaling exponents');
xlabel('Moment m'),ylabel('\zeta(m)');

% Ask if the user wants to do a different regression.

st=input('Enter new START point for linear fit (or return to end): ');

if isempty(st)==0
    carryon=1;
    en=input('Enter new END point for linear fit:');
else
    carryon=0;
end

end

% Line fitting

if fitline>0
    
    SlopeRegress(h);
    
end

end