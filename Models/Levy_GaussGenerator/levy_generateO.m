%To make a walk Y from an alpha stable process X
% 

clear
close

N=1e6 % Number of sample intervals
alpha=1.9
V=pi*rand(1,N)-pi/2; % Random variable V uniformly distributed on (-pi/2,pi/2)
W=exprnd(1,1,N); % W is exponentially distributed random variable with mean 1
fac1=sin(alpha*V)./((cos(V)).^(1/alpha));
fac2=((cos(V-alpha*V))./W).^((1-alpha)/alpha);
X=fac1.*fac2;
t=1*(1:N);

clear V W fac1 fac2
clear N alpha

Y=cumsum(X);

% clear X;

t=t'; 
% X=X'; 
% std_X=std(X)

%Y=sort(abs(X));
% std_Y=std(Y)

Y=Y'; 
% a=[t Y];
% save levywalk10lastTrial.dat -ASCII a;
 
% clear a

% NormHisto(X,30000);