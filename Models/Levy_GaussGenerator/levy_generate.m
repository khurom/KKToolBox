% To make a walk Y from an alpha stable symmetric Levy process X
% 

clear all; close all;

randn('state', sum(100*clock));

N=1e6; % Number of sample intervals
alpha=1.4;
V=pi*rand(1,N)-pi/2; % Random variable V uniformly distributed on (-pi/2,pi/2)
W=exprnd(1,1,N); % W is exponentially distributed random variable with mean 1
fac1=sin(alpha*V)./(cos(V).^(1/alpha));
fac2=(cos(V-alpha*V))./W.^((1-alpha)/alpha);
X=fac1.*fac2;

t=60*(1:N);

figure(1)
plot(t,X)
title('Simulation of \alpha -stable motion')
ylabel('\alpha -stable motion, X')
xlabel('time, seconds')

figure(2)
Y=cumsum(X);
plot(t,Y)
grid on
title(['Simulation of Levy walk, Y, of index \alpha=',num2str(alpha)])
ylabel('Y')
xlabel('time, seconds')


% a=[t;Y];
% alphastr=num2str(alpha,'%4.2f');
% idx=findstr(alphastr,'.');
% fid=fopen(strcat('levyout',alphastr(1:idx-1),'Pt',alphastr(idx+1:(numel(alphastr))),'.dat'),'w');
% fprintf(fid,'%f %f\n',a);
% fclose(fid);
