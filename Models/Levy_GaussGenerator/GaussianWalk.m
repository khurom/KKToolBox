% To make a walk Y from a Brownian process X
% 
 clear all; close all;

N=1e6; % Number of sample intervals

t=60*(1:N);

randn('state', sum(100*clock));

X = normrnd(0,0.01,N,1);

X=X';

figure(1)
plot(t,X)
title('Simulation of Brownian increments')
ylabel('Brownian increments, X')
xlabel('time, seconds')

figure(2)
Y=cumsum(X);
plot(t,Y)
grid on
title('Simulation of Brownian walk, Y')
ylabel('Y')
xlabel('time, seconds')

% figure(3)
% plot(t,OldX)
% title('Simulation of Brownian increments')
% ylabel('Brownian increments, OldX')
% xlabel('time, seconds')

% a=[t;Y];
% 
% fid=fopen('Gaussian.dat','w');
% fprintf(fid,'%f %f\n',a);
% fclose(fid);