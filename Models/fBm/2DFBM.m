
clear all; close all;

randn('state', sum(100*clock));

X = normrnd(0,1,10000,1000);

% try to make the dimensions even as putting odd ones will require extra
% work for me.

Y=fft2(X);

[f1,f2] = freqspace(size(Y));

new1=[f1(end/2+1:end),f1(1:end/2)];
new2=[f2(end/2+1:end),f2(1:end/2)];

[x1,y1] = meshgrid(new1,new2);

Z=1./sqrt(x1.^2 + y1.^2);
%%
[xinf,yinf]=find(Z==inf);

Z(xinf,yinf)=0*Z(xinf,yinf+1);

% surf(x1,y1,Z)

%%

s=Z.*Y;

XBM=ifft2(s);
%%
surfc(XBM(1:100,1:100)); figure(gcf)

% We need to add a background DC field to this field now. Otherwise we
% will not get a spin tone effect.