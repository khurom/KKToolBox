
s=BFGM(:,4);

maxsize=nextpow2(numel(s)) - 1;

s=s(1:(2^maxsize));

leve=3;    % There is a maximum level that you can go up to here and I think you do not want
                   % to go up to this maximum level
%%
tic

[swa,swd] = swt(s,leve,'db4');

swa=swa'; swd=swd';
% 
% figure(1)
% subplot(4,2,1)
% plot(s);
% 
% subplot(4,2,2)
% plot(s);
% 
% subplot(4,2,3)
% plot(swa(:,1));
% 
% subplot(4,2,4)
% plot(swd(:,1));
% 
% subplot(4,2,5)
% plot(swa(:,2));
% 
% subplot(4,2,6)
% plot(swd(:,2));
% 
% subplot(4,2,7)
plot(swa(:,3)); hold all;
% 
% subplot(4,2,8)
% plot(swd(:,3));

toc

%%
% subplot(4,2,7)
plot(s,'.b'); hold all;

tic

p= iswt(swa,zeros(size(swd)),'db2');

toc

% for i=1:1:(leve-1)
%     p= iswt(p,zeros(size(p)),'db4');
% end

plot(p,'-r');

%%
[cA,cD] = dwt(s,'db4','mode','per');
%%
leve=13;
[dwtwc,l] = wavedec(s,leve,'db20');
% p=dwtwc(1:l(1));
% plot(p);
%%
q=dyadup(dwtwc(1:l(1)),1);
q=dyadup(q,1);
q=dyadup(q,1);
%%
plot(q(110:500).*(1/(2*sqrt(2))),'Marker','.','LineStyle','none','Color',[0 0.498 0]); hold all;
hold on; plot(s(1:500),'.b'); % This 100 value is arbitray; just me guessing
%%
q=upsample(dwtwc(1:l(1)),2);
q=upsample(q,2);
q=upsample(q,2);
plot(q(16:end).*(1/(2*sqrt(2))),'Marker','.','LineStyle','none','Color',[0 0.498 0]); hold all;
%%
q=downsample(q,2);
q=downsample(q,2);
q=downsample(q,2);
plot(dwtwc(1:l(1)),'-b'); hold all; 
plot(q,'Marker','.','LineStyle','none','Color',[0 0.498 0]); hold all;
% plot(q(36:end).*(1/(2*sqrt(2))),'Marker','.','LineStyle','none','Color',[
% 0 0.498 0]); hold all;
%% plots the last approximation over the original signal

% figure(2)
% subplot(5,1,1)
% plot(s);

% subplot(4,2,7)
% q=upsample(dwtwc(1:l(1)),8);
q=idwt(dwtwc(1:l(1)),zeros(numel(dwtwc(1:l(1))),1),'db20');
for i=1:1:(leve-1)
    q=idwt(q,zeros(numel(q),1),'db20');
end

plot(s,'.b'); hold all; plot(q,'.r');

% plot(dwtwc(1:l(1)),'.r');

% sum=l(1);
% 
% subplot(5,1,5)
% plot(dwtwc(sum+1:sum+l(2)));
% 
% sum=sum+l(2);
% 
% subplot(5,1,4)
% plot(dwtwc(sum+1:sum+l(3)));
% 
% sum=sum+l(3);
% 
% subplot(5,1,3)
% plot(dwtwc(sum+1:sum+l(4)));
% 
% sum=sum+l(4)

