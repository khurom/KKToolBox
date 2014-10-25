function [sj,fw,a]=WavelettPSDVar(data,ts)

% ts; N(number of points for fft);

fs = 1./ts;   
fN = 0.5*fs;        % Nyquist frequency

% V=[];
% 
% for i=-3.0:0.1:1.0,
%     V=[V 10^(i)];
% end

wname='morl';
A0 = [16 32 64 128 256 512 1024 2048 4096 8192];
% A0=V;
F0 = scal2frq(A0,wname,ts);
j0 = [F0 < fN];
sj = A0(j0);   % scales
fw = F0(j0);

[Xw] = cwt(data,sj,wname);

a=mean(Xw,2);

plot(log10(fw'),log10(a),'.b');

% Try wscalogram