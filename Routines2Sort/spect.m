function [Spec,freq,dSpec] = spect(x,dt,npts,overl,displ,wind)

% [Spec,freq,dSpec] = spect(x,dt,npts,overl,displ,wind)
%
% SPECT calculates and displays the power spectrum of x
% The periodograms are calculated for packets of npts data points
% that are averaged. If x is a matrix, each column is considered
% as an array. If x is complex, the spectrum is extended up to the
% sampling frequency
%	x       input signals (data sets are columns)
%	dt      sampling period, default is 1
%	npts    number of data points in each periodogram; default (0)
%       	takes largest power of 2
%	overl	overlapping factor in [%] between adjacent windows
%           default is 50 %
%	displ   spectra are displayed if displ>0 (default is with display)
%	wind	window type : 0 for Hanning (default), 1 for no window,
%           2 for Welch
%
%	Spec	power spectral density (data sets are columns)
%	freq	corresponding frequency array (column)
%	dSpec	standard deviation of Spec (same size as Spec)
%
%					ThDdW 5/95,  adapted from psd.m  



if nargin<6,	wind = 0;	end
if nargin<5,	displ = 1;	end
if nargin<4,	overl = 50;	end
if nargin<3,	npts = 0;	end
if nargin<2,	dt = 1;		end
if overl<=0,	overl = 0;	end

confint = 95;               % confidence interval in [%]

[n,ncol] = size(x);
if n==1, x = x(:); n = ncol; ncol = 1; end
isimag = any(any(imag(x)~=0));		% check whether x is imaginary

if npts<=0				% default npts is largest power of 2
	npts = 2^fix(log(n)/log(2));
elseif npts>n           % padd x with zeros to increase array
	dn = (npts-n)/2;
	wind = hanning(n);              % use a Hanning filter first
	scale = (n/npts)*norm(wind)^2;
	wind = wind / scale;
	x = [zeros(ceil(dn),ncol); x.*(wind*ones(1,ncol)); ...
             zeros(floor(dn),ncol)];
	n = npts;
end

if wind==2,
	window = welch(npts);
elseif wind==1,
	window = ones(npts,1);
else
	window = hanning(npts);
end


window = window(:);
nwind = length(window);                     % length of window
noverlap = fix(overl/100*nwind);            % overlapping in # of pts
neff = fix((n-noverlap)/(nwind-noverlap)); 	% Number of windows


index = 1:nwind;
KMU = neff*norm(window)^2; 	% Normalizing scale factor
 
Spec = zeros(npts,ncol); 
Spec2 = zeros(npts,ncol);
tr = [(1:nwind)'/nwind  ones(nwind,1)];
for i=1:neff % for each sequence, detrend and then calculate the fft
	xw = (window*ones(1,ncol)).*(x(index,:) - tr*(tr\x(index,:)));	% detrend x
	index = index + (nwind - noverlap);
	Xx = abs(fft(xw,npts)).^2;
	Spec = Spec + Xx;
	Spec2 = Spec2 + abs(Xx).^2;
end
 

% Select first half if the input is real

if isimag				% x is complex
	select = (1:npts)';
else					% if x is not complex
	if rem(npts,2)			% if npts is odd
         select = (1:(npts+1)/2)';
	else
         select = (1:npts/2+1)';
    end
	Spec = Spec(select,:);
	Spec2 = Spec2(select,:);
end
Spec = Spec*(1/KMU);
freq = (select-1)/dt/npts;
 

% find confidence interval 

if nargout>2,					% compute confidence intervals
	dSpec = zeros(size(Spec));
	if neff > 1
		c = (neff.*Spec2-abs(Spec).^2)./(neff-1);
       	c = max(c,zeros(size(Spec)));
        dSpec = sqrt(c);
	end
	ff = sqrt(2)*erfinv(confint/100);	% Equal-tails.
	dSpec = (ff.*dSpec)*(1/KMU);
 
% 	compute variance using Welch's equation (6.112 in Rabiner and Gold)
% 	(assumes underlying Gaussian process)
% 	D = nwind - noverlap;
% 	wcorr = xcorr(window,window);
% 	rho = (wcorr((nwind+D):D:length(wcorr))/wcorr(nwind)).^2
% 	var = abs(Spec).^2/neff*( 1 +  (neff-(1:length(rho)))*rho(:)*(2/neff) );
	if noverlap>0
		disp('** confidence intervals inaccurate for OVERLAP>0 **')
	end
end



if displ,
	clf,
	if neff==1,
		disp(['** using 1 sequence of  ',int2str(npts),...
                '  out of  ',int2str(n),'  **'])
	else
		disp(['** using ',int2str(neff),' sequences of  ',int2str(npts),...
                '  out of  ',int2str(n),'  **'])
	end
	nf = length(freq);
	if isimag
		wpos = (2:nf/2);	fpos = freq(wpos);
		wneg = (nf/2+1:nf);	fneg = flipud(freq(2:length(wneg)+1));
		subplot(211)
		loglog(fpos,Spec(wpos,:)),	grid
		xlabel('Frequency > 0'); ylabel('Power Spectral Density')
		subplot(212)
		loglog(fneg,Spec(wneg,:)),	grid
		xlabel('Frequency < 0'); ylabel('Power Spectral Density')
	else
		loglog(freq(2:nf),Spec(2:nf,:))
		grid
		xlabel('Frequency'), ylabel('Power Spectral Density');
	end
end