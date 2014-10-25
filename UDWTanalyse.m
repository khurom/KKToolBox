%**************************************************************************
function [SWT_parts]=UDWTanalyse(signal,pads,h_0,Lps,Hps,nlevel)
%**************************************************************************
% 
% File name: UDWTanalyse.m
% Author: Khurom H. Kiyani
% File created: 26th November 2013
% Last updated: 27th November 2013
% Updated by: Khurom H. Kiyani
% 
% Description:
% ------------
% This function performs the Undecimated Discrete Wavelet Transform (UDWT) 
% decompositon and reconstruction to calculate the fluctuating and background
% fields at each scale. It uses mex files from the third-party Rice wavelet
% toolbox [http://dsp.rice.edu/software/rice-wavelet-toolbox]. The UDWT is also
% known in the literature as the Stationary Wavelet Transform (SWT) (as in 
% MATLAB), the Translation Invariant Wavelet Transform (TIWT) (as in Ogden's
% book), the a'trous Wavelet (as in Mallet's book), the Nondecimated Discrete
% Wavelet Transform (also in MATLAB), or the Redundant Discrete Wavelet 
% Transform (as in the Rice Wavelet Toolbox).
% 
% Input variables:
% 
% 'signal':     The signal with the first column being the time vector.
% 'pads':       The number of points by which one will have to extend the above
%               signal on either side in order to make the length of the signal a
%               power of two.
% 'h_0':        Wavelet scaling filter.
% 'Lps':        Low pass filter phase shift for level 1.
% 'Hps':        High pass filter phase shift for level 1.
% 'nlevel':     Maximum level of the wavelet decomposition.
% 
% Output variables:
% 
% 'SWT_parts':  A cell containing, for each component of the signal, the
%               decomposition of the signal into fluctuations and
%               background fields at each scale -- decomposed using the
%               UDWT. Edge effects are also removed.
% 
% Dependencies/Auxillary files:
% -----------------------------
% This function uses the mex file based Rice Wavelet Toolbox functions: 
% mrdwt and mirdwt. For the usage of these functions see the 'help' and 
% examples in the Rice Wavelet Toolbox. 
% 
%**************************************************************************

columnnum=size(signal); columnnum=columnnum(2);

SWT_parts = cell(1,4);

SWT_parts{1,1} = signal(:,1);

for m=2:1:columnnum
    
% Gets the data size up to the next power of 2 due to UDWT restrictions
% Although periodic extension is used for the wavelet edge handling 
% automatically by the UDWT routine here, we are getting the data up to the
% next power of 2 here by extending the data sample with a constant value.
    
    y=[signal(1,m).*ones(pads/2,1); signal(:,m); signal(end,m).*ones(pads/2,1)];
    
% Decompose the signal using the UDWT
    
    [swa,swd,L] = mrdwt(y,h_0,nlevel);
    
    clear('y','L');
    
% Reconstruct all the approximations and details at all levels
    
    mzero = zeros(size(swd));
    A = mzero;
    A(:,nlevel) = mirdwt(swa,mzero,h_0,nlevel);
    
    D = mzero;
    
    for i = 1:nlevel
        
        swcfs = mzero;
        swcfs(:,i) = swd(:,i);
        D(:,i) = mirdwt(zeros(length(swa),1),swcfs,h_0,nlevel);
        
    end
    
    clear('swa','i','swcfs');
    
    for j=nlevel-1:-1:1
        
        A(:,j) = A(:,j+1) + D(:,j+1);
    
    end
    
    clear('j');

% *************************************************************************
% VERY IMPORTANT: LINEAR PHASE SHIFT CORRECTION
% *************************************************************************
% Correct for linear phase shift in wavelet coefficients at each level. No
% need to do this for the low-pass filters approximations, A, as they will 
% be reconstructed and the shift will automatically be reversed. The formula
% for the shift has been taken from Walden's paper, or was constructed by
% myself (can't exactly remember) -- but it is verified and correct.
% *************************************************************************

    for j=1:nlevel

        shiftfac=Hps*(2^(j-1));
    
        for l=1:1:j
    
            shiftfac=shiftfac + Lps*(2^(l-2))*((l-2)>=0) ;
    
        end
    
        swd(:,j)=circshift(swd(:,j),shiftfac);
        clear('shiftfac','l');
    
    end
    
    clear('j');
    
% *************************************************************************
    
    A=A';
    D=D';
    swd=swd';
    
    % Put all the files together into a cell structure
    
    SWT_parts{1,m} = cell(1,2);
    SWT_parts{1,m}{1,1} = A(:,((pads/2)+1):(end-(pads/2)));
        clear('A');    
    SWT_parts{1,m}{1,2} = swd(:,((pads/2)+1):(end-(pads/2)));
        clear('swd'); clear('D');
        
    % Replace the above two lines for 'swd' for these following two lines
    % if you desire to have the fluctuations in real physical space, 'D',
    % as with the approximations'A'.
%     SWT_parts{1,m}{1,2} = D(:,((pads/2)+1):(end-(pads/2)));
%         clear('D'); clear('swd');
    
end
    
    clear('m');
    
% Getting rid of the edge effects; to keep edges skip this section

filterlength = length(h_0);

for j=1:1:nlevel
    
%     extra = (2^(j-2))*filterlength; % easy edge removal 
%     extra = ((2^(j-1))-0.5)*filterlength; % middle-ground edge removal 
    extra = ((2^j)-1)*filterlength; % harsh edge removal 
    
    for m=2:1:4
        
        % for approximations
        SWT_parts{1,m}{1,1}(j,1:extra)=NaN;
        SWT_parts{1,m}{1,1}(j,end-extra+1:end)=NaN;
        
        % for details
        SWT_parts{1,m}{1,2}(j,1:extra)=NaN;
        SWT_parts{1,m}{1,2}(j,end-extra+1:end)=NaN;
        
    end
    
    clear('m','extra');
    
end

clear('j','filterlength');

end