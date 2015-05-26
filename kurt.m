% Function to calculate Kurtosis as a function of scale

function [scales,kurt]=kurt(data,scales,displ)


% displ takes the values '1' for plots, or '0' for no plots.

%% Setting up font sizes for plots
set(0,'defaultaxesfontsize',16);
set(0,'defaulttextfontsize',16);

%% Structure functions

kurt=ones(length(scales),1);

for n=1:1:length(scales),
    
        kurt(n)=kurtosis(data{1,n});

end

%% plots

if displ
    
    clf
    
    figure(1);
    semilogx(scales,kurt,'.r');
    xlabel('\tau [secs])',ylabel('Kurtosis');

end