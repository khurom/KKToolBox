
function one_d_iondis(filename)

% Ion 1-D Energy Flux plot commands in matlab

p=cdfread(filename);
imagesc(1:1:811,flipud((p(9,1).Data')./1000),log10(p(8,1).Data'));
% imagesc(1:1:811,flipud((p(8,1).Data')./1000),log10(p(11,1).Data'));