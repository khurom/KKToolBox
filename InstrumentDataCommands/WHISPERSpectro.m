
function WHISPERSpectro(filename)

p=cdfread(filename);
imagesc(1:1:6213,p(6,1).Data(1:480)',log10(p(5,1).Data(:,1:480)'))
set(gca,'YDir','normal')
diary off
