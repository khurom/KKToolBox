%**********************************************************************************
function [SumPar]=SummaryPlots(Bvec,Particles,Electrons,trendscale)
%**********************************************************************************
% 
% File name: SummaryPlots.m
% Author: Khurom H. Kiyani
% File created: 27th September 2013
% Last updated: 3rd November 2013
% Updated by: Khurom H. Kiyani
% 
% Description:
% ------------
% This m-file takes in the various field and plasma parameters for the
% interval of interest and produces summary plots of plasma parameters.
% 
% Input variables/signals:
% 
% 'Bvec':           Magnetic field *vector* [nT]
% 'Particles':      Array of data for the plasma parameters: 
%       
%           'Dens':           Density [cm^-3]
%           'Vvec':           Velocity field *vector* [km s^-1]
%           'T':              Total temperature [MK]
%           'Tpar':           Temperature parallel to the magnetic field [MK]
%           'Tperp':          Temperature perpendicular to the magnetic field [MK]
% 
% 'sc_flag':      Put a flag for spacecraft [cluster, wind, ace]
% 'trendscale':     The scale in seconds at which you want to perform a 
%                   low-pass filter to do a moving average and plot a trend-line
% 
% Output variables/signals:
% 
% 'SumPar':         Summary cell array of mean values and percentage errors
% 
%**********************************************************************************
% 
% Additional Notes:
% -----------------
% 0./ Don't worry if the data has it's time in datenum format. As long as 
%     the real sampling is in units of seconds everything is ok. 
% 1./ Add Whisper spectrogram plot.
% 2./ Add ion spectrogram plot.
% 3./ Add ion temperature anisotropy Vs. ion plasma beta parallel plot with
%     theoretical marginal instability thresholds plotted.
%**********************************************************************************

% Setting sampling frequencies [assumes uniform sampling]

SampRes=NaN(3,1);

SampRes(1)=nanmean(1./(24*60*60.*(diff(Bvec(:,1)))))
SampRes(2)=nanmean(1./(24*60*60.*(diff(Particles(:,1)))))
SampRes(3)=nanmean(1./(24*60*60.*(diff(Electrons(:,1)))))

Instr=cell(3,1);
Instr{1}='Bvec';
Instr{2}='Particles';
Instr{3}='Electrons';

% Choose the instrument with the smallest sampling frequency to be the one
% which everything else gets interpolated on.

InterpRes=min(SampRes);

ind=(SampRes==min(SampRes));
time=eval(strcat(Instr{ind},'(:,1)'));

% switch mode_flag
%     case 'BM'
%         SampRes(1)=67; SampRes(2)=0.241;
%     case 'NM'
%         SampRes(1)=22.5; SampRes(2)=0.241;
% end

otherind=find(ind<1);

for m=1:2
    
    tempvar=eval(Instr{otherind(m)});
    tempvar=irf_filt(tempvar,0,InterpRes,SampRes(otherind(m)),3);
    tempvar=[time,interp1(tempvar(:,1),tempvar(:,2:end),time)];
    
    if otherind(m)>2
        Electrons=tempvar;
    elseif otherind(m)>1
        Particles=tempvar;
    else
        Bvec=tempvar;
    end
    
    clear tempvar;
    
end

% Calculate summary parameters for plotting Ion Plasma Beta

Bmag=sqrt(Bvec(:,2).*Bvec(:,2) + Bvec(:,3).*Bvec(:,3) + Bvec(:,4).*Bvec(:,4));

Beta=(Particles(:,6).*Particles(:,2).*(8.*pi.*1.3807))./(Bmag.^2);

% Set up some figure parameters

scrsz = get(0,'ScreenSize');
figure1 = figure('Position',[scrsz(3)-(scrsz(3)/2) ...
    scrsz(4)-(scrsz(4)/1.1) scrsz(3)/3 scrsz(4)/1.5]);

height=0.12; width=0.775; leftpos=0.138; givespace=0.005;

% Plot magnitude of B

axes1 = axes('Parent',figure1, ...
    'Position',[leftpos 0.079+6*(height+givespace) width height]);
box(axes1,'on'); hold(axes1,'all');

plot(time,Bmag,'Parent',axes1); 
% trend line
out=irf_filt(ConstNan(time,Bmag),0,1/trendscale,InterpRes,3);
out(isnan(Bmag))=NaN;
plot(time,out,'r','Linewidth',2);
clear out;
%
ylabel('|B| [nT]','FontSize',12); datetick('x'); 
set(axes1,'XTickLabel',''); axis tight; 

title('Plasma Summary Parameters','FontWeight','bold','FontSize',16);

% Plot ion plasma beta

axes2 = axes('Parent',figure1, ...
    'Position',[leftpos 0.079+5*(height+givespace) width height]);
box(axes2,'on'); hold(axes2,'all');

plot(time,Beta,'Parent',axes2,'Color','k'); 
% trend line
out=irf_filt(ConstNan(time,Beta),0,1/trendscale,InterpRes,3);
out(isnan(Beta))=NaN;
plot(time,out,'r','Linewidth',2);
clear out;
%
ylabel('plasma \beta_{i}','FontSize',12); datetick('x'); 
set(axes2,'XTickLabel',''); axis tight;

% Plot solar wind speed

speed=sqrt(Particles(:,3).*Particles(:,3)+Particles(:,4).*Particles(:,4) ...
    +Particles(:,5).*Particles(:,5));

axes3 = axes('Parent',figure1, ...
    'Position',[leftpos 0.079+4*(height+givespace) width height]);
box(axes3,'on'); hold(axes3,'all');

plot(time,speed,'Parent',axes3, ...
    'Color',[0.0784313753247261 0.168627455830574 0.549019634723663]); 
% trend line
out=irf_filt(ConstNan(time,speed),0,1/trendscale,InterpRes,3);
out(isnan(speed))=NaN;
plot(time,out,'r','Linewidth',2);
clear out;
%
ylabel('|V| [km s^{-1}]','FontSize',12); datetick('x'); 
set(axes3,'XTickLabel',''); axis tight;


% Plot total temperature

axes4 = axes('Parent',figure1, ...
    'Position',[leftpos 0.079+3*(height+givespace) width height]);
box(axes4,'on'); hold(axes4,'all');

plot(time,Particles(:,6),'Parent',axes4,'Color','r'); 
% trend line
out=irf_filt(ConstNan(time,Particles(:,6)),0,1/trendscale,InterpRes,3);
out(isnan(Particles(:,6)))=NaN;
plot(time,out,'b','Linewidth',2);
clear out;
%
ylabel('T [MK]','FontSize',12); datetick('x'); 
set(axes4,'XTickLabel',''); axis tight;

% Plot ion number density

axes5 = axes('Parent',figure1, ...
    'Position',[leftpos 0.079+2*(height+givespace) width height]);
box(axes5,'on'); hold(axes5,'all');

plot(time,Particles(:,2),'Parent',axes5, ...
    'Color',[0 0.498039215686275 0]);
% trend line
out=irf_filt(ConstNan(time,Particles(:,2)),0,1/trendscale,SampRes(2),3);
out(isnan(Particles(:,2)))=NaN;
plot(time,out,'r','Linewidth',2);
clear out;
%
ylabel('n [cm^{-3}]','FontSize',12); datetick('x'); 
set(axes5,'XTickLabel',''); axis tight;

% Plot ratio of electron to ion temperature

axes6 = axes('Parent',figure1, ...
    'Position',[leftpos 0.079+(height+givespace) width height]);
box(axes6,'on'); hold(axes6,'all');

ElecTemp=sqrt(Electrons(:,6).*Electrons(:,6)+Electrons(:,7).*Electrons(:,7));

plot(time,ElecTemp./Particles(:,6),'Parent',axes6, ...
    'Color',[0.870588235294118 0.490196078431373 0]);
% trend line
out=irf_filt(ConstNan(time,ElecTemp./Particles(:,6)),0,1/trendscale,SampRes(2),3);
out(isnan(Particles(:,2)))=NaN;
plot(time,out,'r','Linewidth',2);
clear out;
%
ylabel('T_{e}/T_{i}','FontSize',12); datetick('x'); 
set(axes6,'XTickLabel',''); axis tight;

% Plot ion temperature anisotropy 

axes7 = axes('Parent',figure1, ...
    'Position',[leftpos 0.079 width height]);
box(axes7,'on'); hold(axes7,'all');

plot(time,Particles(:,7)./Particles(:,8),'Parent',axes7, ...
    'Color',[0.0431372549019608 0.517647058823529 0.780392156862745]); 
% trend line
out=irf_filt(ConstNan(time,Particles(:,7)./Particles(:,8)),0,1/trendscale,SampRes(2),3);
out(isnan(Particles(:,7)./Particles(:,8)))=NaN;
plot(time,out,'r','Linewidth',2);
clear out;
%
xlabel(['Time on ',datestr(Particles(1,1),'dd-mmm-yyyy')], ...
    'FontSize',12); 
ylabel('T_{i ||}/T_{i \perp}','FontSize',12); datetick('x'); axis tight; 

% % Writing mean summary parameters and their uncertainty out to cell array
% 
% SumPar=cell(6,3);
% SumPar{1,1}='Variable'; SumPar{1,2}='Mean'; SumPar{1,3}='PercentageError';
% SumPar{2,1}='Bmagnitude'; SumPar{2,2}=nanmean(Bmag); 
%                           SumPar{2,3}=nanstd(Bmag)*100/nanmean(Bmag);
% SumPar{3,1}='PlasmaBeta'; SumPar{3,2}=nanmean(Beta); 
%                           SumPar{3,3}=nanstd(Beta)*100/nanmean(Beta);
% SumPar{4,1}='Temperature'; SumPar{4,2}=nanmean(Particles(:,6)); 
%                            SumPar{4,3}=nanstd(Particles(:,6))*100/nanmean(Particles(:,6));
% SumPar{5,1}='Density'; SumPar{5,2}=nanmean(Particles(:,2)); 
%                        SumPar{5,3}=nanstd(Particles(:,2))*100/nanmean(Particles(:,2));
% SumPar{6,1}='Temperature e/i'; SumPar{6,2}=nanmean(ElecTemp./Particles(:,6)); 
% SumPar{6,3}=nanstd(ElecTemp./Particles(:,6))*100/nanmean(ElecTemp./Particles(:,6));
% SumPar{7,1}='TemperatureAnisotropy'; SumPar{7,2}=nanmean(Particles(:,7)./Particles(:,8)); 
% SumPar{7,3}=nanstd(Particles(:,7)./Particles(:,8))*100/nanmean(Particles(:,7)./Particles(:,8));

% Writing mean summary parameters and their uncertainty out to cell array without labels

SumPar=cell(1,14); 

SumPar{1,1}=nanmean(Bmag); SumPar{1,2}=nanstd(Bmag)*100/nanmean(Bmag);
SumPar{1,3}=nanmean(Beta); SumPar{1,4}=nanstd(Beta)*100/nanmean(Beta);
SumPar{1,5}=nanmean(Particles(:,6)); SumPar{1,6}=nanstd(Particles(:,6))*100/nanmean(Particles(:,6));
SumPar{1,7}=nanmean(Particles(:,2)); SumPar{1,8}=nanstd(Particles(:,2))*100/nanmean(Particles(:,2));
SumPar{1,9}=nanmean(ElecTemp./Particles(:,6)); 
SumPar{1,10}=nanstd(ElecTemp./Particles(:,6))*100/nanmean(ElecTemp./Particles(:,6));
SumPar{1,11}=nanmean(Particles(:,7)./Particles(:,8)); 
SumPar{1,12}=nanstd(Particles(:,7)./Particles(:,8))*100/nanmean(Particles(:,7)./Particles(:,8));
SumPar{1,13}=nanmean(speed); SumPar{1,14}=nanstd(speed)*100/nanmean(speed);

end

