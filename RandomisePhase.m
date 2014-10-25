function [ data ] = RandomisePhase( data )
%UNTITLED Summary of this function goes here
%phase randomisation of the signals [get rid of]

[~,columns]=size(data);
a=fft(data(:,2:columns));
theta=2*pi.*rand(length(a),columns-1)-pi;
b=ifft(abs(a).*exp(1i.*(theta)),'symmetric');
data=[data(:,1),b];
clear('b','theta','a','columns');

end

