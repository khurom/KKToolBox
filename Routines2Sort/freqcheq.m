%% M=1
clear all;
a=[]; b=1; sf=1;

for i=1:1:11
    b=b/2;
    a=[a b];
end

plot(a.*sf,zeros(length(a),1),'.b');
hold all;

%% M=2
clear all;
a=[]; b=1; sf=1;

for i=1:1:11
    
    if i<2
        b=b/2;
        a=[a b];
    end
    
    if i>1 
        a=[a (a(end)+(a(end)/2)) (a(end)-(a(end)/2))];
    end
    
end

x=zeros(length(a),1);
x=0.2;

plot(a.*sf,x,'.r'); hold all;

%% M=3
clear all;
a=[]; b=1; sf=1;

for i=1:1:11
    
    if i<2
        b=b/2;
        a=[a b];
    end
    
    if i>1 && i<3
        a=[a (a(end)+(a(end)/2)) (a(end)-(a(end)/2))];
    end
    
    if i>2
        a=[a (a(end-1)+(a(end)/2)) (a(end-1)-(a(end)/2)) (a(end)+(a(end)/2)) (a(end)-(a(end)/2))];
    end
    
end

x=zeros(length(a),1);
x=0.4;

plot(a*sf,x,'.g'); hold all;

%%

clear all;
sf=1;
s=1.2.^-[0:1:40]; 
x=zeros(length(s),1);
x=-0.2;
semilogx(s.*sf,x,'.k'); hold all;
