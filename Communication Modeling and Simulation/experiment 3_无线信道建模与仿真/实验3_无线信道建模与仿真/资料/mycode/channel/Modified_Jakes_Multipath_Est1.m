%From Yingbo Li, Modified Jakes` Model for Simulating Multiple Uncorrelated Fading Waveforms Model1, 
%2000 IEEE Communications Letters.

function  [fadingcoeff]=multipath(infrom_length,Fc,V,Fs,path_num)

C=3e8;                                %light speed Km/s
Fm=(V*Fc*1000)/(C*3600);              %Doppler shift
Ts = 1/Fs;                            %Chip duration according to digital bandwidth Fs
Wm=2*pi*Fm;
Ts_adjustment=7.5+Ts*infrom_length*rand(1);
t=Ts+Ts_adjustment:Ts:Ts*infrom_length+Ts_adjustment;   % Time array


N0=40;
N=4*N0+1;
wn=Wm*cos(2*pi*((1:N0)-0.5*ones(1,length(N0)))/N);
phi=pi*(1:N0)/N0;

Xc=zeros(path_num,length(t));
Xs=zeros(path_num,length(t));
fadingcoeff=zeros(path_num,infrom_length);

for pp=1:path_num
    
    for n=1:N0
        Xc(pp,:)=2*cos(phi(n)*pp)*cos(wn(n)*t+pi*pp/2)+Xc(pp,:);
        Xs(pp,:)=2*sin(phi(n)*pp)*cos(wn(n)*t+pi*pp/2)+Xs(pp,:);
    end
    fadingcoeff(pp,:)=(Xc(pp,:)+sqrt(-1)*Xs(pp,:))/sqrt(N0);
    
end