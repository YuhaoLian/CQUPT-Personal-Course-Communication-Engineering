%From Yahong R. Zheng, Improved Models for the Generation of Multiple Uncorrelated Rayleigh Fading Waveforms, 
%2002 IEEE Communications Letters.

function  [fadingcoeff]=multipath(infrom_length,Fc,V,Fs,path_num)

C=3e8;                                %light speed Km/s
Fm=(V*Fc*1000)/(C*3600);              %Doppler shift
Ts = 1/Fs;                            %Chip duration according to digital bandwidth Fs
Wm=2*pi*Fm;

M=40;
N=4*M;
theta=pi*(rand(1,1)-0.5);
phi_Zc=pi*(rand(1,1)-0.5);
phi_Zs=pi*(rand(1,1)-0.5);
alpha=(2*pi*(1:M)-(pi+theta)*ones(1,length(M)))/N;

fadingcoeff=zeros(path_num,infrom_length);

for ii=1:path_num
    
    Ts_adjustment=7.5+Ts*infrom_length*rand(1);
    t=Ts+Ts_adjustment:Ts:Ts*infrom_length+Ts_adjustment;   % Time array
    Zc=zeros(1,length(t));
    Zs=zeros(1,length(t));
    
	for n=1:M
        Zc(1,:)=sqrt(2/M)*cos(Wm*t*cos(alpha(n))+phi_Zc*ones(1,length(t)))+Zc(1,:);
        Zs(1,:)=sqrt(2/M)*cos(Wm*t*sin(alpha(n))+phi_Zs*ones(1,length(t)))+Zs(1,:);
	end
    
	fadingcoeff(ii,:)=(Zc(1,:)+sqrt(-1)*Zs(1,:))/sqrt(M);
    
end