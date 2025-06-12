%The modified Jakes's model, for which the first and second statistics
%of each path waveforms agree well with the theoretical expectations.

function  [fadingcoeff]=multipath(infrom_length,Fc,V,Fs,path_num)

C=3e8;                                    %light speed Km/s
Fm=(V*Fc*1000)/(C*3600);                  %Doppler shift
Ts = 1/Fs;                                %Chip duration according to digital bandwidth Fs
Wm=2*pi*Fm;

N0 = 24;
N  = 4*N0;
M  = path_num;
Ank= 2*pi/(N*M)*((0:N0-1)'*ones(1,M)*M+ones(N0,1)*(0:M-1)+1/4);
Phase = 2*pi*[0.8 0.28];%2*pi*rand(1);							% Random Phase

fadingcoeff=zeros(M,infrom_length);

for kk=1:infrom_length
	fadingcoeff(:,kk)=sqrt(1/N0)*sum(cos(Wm*cos(Ank)*(kk+1e8)*Ts+Phase(1))+sqrt(-1)*sin(Wm*sin(Ank)*(kk+1e8)*Ts+Phase(2))).';
end