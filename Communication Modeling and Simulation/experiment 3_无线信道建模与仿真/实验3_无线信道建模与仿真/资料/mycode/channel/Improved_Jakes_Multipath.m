function  [Output_signal]=multipath(original_signal,Fc,V,Fs)

PowerDelayProfile = [0,-1,-9,-10,-15,-20];
NonZeroPaths = [1,3,5,7,10,14];
Path_Num=length(NonZeroPaths);
infrom_length=length(original_signal);
Ts=1/Fs;                    %Chip duration according to digital bandwidth Fs
t=Ts:Ts:Ts*length(original_signal);
C=3e8;                      %light speed Km/s
fm=V*Fc/C*1000/3600;
wm=2*pi*fm;


M=40;
N=4*M;
theta=pi*(rand(1,1)-0.5);
phi_Zc=pi*(rand(1,1)-0.5);
phi_Zs=pi*(rand(1,1)-0.5);
alpha=(2*pi*(1:M)-(pi+theta)*ones(1,length(M)))/N;

Zc=zeros(Path_Num,length(t));
Zs=zeros(Path_Num,length(t));

for pp=1:Path_Num
    
    for n=1:M
        Zc(pp,:)=sqrt(2/M)*cos(wm*t*cos(alpha(n))+phi_Zc*ones(1,length(t)))+Zc(pp,:);
        Zs(pp,:)=sqrt(2/M)*sin(wm*t*sin(alpha(n))+phi_Zs*ones(1,length(t)))+Zs(pp,:);
    end
    fadingcoeff(pp,:)=(Zc(pp,:)+sqrt(-1)*Zs(pp,:))/sqrt(M);

    
	if NonZeroPaths(1,pp) > 1
	   original_signal_temp(pp,:) = [zeros(1,NonZeroPaths(1,pp)-1) original_signal(1,[1:end-NonZeroPaths(1,pp)+1])];
	else
	   original_signal_temp(pp,:) = original_signal;
	end
    
end

MeanAmp = 10.^(PowerDelayProfile/10); %Mean amplitude values for different paths in decimal format:
MeanAmp = sqrt( MeanAmp/sum(MeanAmp) );	%Nomorlized the power to 1.
Output_signal=MeanAmp*(original_signal_temp.*fadingcoeff);

% hist(abs(fadingcoeff(1,:))*length(fadingcoeff(1,:)),100);
% ylabel('Number of Hits');
% xlabel('Magnitude');
