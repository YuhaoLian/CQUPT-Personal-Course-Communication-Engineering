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


N0=40;
N=4*N0+2;
wn=2*pi*fm*cos(2*pi*(1:N0)/N);
phi=pi*(1:N0)/N0;

Xc=zeros(Path_Num,length(t));
Xs=zeros(Path_Num,length(t));

for pp=1:Path_Num
    
    alpha=pi*pp;
    Xc_initialization=2*sqrt((cos(alpha))^2)*cos(wm*t+pi*pp/2);
    Xs_initialization=2*sin(alpha)*cos(wm*t+pi*pp/2);
    
    for n=1:N0
        if n==1
           Xc(pp,:)=2*cos(phi(n)*pp)*cos(wn(n)*t+pi*pp/2)+Xc_initialization;
           Xs(pp,:)=2*sin(phi(n)*pp)*cos(wn(n)*t+pi*pp/2)+Xs_initialization;
        else
           Xc(pp,:)=2*cos(phi(n)*pp)*cos(wn(n)*t+pi*pp/2)+Xc(pp,:);
           Xs(pp,:)=2*sin(phi(n)*pp)*cos(wn(n)*t+pi*pp/2)+Xs(pp,:);
        end
    end
    fadingcoeff(pp,:)=(Xc(pp,:)+sqrt(-1)*Xs(pp,:))/sqrt(N0);
    
	if NonZeroPaths(1,pp) > 1
	   original_signal_temp(pp,:) = [zeros(1,NonZeroPaths(1,pp)-1) original_signal(1,[1:end-NonZeroPaths(1,pp)+1])];
	else
	   original_signal_temp(pp,:) = original_signal;
	end
    
end

MeanAmp = 10.^(PowerDelayProfile/10); %Mean amplitude values for different paths in decimal format:
MeanAmp = sqrt( MeanAmp/sum(MeanAmp) );	%Nomorlized the power to 1.
Output_signal=MeanAmp*(original_signal_temp.*fadingcoeff);

% hist(abs(Output_signal),100);
% ylabel('Number of Hits');
% xlabel('Magnitude');
