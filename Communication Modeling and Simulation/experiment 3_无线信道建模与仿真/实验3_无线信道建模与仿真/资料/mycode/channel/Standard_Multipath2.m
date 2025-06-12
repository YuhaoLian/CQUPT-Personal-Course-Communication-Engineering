function  [original_signal1]=multipath( original_signal,Fc,V,Fs)
data_length=length(original_signal);
C = 3e8;                      %light speed
%Fc = 900e6;                  %carrier frequency
%V=70;                             %km/h
%Tc = 1/(19.2e3);     %=Ts            %Chip duration according to digital bandwidth
Del_Idx = [1 2 3 4 6];        %Model B with rate = 1
Del_Gain = [0.8084 0.462 0.253 0.259 0.0447];
Sample_Rate = 1;
Ts=1/Fs;

Max_Del= max(Del_Idx);        %Max time delay
Path_Num= length(Del_Idx);    %Time invariant multipath
Max_Num = 4;                  %number of Rake fi

N0=40; 
N=4*N0;
M=Path_Num;
Ank=2*pi/(N*M)*((0:N0-1)'*ones(1,M)*M+ones(N0,1)*(0:M-1)+1/4);
Wm=2*pi*V*Fc/C*1000/3600;
phi=2*pi*[0.8 0.28];
%temp0=trans_signal11;
for kk = 1:data_length
   CH_Data(:,kk)=sqrt(1/N0)*sum(cos(Wm*cos(Ank)*(kk+1e8)*Ts/Sample_Rate+phi(1))+j*sin(Wm*sin(Ank)*(kk+1e8)*Ts/Sample_Rate+phi(2)))';
end
% hist(abs(CH_Data(1,:)),100);
for i = 1 : length(Del_Idx)
	if Del_Idx(1,i) > 1
	Temp_Sig(i,:) = [ zeros(1,Del_Idx(1,i)-1) original_signal(1,[1:end-Del_Idx(1,i)+1])]; 
	% Temp_Sig0(i,:) = [ zeros(1,Del_Idx(1,i)-1) temp0(1,[1:end-Del_Idx(1,i)+1])]; 
	else     
	Temp_Sig(i,:) = original_signal;
	%Temp_Sig0(i,:) = temp0;
	end
end  
original_signal1= Del_Gain*(Temp_Sig.*CH_Data); %pass through the channel
% hist(abs(original_signal1),100);
%seek=Del_Gain*(Temp_Sig0.*CH_Data);
%fre=fft(original_signal )
%plot(abs(original_signal));
return  