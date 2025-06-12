% Simulation of my arithmetic
close all;
clear;
clc;
tic;
fc = 2e9;			           %Carrier frequency
Insert_samples=3;
fs = 20e+6; 			       %System bandwith is 20MHz
V=100;
SubCarrNum = 16;	           %Number of carriers
FrameNum=10;
Info_length=(SubCarrNum)*Insert_samples;
frequency_offset=0.25;
time_initialization=0;
PilotSymbol = (1+j)/sqrt(2);
SymbolNum=10;                          %Symbol Number

PN_code=PN_Signal(2*SubCarrNum);
PN_code=PN_code(1,(SubCarrNum+1):1:end);
% PN_code=AGWN_Generator(Info_length);

%Insert sampls
PN_code_temp=zeros(1,length(PN_code)*Insert_samples);
for i=1:length(PN_code)
    for ii=1:Insert_samples
        PN_code_temp(1,(i-1)*Insert_samples+ii)=PN_code(1,i);
    end
end



for SNR=0:1:10

    error_MSE_time=0;
	error_MSE_angle=0;
    
	for nframe=1:FrameNum
		snr=10^(2*SNR/10);
		E=1;                % 符号能量假设为单位能量.
		Noise_Var = E/(2*snr);
       	for m = 1:1:SymbolNum

            OfdmSymbol_Pilot = PilotSymbol*ones(1,Info_length);
            OfdmSymbol_Pilot_I=real(OfdmSymbol_Pilot).*real(PN_code_temp);   
            OfdmSymbol_Pilot_Q=imag(OfdmSymbol_Pilot).*imag(PN_code_temp);    
            OfdmSymbol_Pilot=OfdmSymbol_Pilot_I+sqrt(-1)*OfdmSymbol_Pilot_Q;     
            
            if m==1
               ChannelInput=OfdmSymbol_Pilot;
            else
               ChannelInput=[OfdmSymbol_Pilot ChannelInput];
            end
 		
		end
        
        %Pass through channel
		ChannelOutput=Frequency_offset_channel(ChannelInput,frequency_offset);
        ChannelOutput=AGWN_channel(ChannelOutput,Noise_Var);
        ChannelOutput=Standard_Multipath1(ChannelOutput,fc,V,fs);  
        ChannelOutput=[zeros(1,time_initialization),ChannelOutput];
%         %Remove sampls
%         Step_cz=Insert_samples;
%         for i=1:Info_length*SymbolNum
%             ChannelOutput_temp(1,i)=ChannelOutput(1,Insert_samples*(i-1)+Step_cz);
%         end
%         ChannelOutput=ChannelOutput_temp;
        
		ChannelOutput2=[zeros(1,Info_length) ChannelOutput zeros(1,Info_length)];        
		for ii=1:(SymbolNum+1)*Info_length;
            Gamam_P_temp=0;
			for jj=1:length(PN_code_temp)
                Gamam_P_temp=Gamam_P_temp+conj(PN_code_temp(1,jj))*ChannelOutput2(1,jj+ii-1);
			end
            Gamam_P(1,ii)=Gamam_P_temp;
		end
		
		Gamam_P=Gamam_P(1,2:1:end);
		Gamam_PP=Gamam_P;
		Gamam_P1=[Gamam_P zeros(1,length(PN_code_temp))];
		Gamam_P2=[zeros(1,length(PN_code_temp)) Gamam_P];
		Gamam_P=Gamam_P1.*conj(Gamam_P2);
		Gamam_P=Gamam_P(1,1:1:(end-Info_length+1));
		argmax=abs(Gamam_P);
		       

		argMax=max(argmax(1,(Info_length*SymbolNum-10):1:(Info_length*SymbolNum+10)));
        for jj=(Info_length*SymbolNum-10):1:(Info_length*SymbolNum+10)
            if argmax(1,jj)==argMax 
               angle_offset=angle(Gamam_P(1,jj))*SymbolNum/(2*pi)-frequency_offset;
               time_offset=jj-Info_length*SymbolNum-time_initialization;
               error_MSE_time=abs(time_offset)^2+error_MSE_time;
               error_MSE_angle=abs(angle_offset)^2+error_MSE_angle;
            end
		end
    end
    MSE_time(1,SNR+1)=sqrt(error_MSE_time/FrameNum);
    MSE_angle(1,SNR+1)=sqrt(error_MSE_angle/FrameNum);
 
end

%Output figures

figure(1)
semilogy(2*(0:1:10),MSE_time,'x');
grid on;
hold on;

figure(2)
semilogy(2*(0:1:10),MSE_angle,'x');
grid on;
hold on;

figure(3)
correlation_flag=length(PN_code_temp); 
[R_xx,lags] = xcorr(PN_code_temp,correlation_flag,'coeff');
subplot(211);
plot(1:length(PN_code_temp),real(PN_code_temp));
grid on;

subplot(212);
plot(lags,real(R_xx))
grid on;

figure(4)
subplot(311);
plot(1:length(Gamam_PP),abs(Gamam_PP));
grid on;

subplot(312);
plot(1:length(argmax)-Info_length,argmax(1,Info_length+1:end));
grid on;

subplot(313);
for m = 2:SymbolNum
    plot((-Info_length/2)+1:1:Info_length/2,argmax(1,m*Info_length-(Info_length/2)+1:(m+1)*Info_length-(Info_length/2)));  
    hold on
end
grid on;

TimeElapse = toc/60;	% Simulation minutes
fprintf('****** The total time used is : %3.4f minutes ********\n', TimeElapse);