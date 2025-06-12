% Simulation of《Symbol time offset estimation in choerent OFDM systems》
close all;
clear;
clc;
tic;
fc = 2.4e9;			%Carrier frequency
fs = 5e+6; 			%System bandwith is 5MHz
V=5;
SubCarrNum = 256;	           %Number of carriers
GuardNum = 8;       %Number of cyclic prefix
FrameNum=10;
Info_length=SubCarrNum+GuardNum;
frequency_offset=0.1;
PilotSymbol = (1+j)/sqrt(2);
M=2;

PN_code=PN_Signal(Info_length);
% PN_code=AGWN_Generator(Info_length);

kk=0;
for SNR=0:2:20
    kk=kk+1;
    error_MSE_time=0;
	error_MSE_angle=0;
	for nframe=1:FrameNum
		snr=10^(SNR/10);
		E=1;                % 符号能量假设为单位能量.
		Noise_Var = E/(2*snr);
       	for m = 1:1:M
		
            OfdmSymbol_Pilot = PilotSymbol*ones(1,Info_length);
            OfdmSymbol_Pilot_I=real(OfdmSymbol_Pilot).*real(PN_code);   
            OfdmSymbol_Pilot_Q=imag(OfdmSymbol_Pilot).*imag(PN_code);    
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
        
		ChannelOutput2=[zeros(1,Info_length) ChannelOutput zeros(1,Info_length)];        
		for ii=1:(M+1)*Info_length;
            Gamam_P_temp=0;
			for jj=1:length(PN_code)
                Gamam_P_temp=Gamam_P_temp+conj(PN_code(1,jj))*ChannelOutput2(1,jj+ii-1);
			end
            Gamam_P(1,ii)=Gamam_P_temp;
		end
		
		Gamam_P=Gamam_P(1,2:1:end);
		Gamam_PP=Gamam_P;
		Gamam_P1=[Gamam_P zeros(1,length(PN_code))];
		Gamam_P2=[zeros(1,length(PN_code)) Gamam_P];
		Gamam_P=Gamam_P1.*conj(Gamam_P2);
		Gamam_P=Gamam_P(1,1:1:(end-Info_length+1));
		argmax=abs(Gamam_P);
		       
		argMax=max(argmax(1,Info_length*M-Info_length/4:1:Info_length*M+Info_length/4));
        for jj=Info_length*M-Info_length/4:1:Info_length*M+Info_length/4
            if argmax(1,jj)==argMax 
               angle_offset=angle(Gamam_P(1,jj))*M/(2*pi)-frequency_offset
               time_offset=jj-Info_length*M;
               if time_offset==0
               else
                  error_MSE_time=abs(time_offset)^2+error_MSE_time;
               end
               if angle_offset==0
               else
                  error_MSE_angle=abs(angle_offset)^2+error_MSE_angle;
               end 
            end
		end
    end
    MSE_time(1,kk)=sqrt(error_MSE_time/FrameNum);
    MSE_angle(1,kk)=sqrt(error_MSE_angle/FrameNum);
  
end

%Output figures

figure(1)
semilogy(1:kk,MSE_time,'x');
grid on;
hold on;

figure(2)
semilogy(1:kk,MSE_angle,'x');
grid on;
hold on;

figure(3)
correlation_flag=length(PN_code); 
[R_xx,lags] = xcorr(PN_code,correlation_flag,'coeff');
subplot(211);
plot(1:length(PN_code),real(PN_code));
grid on;

subplot(212);
plot(lags,real(R_xx))
grid on;

figure(4)
subplot(211);
plot(1:length(Gamam_PP),Gamam_PP);
grid on;

subplot(212);
for m = 1:M
    plot(0:(Info_length-1),Gamam_PP(1,(m-1)*Info_length+1:m*Info_length));  
    hold on
end
grid on;

figure(5)
subplot(211);
plot(1:length(argmax),argmax);
grid on;

subplot(212);
for m = 1:M
    plot(0:(Info_length-1),argmax(1,(m-1)*Info_length+1:m*Info_length));  
    hold on
end
grid on;

TimeElapse = toc/60;	% Simulation minutes
fprintf('****** The total time used is : %3.4f minutes ********\n', TimeElapse);