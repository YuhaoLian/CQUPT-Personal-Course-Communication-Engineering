% Simulation of《Symbol time offset estimation in choerent OFDM systems》
close all;
clear;
clc;
tic;
fc = 2e9;			%Carrier frequency
fs = 20e+6; 		%System bandwith is 20MHz
V=0;
SubCarrNum = 16;	%Number of carriers
GuardNum = 0;       %Number of cyclic prefix
FrameNum=10;
Info_length=SubCarrNum+GuardNum;
frequency_offset=0.1;
PilotSymbol = (1+j)/sqrt(2);
M=10;

PN_code=PN_Signal(Info_length*2);
PN_code=PN_code(1,SubCarrNum+1:end);
% PN_code=AGWN_Generator(Info_length);

kk=0;
for SNR=0:10:20
    kk=kk+1;
    error_MSE_time=0;
	error_MSE_angle=0;
% 	for nframe=1:FrameNum
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
        
        TxBits = randint(1,Info_length*20,2);		
		ModData = data_modulate(TxBits, 'QPSK');
        OfdmSymbol_Data1 = sqrt(Info_length*20) * ifft(ModData,Info_length*20);
        
        TxBits = randint(1,Info_length*20,2);		
		ModData = data_modulate(TxBits, 'QPSK');
        OfdmSymbol_Data2 = sqrt(Info_length*20) * ifft(ModData,Info_length*20);
        
        ChannelInput=[OfdmSymbol_Data1 ChannelInput OfdmSymbol_Data2];
        
        %Pass through channel
		ChannelOutput=Frequency_offset_channel(ChannelInput,frequency_offset);
        ChannelOutput=AGWN_channel(ChannelOutput,Noise_Var);
%         ChannelOutput=Standard_Multipath1(ChannelOutput,fc,V,fs);

		for ii=1:length(ChannelOutput)-Info_length
            Gamam_P_temp=0;
			for jj=1:length(PN_code)
                Gamam_P_temp=Gamam_P_temp+PN_code(1,jj)*conj(ChannelOutput(1,jj+ii-1));
			end
            Gamam_P(1,ii)=Gamam_P_temp;
		end

        Gamam_P2=[Gamam_P zeros(1,Info_length)];
        Gamam_P3=[zeros(1,Info_length) Gamam_P];
        
        argmax=Gamam_P2.*conj(Gamam_P3);
        
		argmax2=abs(Gamam_P);
        
%         for mm=1:16
%             argmax_t(mm,:)=[argmax(1,end-m:end) argmax(1,[1:end-m])];
%         end
argmax_t(1,:)=[argmax(1,end-1:end) argmax(1,[1:end-1])];
argmax_t(2,:)=[argmax(1,end-2:end) argmax(1,[1:end-2])];
argmax_t(3,:)=[argmax(1,end-3:end) argmax(1,[1:end-3])];

argmax_t(11,:)=[argmax(1,end-11:end) argmax(1,[1:end-11])];
argmax_t(12,:)=[argmax(1,end-12:end) argmax(1,[1:end-12])];
argmax_t(13,:)=[argmax(1,end-13:end) argmax(1,[1:end-13])];

argmax_t(4,:)=[argmax(1,end-4:end) argmax(1,[1:end-4])];
argmax_t(5,:)=[argmax(1,end-5:end) argmax(1,[1:end-5])];
argmax_t(6,:)=[argmax(1,end-6:end) argmax(1,[1:end-6])];

argmax_t(14,:)=[argmax(1,end-14:end) argmax(1,[1:end-14])];
argmax_t(15,:)=[argmax(1,end-15:end) argmax(1,[1:end-15])];
argmax_t(16,:)=[argmax(1,end-16:end) argmax(1,[1:end-16])];

argmax_t(7,:)=[argmax(1,end-7:end) argmax(1,[1:end-7])];
argmax_t(8,:)=[argmax(1,end-8:end) argmax(1,[1:end-8])];
argmax_t(9,:)=[argmax(1,end-9:end) argmax(1,[1:end-9])];

argmax_t(10,:)=[argmax(1,end-10:end) argmax(1,[1:end-10])];


% argmax_t(1,:)=[Gamam_P(1,end-1:end) Gamam_P(1,[1:end-1])];
% argmax_t(2,:)=[Gamam_P(1,end-2:end) Gamam_P(1,[1:end-2])];
% argmax_t(3,:)=[Gamam_P(1,end-3:end) Gamam_P(1,[1:end-3])];
% 
% argmax_t(11,:)=[Gamam_P(1,end-11:end) Gamam_P(1,[1:end-11])];
% argmax_t(12,:)=[Gamam_P(1,end-12:end) Gamam_P(1,[1:end-12])];
% argmax_t(13,:)=[Gamam_P(1,end-13:end) Gamam_P(1,[1:end-13])];
% 
% argmax_t(4,:)=[Gamam_P(1,end-4:end) Gamam_P(1,[1:end-4])];
% argmax_t(5,:)=[Gamam_P(1,end-5:end) Gamam_P(1,[1:end-5])];
% argmax_t(6,:)=[Gamam_P(1,end-6:end) Gamam_P(1,[1:end-6])];
% 
% argmax_t(14,:)=[Gamam_P(1,end-14:end) Gamam_P(1,[1:end-14])];
% argmax_t(15,:)=[Gamam_P(1,end-15:end) Gamam_P(1,[1:end-15])];
% argmax_t(16,:)=[Gamam_P(1,end-16:end) Gamam_P(1,[1:end-16])];
% 
% argmax_t(7,:)=[Gamam_P(1,end-7:end) Gamam_P(1,[1:end-7])];
% argmax_t(8,:)=[Gamam_P(1,end-8:end) Gamam_P(1,[1:end-8])];
% argmax_t(9,:)=[Gamam_P(1,end-9:end) Gamam_P(1,[1:end-9])];
% 
% argmax_t(10,:)=[Gamam_P(1,end-10:end) Gamam_P(1,[1:end-10])];

%         argmax_t
%         size(argmax_t);
        haha=sum(argmax_t);
%         size(haha)
        
% 		argMax=max(argmax(1,Info_length*M-Info_length/4:1:Info_length*M+Info_length/4));
%         for jj=Info_length*M-Info_length/4:1:Info_length*M+Info_length/4
%             if argmax(1,jj)==argMax 
%                angle_offset=angle(Gamam_P(1,jj))*M/(2*pi)-frequency_offset+1;
%                time_offset=jj-Info_length*M-1;
%                if time_offset==0
%                else
%                   error_MSE_time=abs(time_offset)^2+error_MSE_time;
%                end
%                if angle_offset==0
%                else
%                   error_MSE_angle=abs(angle_offset)^2+error_MSE_angle;
%                end 
%             end
% 		end
% %     end
%     MSE_time(1,kk)=sqrt(error_MSE_time/FrameNum);
%     MSE_angle(1,kk)=sqrt(error_MSE_angle/FrameNum);
  
end

%Output figures

% figure(1)
% semilogy(1:kk,MSE_time,'x');
% grid on;
% hold on;
% 
% figure(2)
% semilogy(1:kk,MSE_angle,'x');
% grid on;
% hold on;

% figure(3)
% correlation_flag=length(ChannelOutput); 
% [R_xx,lags] = xcorr(ChannelOutput,correlation_flag,'coeff');
% subplot(211);
% plot(1:length(ChannelOutput),real(ChannelOutput));
% grid on;
% 
% subplot(212);
% plot(lags,real(R_xx))
% grid on;

% figure(4)
% subplot(211);
% plot(1:length(Gamam_PP),Gamam_PP);
% grid on;
% 
% subplot(212);
% for m = 1:M
%     plot(0:(Info_length-1),Gamam_PP(1,(m-1)*Info_length+1:m*Info_length));  
%     hold on
% end
% grid on;

% figure(5)
% subplot(211);
% plot(1:length(argmax),argmax);
% grid on;
% 
% subplot(212);
% for m = 2:M
%     plot(0:(Info_length-1),argmax(1,(m-1)*Info_length+1:m*Info_length));  
%     hold on
% end
% grid on;

figure(6)
subplot(211);
plot(1:length(Gamam_P),Gamam_P);
grid on;

subplot(212);
plot(1:length(argmax2),argmax2);
grid on;

figure(7)
subplot(211);
plot(1:length(argmax),abs(argmax));
grid on;

subplot(212);
plot(1:length(haha),abs(haha));
grid on;

TimeElapse = toc/60;	% Simulation minutes
fprintf('****** The total time used is : %3.4f minutes ********\n', TimeElapse);