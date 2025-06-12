% Simulation of《Symbol time offset estimation in choerent OFDM systems》
close all;
clear;
clc;
tic;
fc = 2.4e9;			%Carrier frequency
fs = 5e+6; 			%System bandwith is 5MHz
V=5;
FrameNum=10;
SubCarrNum = 1024;	%Number of carriers
GuardNum = 0;       %Number of cyclic prefix
Rho_power=50;
frequency_offset=3.2;

MF_container=[1024];
for Mf_id = 1:1:length(MF_container)
    
    %Pilot information
    Mf=MF_container(1,Mf_id);
	PilotCarrNum = ceil( SubCarrNum/Mf );
	PilotSymbol = (1+j)/sqrt(2);
	PilotIndex = 1:Mf:SubCarrNum;
	TmpIndex = zeros(1,SubCarrNum);
	TmpIndex(PilotIndex) = PilotIndex;
	DataIndex = find((1:SubCarrNum) ~= TmpIndex);
	clear('TmpIndex');
    
    TransSymbol_Pilot = ones(1,SubCarrNum);
    TransSymbol_Pilot(DataIndex) = zeros(1,length(SubCarrNum));
	TransSymbol_Pilot(PilotIndex) = TransSymbol_Pilot(PilotIndex)*PilotSymbol;
    
	% IFFT of pilot datas
    OfdmSymbol_Pilot = sqrt(SubCarrNum) * ifft(TransSymbol_Pilot,SubCarrNum);
    
    %Add the cyclic prefix to the pilot symbol
	OfdmSymbol_Pilot = [OfdmSymbol_Pilot(:,SubCarrNum-GuardNum+1:SubCarrNum) OfdmSymbol_Pilot];
            
    PN_code1=PN_Signal(length(find(OfdmSymbol_Pilot(1:end) ~= 0))/4);
    length(PN_code1)
    PN_code=PN_code1;
    for ii=1:4
        PN_code=[PN_code PN_code];
    end
%       PN_code=AGWN_Generator(length(find(OfdmSymbol_Pilot(1:end) ~= 0)));
        
    jj=1;
    for ii=1:length(OfdmSymbol_Pilot)
        if OfdmSymbol_Pilot(1,ii) ~= 0
           OfdmSymbol_Pilot_I(1,ii)=real(OfdmSymbol_Pilot(1,ii))*real(PN_code(1,jj));
           OfdmSymbol_Pilot_Q(1,ii)=imag(OfdmSymbol_Pilot(1,ii))*imag(PN_code(1,jj));
           OfdmSymbol_Pilot(1,ii)=OfdmSymbol_Pilot_I(1,ii)+sqrt(-1)*OfdmSymbol_Pilot_Q(1,ii);
           jj=1+jj;
        end       
    end
    
    ChannelInput_Pilot=sqrt(Rho_power)*OfdmSymbol_Pilot;
    
    ModType=2;          %For QPSK
	Info_length=length(DataIndex) * ModType;
        
	for SNR=1:1:1
		Rho=SNR/(SNR+1);
		snr=10^(SNR/10);
		E=1;                % 符号能量假设为单位能量.
		Noise_Var = E/(2*ModType*snr);
		       
		for nframe = 1:FrameNum
		
			TxBits = randint(1,Info_length,2);
			
			ModData = zeros(1,length(TxBits)/2);%data_modulate(TxBits, 'QPSK');
			
			% Insert pilot data
			TransSymbol_Data = zeros(1,SubCarrNum);
			TransSymbol_Data(DataIndex) = ModData;
            			
			% IFFT
			OfdmSymbol_Data = sqrt(SubCarrNum) * ifft(TransSymbol_Data,SubCarrNum);
            
			% Add the cyclic prefix to the ofdm symbol
			ChannelInput_Data = [OfdmSymbol_Data(:,SubCarrNum-GuardNum+1:SubCarrNum) OfdmSymbol_Data];
			                                 
            ChannelInput=ChannelInput_Pilot;% + ChannelInput_Data;
            
            ChannelOutput_temp=Frequency_offset_channel(ChannelInput,frequency_offset);
            
%             size(ChannelOutput_temp.*ChannelInput)
%             sum(angle(conj(ChannelOutput_temp).*ChannelInput))/(pi*length(ChannelOutput_temp))
            
			if nframe==1 
               ChannelOutput=ChannelOutput_temp;
			else
               ChannelOutput=[ChannelOutput ChannelOutput_temp];    
			end
		end 
		
		
		%Pass through channel
% 		ChannelOutput=AGWN_channel(ChannelInput,Noise_Var);
%         ChannelOutput=ChannelInput;
%         ChannelOutput=Frequency_offset_channel(ChannelInput,frequency_offset);
		%ChannelOutput_temp=Standard_Multipath1( ChannelInput,fc,V,fs);
		
		
		%信号移位相乘
% 		ChannelOutput_delayN=[zeros(1,SubCarrNum) ChannelOutput(:,1:end-SubCarrNum)];
% 		before_Feim_CP=ChannelOutput_delayN.*conj(ChannelOutput);
		
% 		for i=1:length(before_Feim_CP)-GuardNum
% 			Gamam_PP_temp=0;
% 			Feim_CP_temp=0;
% 			Feim_P_temp=0;
% 			for j=1:GuardNum
%                 Gamam_PP_temp=Gamam_PP_temp+before_Feim_CP(:,(i+j-1));
%                 Feim_CP_temp=Feim_CP_temp+(abs(ChannelOutput(:,(i+j-1))))^2+(abs(ChannelOutput_delayN(:,(i+j-1))))^2;
%                 Feim_P_temp=Feim_P_temp+conj(ChannelOutput(:,(i+j-1))+ChannelOutput_delayN(:,(i+j-1)))*ChannelInput_Pilot(1,j);
% 			end
% 			Gamam_PP(1,i)=Gamam_PP_temp;
% 			Feim_CP(1,i)=Feim_CP_temp;
% 			Feim_P(1,i)=Feim_P_temp;
% 		end
			
		ChannelOutput2=[zeros(1,length(ChannelInput_Pilot)) ChannelOutput];
			
%         length(PN_code1)
%          length(ChannelInput_Pilot/4*8)       
		for ii=1:length(ChannelOutput2)-length(ChannelInput_Pilot);
            Gamam_P_temp=0;
            for mm=1:4
				for jj=1:length(PN_code1)
                    Gamam_P_temp=Gamam_P_temp+conj(PN_code1(1,jj))*ChannelOutput2(1,jj*mm+ii-1);
				end
				Gamam_P(1,ii)=Gamam_P_temp;
            end
		end
		
        Gamam_P=Gamam_P(1,GuardNum+1:end);
        
        Gamam_P1=Gamam_P;
        Gamam_P2 =[zeros(1,length(PN_code1)) Gamam_P(1,1:end-(length(PN_code1)))];
        
        Gamam_P=Gamam_P1.*conj(Gamam_P2);
        
		argmax_before=abs(Gamam_P);
        
% 		argmax_before_CP=real(Gamam_PP)-(Rho*Feim_CP)/2;
% 		argmax_before_P=(1+Rho)*real(Gamam_P)-Rho*real(Feim_P);
% 		argmax_before=(1-Rho)*argmax_before_P+Rho*argmax_before_CP;
        
		error_time=0;
		error_ratio_time=0;
		for i=2:FrameNum-1
			[argmax jj]=max(argmax_before(1,(i*(SubCarrNum+GuardNum)-GuardNum-10):1:(i*(SubCarrNum+GuardNum)-GuardNum+10)));
            for jj=(i*(SubCarrNum+GuardNum)-GuardNum-50):1:(i*(SubCarrNum+GuardNum)-GuardNum+50)
                if argmax_before(1,jj)==argmax 
                   haha=(SubCarrNum+GuardNum)/length(PN_code1)
                   jj-256*3
                   angle(Gamam_P(1,jj-256*3))*haha/(2*pi)
                   time_offset(1,i-1)=jj-1-i*(SubCarrNum+GuardNum)+GuardNum;
                   error_time=abs(time_offset(1,i-1))^2+error_time;
                   if time_offset(1,i-1)==0
                   else
                      error_ratio_time=error_ratio_time+1;
                   end
                end
			end
		end
		
		MSE_time(1,SNR+1)=sqrt(error_time/(FrameNum-2));
		Error_ratio_time(1,SNR+1)=error_ratio_time/(FrameNum-2)
	end
	
    if Mf_id==1
       picture_style='-o';
    elseif Mf_id==2
       picture_style='-^';
%     elseif Mf_id==3
%        picture_style='-v';
    elseif Mf_id==3
       picture_style='-*';
    elseif Mf_id==4
       picture_style='-x';
    end
    
% 	figure(1)
% 	semilogy(0:1:5,Error_ratio_time,picture_style);
% 	grid on;
% 	hold on;
   
end

%Output figures
figure(2)
correlation_flag=length(PN_code); 
[R_xx,lags] = xcorr(PN_code,correlation_flag,'coeff');
subplot(211);
plot(1:length(PN_code),real(PN_code));
grid on;

subplot(212);
plot(lags,real(R_xx))
grid on;

figure(3)
correlation_flag=length(ChannelInput_Pilot); 
[R_xx,lags] = xcorr(ChannelInput_Pilot,correlation_flag,'coeff');
subplot(211);
plot(1:length(ChannelInput_Pilot),real(ChannelInput_Pilot));
grid on;

subplot(212);
plot(lags,real(R_xx))
grid on;

% figure(4)
% subplot(211);
% plot(1:length(argmax_before_CP),argmax_before_CP);
% grid on;
% 
% subplot(212);
% for nframe = 2:FrameNum-1
%     plot(0:(SubCarrNum+GuardNum-1),argmax_before_CP(1,(nframe-1)*(SubCarrNum+GuardNum)+1:nframe*(SubCarrNum+GuardNum)));  
%     hold on
% end
% grid on;
% 
% figure(5)
% subplot(211);
% plot(1:length(argmax_before_P),argmax_before_P);
% grid on;
% 
% subplot(212);
% for nframe = 2:FrameNum-1
%     plot(0:(SubCarrNum+GuardNum-1),argmax_before_P(1,(nframe-1)*(SubCarrNum+GuardNum)+1:nframe*(SubCarrNum+GuardNum)));  
%     hold on
% end
% grid on;

figure(6)
subplot(211);
plot(1:length(argmax_before),argmax_before);
grid on;

subplot(212);
for nframe = 2:FrameNum-1
    plot(0:(SubCarrNum+GuardNum-1),argmax_before(1,(nframe-1)*(SubCarrNum+GuardNum)+1:nframe*(SubCarrNum+GuardNum)));  
    hold on
end
grid on;

% figure(7)
% subplot(311);
% for nframe = 2:FrameNum-1
%     plot(0:(SubCarrNum+GuardNum-1),argmax_before_CP(1,(nframe-1)*(SubCarrNum+GuardNum)+1:nframe*(SubCarrNum+GuardNum)));  
%     hold on
% end
% grid on;
% 
% subplot(312);
% for nframe = 2:FrameNum-1
%     plot(0:(SubCarrNum+GuardNum-1),argmax_before_P(1,(nframe-1)*(SubCarrNum+GuardNum)+1:nframe*(SubCarrNum+GuardNum)));  
%     hold on
% end
% grid on;
% 
% subplot(313);
% for nframe = 2:FrameNum-1
%     plot(0:(SubCarrNum+GuardNum-1),argmax_before(1,(nframe-1)*(SubCarrNum+GuardNum)+1:nframe*(SubCarrNum+GuardNum)));  
%     hold on
% end
% grid on;

TimeElapse = toc/60;	% Simulation minutes
fprintf('****** The total time used is : %3.4f minutes ********\n', TimeElapse);