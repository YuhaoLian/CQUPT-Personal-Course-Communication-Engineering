% Simulation of《ML Estimation of Time and Frequency Offset in OFDM Systems》
close all;
clear;
clc;
tic;
fc = 2.4e9;			%Carrier frequency
fs = 5e+6; 			%System bandwith is 5MHz
V=5;
FrameNum=10;
SubCarrNum = 128;	%Number of carriers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mf = 5;
PilotCarrNum = ceil( SubCarrNum/Mf );
PilotSymbol = (1+j)/sqrt(2);
PilotIndex = 1:Mf:SubCarrNum;
TmpIndex = zeros(1,SubCarrNum);
TmpIndex(PilotIndex) = PilotIndex;
DataIndex = find((1:SubCarrNum) ~= TmpIndex);
clear('TmpIndex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ModType=2;          %For QPSK
Info_length=length(DataIndex) * ModType;
frequency_offset=0.25;

% for ss=1:1:4      %dB
    
    SNR=6*ss;
    Rho=1;SNR/(SNR+1);
% 	snr=10^(SNR/10);
    E=1;                      % 符号能量假设为单位能量.
% 	Noise_Var = E/(2*ModType*snr);
	
	est_angle=zeros(1,17);
	est_time=zeros(1,17);
    
% 	for gg = 0:1:16    %Guard time in samples.
        
        GuardNum=16;%gg;
		
        for nframe = 1:FrameNum
		
            TxBits = randint(1,Info_length,2);
            
            ModData = data_modulate(TxBits, 'QPSK');
            
            % Insert pilot data
		    TransSymbol = ones(1,SubCarrNum);
		    TransSymbol(DataIndex) = ModData;
		    TransSymbol(PilotIndex) = TransSymbol(PilotIndex)*PilotSymbol;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Reference_sinal=ones(1,SubCarrNum);
		    Reference_sinal(DataIndex) = zeros(1,length(ModData));
		    Reference_sinal(PilotIndex) = Reference_sinal(PilotIndex)*PilotSymbol;
            
            % IFFT
            OfdmSymbol = sqrt(SubCarrNum) * ifft(TransSymbol,SubCarrNum);
            M_sinal_Input_temp = sqrt(SubCarrNum) * ifft(Reference_sinal,SubCarrNum);
		
            % Add the cyclic prefix to the ofdm symbol
            ChannelInput = [OfdmSymbol(:,SubCarrNum-GuardNum+1:SubCarrNum) OfdmSymbol];
            M_sinal_Input= [M_sinal_Input_temp(:,SubCarrNum-GuardNum+1:SubCarrNum) M_sinal_Input_temp];
            
            %Pass through channel
%             ChannelOutput_temp=Frequency_offset_channel(ChannelInput,frequency_offset);
%             ChannelOutput_temp=AGWN_channel(ChannelOutput_temp,Noise_Var);
%             ChannelOutput_temp=Standard_Multipath1( ChannelInput,fc,V,fs);
            ChannelOutput_temp=ChannelInput;
            M_sinal_bb= M_sinal_Input;
            
            if nframe==1 
               ChannelOutput=ChannelOutput_temp;
               M_sinal=M_sinal_bb;
            else
               ChannelOutput=[ChannelOutput ChannelOutput_temp];    
               M_sinal=[M_sinal M_sinal_bb];
            end
		end 
		
		%信号移位相乘
		ChannelOutput_delayN=[zeros(1,SubCarrNum) ChannelOutput(:,1:end-SubCarrNum)];
		before_Feim=ChannelOutput_delayN.*conj(ChannelOutput);
		
		%求相关
        Feim=zeros(1,length(before_Feim)-GuardNum);
		Gamam=zeros(1,length(before_Feim)-GuardNum);
        hehe=zeros(1,length(before_Feim)-GuardNum);
		for i=1:length(before_Feim)-GuardNum
            Gamam_temp=0;
            Feim_temp=0;
            hehe_temp=0;
            for j=1:GuardNum
                Gamam_temp=Gamam_temp+before_Feim(:,(i+j-1));
                Feim_temp=Feim_temp+(abs(ChannelOutput(:,(i+j-1))))^2+(abs(ChannelOutput_delayN(:,(i+j-1))))^2;
                hehe_temp=hehe_temp+conj(ChannelOutput(:,(i+j-1))+ChannelOutput_delayN(:,(i+j-1)))*M_sinal(1,j);
            end
            Gamam(1,i)=Gamam_temp;
            Feim(1,i)=Feim_temp;
            hehe(1,i)=hehe_temp;
		end
        
M_sinal1=M_sinal(1,1:SubCarrNum+GuardNum);
correlation_flag=length(M_sinal1); 
[R_xx,lags] = xcorr(M_sinal1,correlation_flag,'coeff');

M_sinal2=[zeros(1,length(M_sinal1)) M_sinal];
haha=zeros(1,length(M_sinal2)-length(M_sinal1));

for ii=1:length(M_sinal2)-length(M_sinal1);
    kkk_temp=0;
    for jj=1:length(M_sinal1)
        kkk_temp=kkk_temp+M_sinal1(1,jj)*conj(M_sinal2(1,jj+ii-1));
    end
    haha(1,ii)=kkk_temp;
end
        
ChannelOutput2=[zeros(1,length(M_sinal1)) ChannelOutput];
haha2=zeros(1,length(ChannelOutput2)-length(M_sinal1));

for ii=1:length(ChannelOutput2)-length(M_sinal1);
    kkk_temp=0;
    for jj=1:length(M_sinal1)
        kkk_temp=kkk_temp+M_sinal1(1,jj)*conj(ChannelOutput2(1,jj+ii-1));
    end
    haha2(1,ii)=kkk_temp;
end

subplot(211);
plot(1:length(M_sinal),real(M_sinal));
grid on;

subplot(212);
plot(lags,R_xx)
grid on;

figure
subplot(211);
plot(1:length(haha2),haha2)
grid on;
grid on;
		
subplot(212);
plot(1:length(haha),haha)
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%dd
% %         hoho=zeros(1,length(before_Feim)-GuardNum-SubCarrNum);
% % 		for i=1:length(before_Feim)-SubCarrNum-GuardNum
% %             hoho_temp=0;
% %             for j=1:GuardNum+SubCarrNum
% %                 hoho_temp=hoho_temp+conj(ChannelOutput(:,(i+j-1)))*M_sinal(1,j);
% %             end
% %             hoho(1,i)=hoho_temp;
% % 		end
%         
%         
%         
% 		argmax_before=real(Gamam)-(Rho*Feim)/2;
%         argmax_before2=(1+Rho)*real(hoho)-Rho*real(hehe);
% 		
%         frequency_adjustment=(SubCarrNum+GuardNum)/SubCarrNum;
% 		error_time=0;
% 		error_angle=0;
% 		time_offset=zeros(1,FrameNum-2);
% 		angel_offset=zeros(1,FrameNum-2);
% 		for i=2:FrameNum-1
%             argmax=max(argmax_before(1,(i*(SubCarrNum+GuardNum)-100):1:(i*(SubCarrNum+GuardNum)+100)));
%             for jj=(i*(SubCarrNum+GuardNum)-100):1:(i*(SubCarrNum+GuardNum)+100)
%                 if argmax_before(1,jj)==argmax
%                    error_time=abs(jj-1-SubCarrNum-(i-1)*(SubCarrNum+GuardNum))^2+error_time;
%                    
%                    error_angle=abs(frequency_adjustment*angle(Gamam(jj))/(2*pi)-frequency_offset)^2+error_angle;
%                    
%                    time_offset(1,i-1)=jj-1-SubCarrNum-(i-1)*(SubCarrNum+GuardNum);
%                    
%                    angel_offset(1,i-1)=frequency_adjustment*angle(Gamam(jj))/(2*pi)-frequency_offset;
%                 end
%             end
% 		end
% 		MSE_time=sqrt(error_time/(FrameNum-2));
% 		MSE_angle=sqrt(error_angle/(FrameNum-2));
% %         est_angle(1,1+gg)=MSE_angle;
% %         est_time(1,1+gg)=MSE_time;
% % 	end
% % 	
% %     if ss==1
% %        picture_style='-x';
% %     elseif ss==2
% %        picture_style='-*';
% %     elseif ss==3
% %        picture_style='-^';
% %     elseif ss==4
% %        picture_style='-o';
% %     end
% % 	
% %     figure(1)
% % 	subplot(211);
% % 	semilogy(0:1:16,est_time,picture_style);
% % 	grid on;
% % 	hold on;
% %     
% % 	subplot(212);
% % 	semilogy(0:1:16,est_angle,picture_style);
% % 	grid on;
% %     hold on;
% %     
% %     save est_time1.mat est_time;
% %     save est_angle1.mat est_angle;
% % end
% figure(1)
% subplot(211);
% plot(1:length(argmax_before),argmax_before);
% grid on;
% 
% % subplot(212);
% % plot(1:length(Gamam),angle(Gamam));
% % grid on;
% length(argmax_before)
% % argmax_before1=zeros(1)
% % argmax_before1=argmax_before(1:end);
% % figure(3)
% subplot(212);
% for nframe = 1:FrameNum-1
%     plot(0:(SubCarrNum+GuardNum-1),argmax_before(1,(nframe-1)*(SubCarrNum+GuardNum)+1:nframe*(SubCarrNum+GuardNum)));  
%     hold on
% end
% grid on;
% 
% figure(2)
% subplot(211);
% plot(1:length(argmax_before2),argmax_before2);
% grid on;
% 
% % subplot(212);
% % plot(1:length(Gamam),angle(Gamam));
% % grid on;
% length(argmax_before2)
% % argmax_before1=zeros(1)
% % argmax_before1=argmax_before2(1:end);
% % figure(3)
% subplot(212);
% for nframe = 1:FrameNum-1
%     plot(0:(SubCarrNum+GuardNum-1),argmax_before2(1,(nframe-1)*(SubCarrNum+GuardNum)+1:nframe*(SubCarrNum+GuardNum)));  
%     hold on
% end
% grid on;
% % figure(3)
% % subplot(211);
% % hist(time_offset,100);
% % axis auto
% % grid on;
% % 
% % subplot(212);
% % hist(angel_offset,100);
% % axis auto
% % grid on;
% % TimeElapse = toc/60	% Simulation minutes
% % fprintf('****** The total time that used is : %3.4f minute ********\n', TimeElapse);