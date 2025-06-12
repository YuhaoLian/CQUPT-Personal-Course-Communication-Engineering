% Simulation of《ML Estimation of Time and Frequency Offset in OFDM Systems》
close all;
clear;
clc;
tic;
fc = 2.0e9;			%Carrier frequency
fs = 20e+6; 			%System bandwith is 5MHz
V=100;
FrameNum=10000;
SubCarrNum = 128;	%Number of carriers
GuardNum=16;        %Guard time in samples.
ModType=2;          %For QPSK
Info_length=SubCarrNum*ModType;
frequency_offset=0.25;

for mutipath_yn=0:1:1
    
	est_angle=zeros(1,11);
	est_time=zeros(1,11);
 
    for ss=0:1:10      %dB
        
    SNR=2*ss;
    Rho=SNR/(SNR+1);
	snr=10^(SNR/10);
    E=1;                      % 符号能量假设为单位能量.
	Noise_Var = E/(2*ModType*snr);
		
        for nframe = 1:FrameNum
		
            TxBits = randint(1,Info_length,2);
            
            ModData = data_modulate(TxBits, 'QPSK');
            
            % IFFT
            OfdmSymbol = sqrt(SubCarrNum) * ifft(ModData,SubCarrNum);
		
            % Add the cyclic prefix to the ofdm symbol
            ChannelInput = [OfdmSymbol(:,SubCarrNum-GuardNum+1:SubCarrNum) OfdmSymbol];
            
            %Add frequency offset
            ChannelOutput_temp=Frequency_offset_channel(ChannelInput,frequency_offset);
                        
            if nframe==1 
               ChannelOutput=ChannelOutput_temp;
            else
               ChannelOutput=[ChannelOutput ChannelOutput_temp];    
            end
		end 
        
        %Pass through channel
        ChannelOutput=AGWN_channel(ChannelOutput,Noise_Var);
        if mutipath_yn==1
           ChannelOutput=Standard_Multipath1(ChannelOutput,fc,V,fs);    
        end
		
		%信号移位相乘
		ChannelOutput_delayN=[zeros(1,SubCarrNum) ChannelOutput(:,1:end-SubCarrNum)];
		before_Feim=conj(ChannelOutput_delayN).*ChannelOutput;
		
		%求相关
        Feim=zeros(1,length(before_Feim)-GuardNum);
		Gamam=zeros(1,length(before_Feim)-GuardNum);
		for i=1:length(before_Feim)-GuardNum
            Gamam_temp=0;
            Feim_temp=0;
            for j=1:GuardNum
                Gamam_temp=Gamam_temp+before_Feim(:,(i+j-1));
                Feim_temp=Feim_temp+(abs(ChannelOutput(:,(i+j-1))))^2+(abs(ChannelOutput_delayN(:,(i+j-1))))^2;
            end
            Gamam(1,i)=Gamam_temp;
            Feim(1,i)=Feim_temp;
		end
		
		argmax_before=abs(Gamam)-(Rho*Feim)/2;
		
        frequency_adjustment=(SubCarrNum+GuardNum)/SubCarrNum;
		error_time=0;
		error_angle=0;
		time_offset=zeros(1,FrameNum-2);
		angel_offset=zeros(1,FrameNum-2);
		for i=2:FrameNum-1
            argmax=max(argmax_before(1,(i*(SubCarrNum+GuardNum)-50):1:(i*(SubCarrNum+GuardNum)+50)));
            for jj=(i*(SubCarrNum+GuardNum)-50):1:(i*(SubCarrNum+GuardNum)+50)
                if argmax_before(1,jj)==argmax
                   error_time=abs(jj-1-SubCarrNum-(i-1)*(SubCarrNum+GuardNum))^2+error_time;
                   
                   error_angle=abs(frequency_adjustment*angle(Gamam(jj))/(2*pi)-frequency_offset)^2+error_angle;
                   
                   time_offset(1,i-1)=jj-1-SubCarrNum-(i-1)*(SubCarrNum+GuardNum);
                   
                   angel_offset(1,i-1)=frequency_adjustment*angle(Gamam(jj))/(2*pi)-frequency_offset;
                end
            end
		end
		MSE_time=sqrt(error_time/(FrameNum-2));
		MSE_angle=sqrt(error_angle/(FrameNum-2));
        est_angle(1,ss+1)=MSE_angle;
        est_time(1,ss+1)=MSE_time;
	end

    if mutipath_yn==0
       picture_style='-o';
    else
       picture_style='--x';
    end
    
    figure(1)
	semilogy(0:2:20,est_time,picture_style);
	grid on;
	hold on;
    
    figure(2)
	semilogy(0:2:20,est_angle,picture_style);
	grid on;
    hold on;
    
end

% figure(3)
% subplot(211);
% plot(1:length(argmax_before),argmax_before);
% grid on;
% 
% subplot(212);
% plot(1:length(Gamam),angle(Gamam));
% grid on;
% 
% figure(4)
% subplot(211);
% hist(time_offset,100);
% grid on;
% 
% subplot(212);
% hist(angel_offset,100);
% grid on;

TimeElapse = toc/60;	% Simulation minutes
fprintf('****** The total time used is : %3.4f minutes ********\n', TimeElapse);