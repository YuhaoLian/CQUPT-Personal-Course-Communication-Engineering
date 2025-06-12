close all;
clear;
clc;
fc = 2e9;			%Carrier frequency
fs = 20e+6; 		%System bandwith is 20MHz
V=0;
FrameNum=100;
SubCarrNum = 128;	%size fourier transform to generate signal. It
					%is equal to the number of samples in the OFDM
					%symbol, and also it is the size of FFT.
GuardNum = 16;	    %Total guard time in samples.
ModType=2;          %for QPSK
Info_length=SubCarrNum*ModType;
frequency_offset=0.25;
frequency_adjustment=(SubCarrNum+GuardNum)/SubCarrNum;

est_time=zeros(1,20);
est_angle=zeros(1,20);
    
for snr=1:1:20

	SNR=10^(2*snr/10);
	E=1;                      % 符号能量假设为单位能量.
	Noise_Var = E/(ModType*SNR);

	for nframe = 1:FrameNum
	
        TxBits = randint(1,Info_length,2);
        ModData = data_modulate(TxBits, 'QPSK');
        % IFFT
        OfdmSymbol = sqrt(SubCarrNum) * ifft(ModData,SubCarrNum);
	
        % Add the cyclic prefix to the ofdm symbol
        OfdmSymbol = [OfdmSymbol(:,SubCarrNum-GuardNum+1:SubCarrNum) OfdmSymbol];
        
        %Add frequency offset
        OfdmSymbol=Frequency_offset_channel(OfdmSymbol,frequency_offset);
        
        if nframe==1 
           ChannelInput=OfdmSymbol;
        else
           ChannelInput=[ChannelInput OfdmSymbol];    
        end
        
	end 
	
	%pass through channel
	ChannelOutput=AGWN_channel(ChannelInput,Noise_Var);
%     if multipath_yn==1
%        ChannelOutput=Standard_Multipath1(ChannelOutput,fc,V,fs);
%     end
    
	%信号移位相乘
	ChannelOutput2=[zeros(1,SubCarrNum) ChannelOutput(:,1:end-SubCarrNum)];
	for i=1:length(ChannelOutput2)
        if i==1
           dd=ChannelOutput2(:,i)*conj(ChannelOutput(:,i));
        else
           dd=[dd ChannelOutput2(:,i)*conj(ChannelOutput(:,i));];
        end
	end
	
	%求相关
	for i=1:length(dd)-GuardNum
        s_temp=0;
        for j=1:GuardNum
            s_temp=s_temp+dd(:,(i+j-1));
        end
        if i==1
           ddd=s_temp;
        else
           ddd=[ddd s_temp];
        end
	end
    ddd1=abs(ddd);
    
    error_time=0;
    error_angle=0;
    for i=1:FrameNum-1
        argmax=max(ddd1(1,((i-1)*(SubCarrNum+GuardNum)+SubCarrNum-10):1:((i-1)*(SubCarrNum+GuardNum)+SubCarrNum+10)));
        for jj=((i-1)*(SubCarrNum+GuardNum)+SubCarrNum-10):1:((i-1)*(SubCarrNum+GuardNum)+SubCarrNum+10)
            if ddd1(1,jj)==argmax             
               time_offset=jj-1-SubCarrNum-(i-1)*(SubCarrNum+GuardNum);
               angle_offset=frequency_adjustment*angle(ddd(1,jj))/(2*pi)+frequency_offset;
               error_time=abs(time_offset)^2+error_time;
               error_angle=abs(angle_offset)^2+error_angle;
            end
        end
	end
    
	MSE_time=sqrt(error_time/(FrameNum-1));
    est_time(1,snr)=MSE_time;
	
    MSE_angle=sqrt(error_angle/(FrameNum-1));
    est_angle(1,snr)=MSE_angle;
    
%     figure(1)
% 	subplot(211);
% 	plot(1:length(ddd1(1,SubCarrNum+1+GuardNum:end-GuardNum)),abs(ddd1(1,SubCarrNum+1+GuardNum:end-GuardNum)));
% 	grid on;
% 	
% 	subplot(212);
% 	plot(1:length(ddd(1,SubCarrNum+1+GuardNum:end-GuardNum)),angle(ddd(1,SubCarrNum+1+GuardNum:end-GuardNum)));
% 	grid on;

end

figure(2)
semilogy(2*(1:1:length(est_time)),est_time,'-o');
grid on;
hold on;

semilogy(2*(1:1:length(est_angle)),est_angle,'--s');
grid on;
hold on;