close all;
clear;
clc;
fc = 2e9;			%Carrier frequency
fs = 20e+6; 		%System bandwith is 20MHz
V=0;
FrameNum=10000;
SubCarrNum = 128;	%size fourier transform to generate signal. It
					%is equal to the number of samples in the OFDM
					%symbol, and also it is the size of FFT.
GuardNum = 32;	    %Total guard time in samples.
ModType=4; %for 16QAM
Info_length=SubCarrNum*ModType;

est_time=zeros(1,20);

for mutilpath_yn=0:1:1

	for snr=1:1:20
	
		SNR=10^(2*snr/10);
		E=1;                      % 符号能量假设为单位能量.
		Noise_Var = E/(2*ModType*SNR);
		  
		for nframe = 1:FrameNum
		
            TxBits = randint(1,Info_length,2);
            ModData = data_modulate(TxBits, '16QAM');
            % IFFT
            OfdmSymbol = sqrt(SubCarrNum) * ifft(ModData,SubCarrNum);
		
            % Add the cyclic prefix to the ofdm symbol
            OfdmSymbol = [OfdmSymbol(:,SubCarrNum-GuardNum+1:SubCarrNum) OfdmSymbol];
            
           if mutilpath_yn==1
              OfdmSymbol=Standard_Multipath1(OfdmSymbol,fc,V,fs);
           end
            
            if nframe==1 
               ChannelInput=OfdmSymbol;
            else
               ChannelInput=[ChannelInput OfdmSymbol];    
            end
		end 
		
		%pass through channel
		ChannelOutput=AGWN_channel(ChannelInput,Noise_Var);
		
        %信号移位相减
		ChannelOutput2=[zeros(1,SubCarrNum) ChannelOutput(:,1:end-SubCarrNum)];
		dd=ChannelOutput2-ChannelOutput;
%         figure(1)
% 		subplot(311)
% 		plot(1:length(abs(dd)),abs(dd));
		
		%滑动平均
		Q=10;
		for i = 1:Q
            if i>1
                Temp_Sig = [ abs(dd(1,[i:end])) zeros(1,i-1)];  
                dd_temp= dd_temp+Temp_Sig;
            else     
                dd_temp = abs(dd);
            end
		end  
		ddd_temp=dd_temp/Q;
% 		subplot(312)
% 		plot(1:length(ddd_temp),ddd_temp);
		
		%求突变点
		for i=1:length(ddd_temp)-1
            if i==1
            mmm=abs(ddd_temp(1,i+1))/abs(ddd_temp(1,i));
            else
            mmm=[mmm abs(ddd_temp(1,i+1))/abs(ddd_temp(1,i))];
            end
		end       
% 		subplot(313)
% 		plot(1:length(mmm),mmm);
        
        error_time=0;
        for i=1:FrameNum-1
            argmax=max(mmm(1,(i*(SubCarrNum+GuardNum)-Q-8):1:(i*(SubCarrNum+GuardNum)-Q+8)));
            for jj=(i*(SubCarrNum+GuardNum)-Q-8):1:(i*(SubCarrNum+GuardNum)-Q+8)
                if mmm(1,jj)==argmax             
                   time_offset=jj+Q-i*(SubCarrNum+GuardNum)-1;
                   error_time=abs(time_offset)^2+error_time;
                end
            end
		end
        
		MSE_time=sqrt(error_time/(FrameNum-1));
        est_time(1,snr)=MSE_time;
	
	end

    if mutilpath_yn==0
       picture_style='-o';
    else
       picture_style='--*';
    end
    
	figure(2)
	semilogy(2*(1:1:length(est_time)),est_time,picture_style);
	grid on;
	hold on;
    
end