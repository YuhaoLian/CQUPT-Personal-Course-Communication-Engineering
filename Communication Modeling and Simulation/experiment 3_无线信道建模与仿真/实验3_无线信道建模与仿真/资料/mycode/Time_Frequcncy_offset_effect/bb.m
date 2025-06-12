close all;
clear;
clc;
fc = 2.4e9;			%Carrier frequency
fs = 5e+6; 			%System bandwith is 5MHz
V=88;
EbN0db = 0:1:5;		% Eb/N0
FrameNum=500;
SubCarrNum = 128;	%size fourier transform to generate signal. It
					%is equal to the number of samples in the OFDM
					%symbol, and also it is the size of FFT.
GuardNum = 16;	%Total guard time in samples.
ModType=2;      %for QPSK
Info_length=SubCarrNum*ModType;
errs_rate=zeros(1,length(EbN0db));

for frequcncy_offset_temp=1:1:4
    
    if frequcncy_offset_temp==0
       frequcncy_offset=0;
    else
       frequcncy_offset=0.1*(frequcncy_offset_temp-1);
    end   
    
	for nEN = 1:length(EbN0db)
        E=1;     % 符号能量假设为单位能量.
        snr = 10^(EbN0db(nEN)/10);% en为bit信噪比；ModType*en为符号信噪比.
		Noise_Var = E/(2*ModType*snr);
        
		for nframe = 1:FrameNum
		
            TxBits = randint(1,Info_length,2);
	
            ModData = data_modulate(TxBits, 'QPSK');
           
            % IFFT
            OfdmSymbol = sqrt(SubCarrNum) * ifft(ModData,SubCarrNum);
	
            % Add the cyclic prefix to the ofdm symbol
            ChannelInput = [OfdmSymbol(:,SubCarrNum-GuardNum+1:SubCarrNum) OfdmSymbol];
            
            if nframe==1 
               TxBits_temp=TxBits;
            else
               TxBits_temp=[TxBits_temp TxBits];
            end
        
            %Add frequency offset
		    ChannelOutput=Frequency_offset_channel(ChannelInput,frequcncy_offset);
            
            %Pass through channel
            ChannelOutput=AGWN_channel(ChannelOutput,Noise_Var);
	
		    RecTimeSymbol = ChannelOutput(:,GuardNum+1:end);
            
            % The FFT of the time domain signal after the removal of cyclic prefix
            RecSymbol= (1/sqrt(SubCarrNum)) * fft(RecTimeSymbol,SubCarrNum);
            
            %phase_adjustment=exp(-1*sqrt(-1)*pi*frequcncy_offset*(length(ChannelInput)-1)/length(ChannelInput));
            if nframe==1 
               RecSymbol_temp=RecSymbol;
            else
               RecSymbol_temp=[RecSymbol_temp RecSymbol];
            end
		end 
        
		%Decide
        RecSymbol_last=mpsk_deciding(RecSymbol_temp,'QPSK');
        
        %Demodulate
		Demodulate_data=data_demodulate(RecSymbol_last,'QPSK');
        
		errs=0;
		errs=errs+length(find(Demodulate_data~=TxBits_temp));
		errs_rate(:,nEN)=errs/length(TxBits_temp);   
	end
    
    figure(1)
    if frequcncy_offset_temp==1
	   subplot(221);
	   plot(real(RecSymbol_temp),imag(RecSymbol_temp),'.');
       pictuer_style='-o';
    elseif frequcncy_offset_temp==2
	   subplot(222);
	   plot(real(RecSymbol_temp),imag(RecSymbol_temp),'.');
       pictuer_style='-^';
	elseif frequcncy_offset_temp==3
	   subplot(223);
	   plot(real(RecSymbol_temp),imag(RecSymbol_temp),'.');
       pictuer_style='-*';
	elseif frequcncy_offset_temp==4
	   subplot(224);
	   plot(real(RecSymbol_temp),imag(RecSymbol_temp),'.');
       pictuer_style='-x';
    end
	figure(2)
	semilogy(EbN0db,errs_rate,pictuer_style);
    grid on;
    hold on;    
end
