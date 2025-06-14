clear;
clc;
fc = 2.4e9;			%Carrier frequency
fs = 5e+6; 			%System bandwith is 5MHz
V=88;
EbN0db = 0;			% Eb/N0
FrameNum=1000;
SubCarrNum = 128;	%size fourier transform to generate signal. It
					%is equal to the number of samples in the OFDM
					%symbol, and also it is the size of FFT.
GuardNum = 8;	%Total guard time in samples.
ModType=2;      %for QPSK
Info_length=SubCarrNum*ModType;
errs_rate=zeros(1,length(EbN0db));

for fT=1:10;
    
    epsilon=0.1*fT;
    
    E=1;     % 符号能量假设为单位能量.
    snr = 10^(EbN0db/10);% en为bit信噪比；ModType*en为符号信噪比.
	Noise_Var = E/(2*ModType*snr);
    
	for nframe = 1:FrameNum
	
        TxBits = randint(1,Info_length,2);

        ModData = data_modulate(TxBits, 'QPSK');
        
        %ICI coding
        ICI_signal=zeros(1,3*length(ModData));
        for i=1:length(ModData)
            ICI_signal(:,(i-1)*3+1)=ModData(:,i);
            ICI_signal(:,(i-1)*3+2)=-2*ModData(:,i);
            ICI_signal(:,(i-1)*3+3)=ModData(:,i);
        end
        
%         figure(3);
%         subplot(221);
%         plot(real(ModData),imag(ModData),'.');
        
        % IFFT
        OfdmSymbol = sqrt(SubCarrNum*3) * ifft(ICI_signal,SubCarrNum*3);

        % Add the cyclic prefix to the ofdm symbol
        ChannelInput = [OfdmSymbol(:,SubCarrNum-GuardNum+1:SubCarrNum) OfdmSymbol];
        
        if nframe==1 
           TxBits_temp=TxBits;
        else
           TxBits_temp=[TxBits_temp TxBits];
        end

    
        %Pass through channel
	    ChannelOutput=Improved_Jakes_Multipath2(ChannelInput,epsilon);
        ChannelOutput=AGWN_channel(ChannelOutput,Noise_Var);

	    RecTimeSymbol = ChannelOutput(:,GuardNum+1:end);
        
        % The FFT of the time domain signal after the removal of cyclic prefix
        RecSymbol= (1/sqrt(SubCarrNum*3)) * fft(RecTimeSymbol,SubCarrNum*3);

        %ICI decoding
        deICI_signal=zeros(1,length(RecSymbol)/3);
        for i=1:length(RecSymbol)/3
            deICI_signal(:,i)=(RecSymbol(:,(i-1)*3+1)-RecSymbol(:,(i-1)*3+2)/2+RecSymbol(:,(i-1)*3+3))/3;
        end
        
        phase_adjustment=exp(-1*sqrt(-1)*pi*epsilon*(length(ChannelInput)-1)/length(ChannelInput));
        if nframe==1 
           RecSymbol_temp=phase_adjustment*deICI_signal;
        else
           RecSymbol_temp=[RecSymbol_temp phase_adjustment*deICI_signal];
        end
	end 
	
% 	subplot(222);
% 	plot(real(RecSymbol_temp),imag(RecSymbol_temp),'.');
    
	%Decide
	RecSymbol_last=mpsk_deciding(RecSymbol_temp,'QPSK');
	
% 	subplot(223);
% 	plot(real(RecSymbol_last),imag(RecSymbol_last),'.');
    
    %Demodulate
	Demodulate_data=data_demodulate(RecSymbol_last,'QPSK');
    
	errs=0;
	errs=errs+length(find(Demodulate_data~=TxBits_temp));
	errs_rate(:,fT)=errs/length(TxBits_temp);
 
end
figure(2)
% subplot(224)
hold on;
semilogy(1:10,errs_rate,'-s');