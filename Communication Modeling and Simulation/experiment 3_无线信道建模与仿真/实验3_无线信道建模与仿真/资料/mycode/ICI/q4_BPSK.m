clear;
clc;
fT=3/10;
fc = 2.4e9;			%Carrier frequency
fs = 5e+6; 			%System bandwith is 5MHz
V=88;
EbN0db = 0:1:5;		% Eb/N0
FrameNum=10;
SubCarrNum = 128;	%size fourier transform to generate signal. It
					%is equal to the number of samples in the OFDM
					%symbol, and also it is the size of FFT.
GuardNum = 8;	%Total guard time in samples.
ModType=1;      %for QPSK
Info_length=SubCarrNum*ModType;
errs_rate=zeros(1,length(EbN0db));

for nEN = 1:length(EbN0db)
    E=1;     % 符号能量假设为单位能量.
    snr = 10^(EbN0db(nEN)/10);% en为bit信噪比；ModType*en为符号信噪比.
	Noise_Var = E/(2*ModType*snr);
    
	for nframe = 1:FrameNum
	
        TxBits = randint(1,Info_length,2);

        ModData = data_modulate(TxBits, 'BPSK');
        
        %ICI coding
        ICI_signal=zeros(1,4*length(ModData));
        for i=1:length(ModData)
            ICI_signal(:,(i-1)*4+1)=ModData(:,i);
            ICI_signal(:,(i-1)*4+2)=-3*ModData(:,i);
            ICI_signal(:,(i-1)*4+3)=3*ModData(:,i);
            ICI_signal(:,(i-1)*4+4)=-1*ModData(:,i);
        end
        
        figure(4);
        subplot(221);
        plot(real(ModData),imag(ModData),'.');
        
        % IFFT
        OfdmSymbol = sqrt(SubCarrNum*4) * ifft(ICI_signal,SubCarrNum*4);

        % Add the cyclic prefix to the ofdm symbol
        ChannelInput = [OfdmSymbol(:,SubCarrNum-GuardNum+1:SubCarrNum) OfdmSymbol];
        
        if nframe==1 
           TxBits_temp=TxBits;
        else
           TxBits_temp=[TxBits_temp TxBits];
        end

    
        %Pass through channel
	    ChannelOutput=Improved_Jakes_Multipath2(ChannelInput,fT);
        ChannelOutput=AGWN_channel(ChannelOutput,Noise_Var);

	    RecTimeSymbol = ChannelOutput(:,GuardNum+1:end);
        
        % The FFT of the time domain signal after the removal of cyclic prefix
        RecSymbol= (1/sqrt(SubCarrNum*4)) * fft(RecTimeSymbol,SubCarrNum*4);

        %ICI decoding
        deICI_signal=zeros(1,length(RecSymbol)/4);
        for i=1:length(RecSymbol)/4
            deICI_signal(:,i)=(RecSymbol(:,(i-1)*4+1)+RecSymbol(:,(i-1)*4+2)/(-3)+RecSymbol(:,(i-1)*4+3)/3-RecSymbol(:,(i-1)*4+4))/4;
        end
        
        phase_adjustment=exp(-1*sqrt(-1)*pi*fT*(length(ChannelInput)-1)/length(ChannelInput));
        if nframe==1 
           RecSymbol_temp=phase_adjustment*deICI_signal;
        else
           RecSymbol_temp=[RecSymbol_temp phase_adjustment*deICI_signal];
        end
	end 
	
	subplot(222);
	plot(real(RecSymbol_temp),imag(RecSymbol_temp),'.');
    
	%Decide
	RecSymbol_last=mpsk_deciding(RecSymbol_temp,'BPSK');
	
	subplot(223);
	plot(real(RecSymbol_last),imag(RecSymbol_last),'.');
    
    %Demodulate
	Demodulate_data=data_demodulate(RecSymbol_last,'BPSK');
    
	errs=0;
	errs=errs+length(find(Demodulate_data~=TxBits_temp));
	errs_rate(:,nEN)=errs/length(TxBits_temp);
 
end
figure(5)
% subplot(224)
hold on;
semilogy(EbN0db,errs_rate,'-o');