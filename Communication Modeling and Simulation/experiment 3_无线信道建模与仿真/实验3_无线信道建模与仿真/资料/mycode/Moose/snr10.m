%by P.Moose
clear;
clc;
EbN0db = 10;			% Eb/N0
FrameNum=1000;
SubCarrNum = 128;	%size fourier transform to generate signal. It
					%is equal to the number of samples in the OFDM
					%symbol, and also it is the size of FFT.
GuardNum = 8;	%Total guard time in samples.
ModType=2;      %for QPSK
Info_length=SubCarrNum*ModType;

est_MSE_error=zeros(1,10);

for fT=1:10;
    
    epsilon=0.1*fT;
    
    E=1;     % 符号能量假设为单位能量.
    snr = 10^(EbN0db/10);% en为bit信噪比；ModType*en为符号信噪比.
	Noise_Var = E/(2*ModType*snr);
    
    est_error=0;
    
	for nframe = 1:FrameNum
	
        TxBits = randint(1,Info_length,2);

        ModData = data_modulate(TxBits, 'QPSK');
        
        % IFFT
        OfdmSymbol = sqrt(SubCarrNum) * ifft(ModData,SubCarrNum);
        
        OfdmSymbol=[OfdmSymbol OfdmSymbol];
        
        OfdmSymbol= Frequency_offset_channel(OfdmSymbol,epsilon);
        
        % Add the cyclic prefix to the ofdm symbol
        ChannelInput = [OfdmSymbol(:,SubCarrNum*2-GuardNum+1:SubCarrNum*2) OfdmSymbol];
    
        %Pass through channel
	    ChannelOutput=AGWN_channel(ChannelInput,Noise_Var);

	    RecTimeSymbol = ChannelOutput(:,GuardNum+1:end);
        
        xxx1= RecTimeSymbol(1,1:SubCarrNum);
        % The FFT of the time domain signal after the removal of cyclic prefix
        xxx1= (1/sqrt(SubCarrNum)) * fft(xxx1,SubCarrNum);
        
        xxx2= RecTimeSymbol(1,1+SubCarrNum:end);
        % The FFT of the time domain signal after the removal of cyclic prefix
        xxx2= (1/sqrt(SubCarrNum)) * fft(xxx2,SubCarrNum);
        
        yyy=xxx2.*conj(xxx1);
        epsilon_est=atan(sum(imag(yyy))/sum(real(yyy)))/(pi);
        est_error=(epsilon_est-epsilon)^2+est_error;

	end 
    
    est_MSE_error(1,fT)=est_error/FrameNum;
    
end

figure(1)
hold on;
semilogy((1:10)/10,est_MSE_error,'-s');