% Simulation of《Symbol time offset estimation in choerent OFDM systems》
% close all;
clear;
clc;
tic;
fc = 2e9;			%Carrier frequency
fs = 20e+6; 		%System bandwith is 20MHz
V=0;
SubCarrNum = 16;	%Number of carriers
FrameNum=10;
M=10;
Info_length=SubCarrNum;
frequency_offset=0.3;
PilotSymbol = (1+j)/sqrt(2);

PN_code=PN_Signal(Info_length*2);
PN_code=PN_code(1,SubCarrNum+1:end);
% PN_code=AGWN_Generator(Info_length);

kkkk=0:10:100;

kk2=0;

for SNR=kkkk
    SNR
    kk2=kk2+1;
    error_time=0;
	error_angle=0;

	snr=10^(SNR/10);
	E=1;                % 符号能量假设为单位能量.
	Noise_Var = E/(2*snr);
    for nFrame=1:FrameNum
        
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
        
        ChannelInput=Frequency_offset_channel(ChannelInput,frequency_offset);
        
        TxBits = randint(1,Info_length*20,2);		
		ModData = data_modulate(TxBits, 'QPSK');
        OfdmSymbol_Data1 = sqrt(Info_length*20) * ifft(ModData,Info_length*20);
        
        TxBits = randint(1,Info_length*20,2);		
		ModData = data_modulate(TxBits, 'QPSK');
        OfdmSymbol_Data2 = sqrt(Info_length*20) * ifft(ModData,Info_length*20);
        
        ChannelInput=[OfdmSymbol_Data1 ChannelInput OfdmSymbol_Data2];
        
        %Pass through channel
		ChannelOutput=AGWN_channel(ChannelInput,Noise_Var);
	    ChannelOutput=Standard_Multipath1(ChannelOutput,fc,V,fs);
	
		for ii=1:length(ChannelOutput)-Info_length
            Gamam_P_temp=0;
			for jj=1:length(PN_code)
                Gamam_P_temp=Gamam_P_temp+PN_code(1,jj)*conj(ChannelOutput(1,jj+ii-1));
			end
            Gamam_P(1,ii)=Gamam_P_temp;
		end
	
        Gamam_P2=[Gamam_P zeros(1,Info_length)];
        Gamam_P3=[zeros(1,Info_length) Gamam_P];
        
        argmax2=Gamam_P3.*conj(Gamam_P2);
        argmax2=argmax2(1,1+Info_length:end-Info_length);
        
		argmax=abs(Gamam_P);
	
		for ii=1:length(argmax2)-Info_length
            Gamam_P2_temp=0;
			for jj=1:Info_length
                Gamam_P2_temp=Gamam_P2_temp+argmax2(1,jj+ii-1);
            end
            Gamam_P2(1,ii)=Gamam_P2_temp;
		end
        
        Gamam_P2=Gamam_P2/max(Gamam_P2);
        Gamam_P2=abs(Gamam_P2);
        
        for ii=1:length(Gamam_P2)-Info_length*4
            sample_count=0;
            if Gamam_P2(1,ii)>0.5
               for kk=1:Info_length*4
                   if Gamam_P2(1,ii+kk)>0.5
                      sample_count=sample_count+1;
                   end
               end  
            end
            if sample_count==Info_length*4
               break;
            end
         end
         mmm=ii+SubCarrNum*M-17;
         if (ii+mmm+4)<length(argmax2)
             argMax=max(argmax2(1,mmm-4:1:mmm+4));
             for jj=mmm-4:1:mmm+4
                 if argmax2(1,jj)==argMax 
                    break
                 end
			 end
             jj=449;
             angle(argmax2(1,jj-2*Info_length))*M/(pi*2);
             error_angle=(angle(argmax2(1,jj-2*Info_length))*M/(pi*2)-frequency_offset)^2+error_angle;
         
             if (jj-449)==0
             else
                error_time=error_time+1;       
             end
         else
           error_time=error_time+1;
           error_angle=1+error_angle;
         end 
     end
     
     MSE_angle(1,kk2)=sqrt(error_angle/FrameNum);
     RrrorRate_time(1,kk2)=sqrt(error_time/FrameNum);
end

%Output figures

figure(1)
semilogy(kkkk,MSE_angle,'--x');
grid on;
hold on;

figure(2)
semilogy(kkkk,RrrorRate_time,'--x');
grid on;
hold on;

figure(3)
subplot(211);
plot(1:length(argmax2),argmax2);
grid on;

subplot(212);
plot(1:length(Gamam_P2),Gamam_P2);
grid on;

figure(4)
plot(1:length(Gamam_P2),Gamam_P2);
grid on;

TimeElapse = toc/60;	% Simulation minutes
fprintf('****** The total time used is : %3.4f minutes ********\n', TimeElapse);