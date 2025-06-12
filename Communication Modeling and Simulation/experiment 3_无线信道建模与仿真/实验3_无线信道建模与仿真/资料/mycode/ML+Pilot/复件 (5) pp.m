% Simulation of《Symbol time offset estimation in choerent OFDM systems》
close all;
clear;
clc;
tic;
fc = 2.4e9;			%Carrier frequency
fs = 5e+6; 			%System bandwith is 5MHz
V=5;
FrameNum=100;
SubCarrNum = 128;	%Number of carriers
GuardNum = 8;       %Number of cyclic prefix

MF_container=[4 8 16 32 64];
for Mf_id = 1:2:length(MF_container)
    
    %Pilot information
    Mf=MF_container(1,Mf_id);
	PilotCarrNum = ceil( SubCarrNum/Mf );
	PilotSymbol = (1+j)/sqrt(2);
	PilotIndex = 1:Mf:SubCarrNum;
	TmpIndex = zeros(1,SubCarrNum);
	TmpIndex(PilotIndex) = PilotIndex;
	DataIndex = find((1:SubCarrNum) ~= TmpIndex);
	clear('TmpIndex');
	
	ModType=2;          %For QPSK
	Info_length=length(DataIndex) * ModType;
        
	%SNR=1;
	for SNR=0:5
		Rho=SNR/(SNR+1);
		snr=10^(SNR/10);
		E=1;                % 符号能量假设为单位能量.
		Noise_Var = E/(2*ModType*snr);
		
		% est_angle=zeros(1,17);
		% est_time=zeros(1,17);
		
		Reference_sinal=ones(1,SubCarrNum);
		Reference_sinal(DataIndex) = zeros(1,length(SubCarrNum));
		Reference_sinal(PilotIndex) = Reference_sinal(PilotIndex)*PilotSymbol;
		% IFFT of pilot datas
		M_sinal_temp = sqrt(SubCarrNum) * ifft(Reference_sinal,SubCarrNum);
		%Add the cyclic prefix to the pilot symbol
		M_sinal=[M_sinal_temp(:,SubCarrNum-GuardNum+1:SubCarrNum) M_sinal_temp];
		
		for nframe = 1:FrameNum
		
			TxBits = randint(1,Info_length,2);
			
			ModData = data_modulate(TxBits, 'QPSK');
			
			% Insert pilot data
			TransSymbol = ones(1,SubCarrNum);
			TransSymbol(DataIndex) = ModData;
			TransSymbol(PilotIndex) = TransSymbol(PilotIndex)*PilotSymbol;
			
			% IFFT
			OfdmSymbol = sqrt(SubCarrNum) * ifft(TransSymbol,SubCarrNum);
			
			% Add the cyclic prefix to the ofdm symbol
			ChannelInput_temp = [OfdmSymbol(:,SubCarrNum-GuardNum+1:SubCarrNum) OfdmSymbol];
			
			if nframe==1 
               ChannelInput=ChannelInput_temp;
			else
               ChannelInput=[ChannelInput ChannelInput_temp];    
			end
		end 
		
		
		%Pass through channel
		ChannelOutput=AGWN_channel(ChannelInput,Noise_Var);
		%ChannelOutput_temp=Standard_Multipath1( ChannelInput,fc,V,fs);
		
		
		%信号移位相乘
		ChannelOutput_delayN=[zeros(1,SubCarrNum) ChannelOutput(:,1:end-SubCarrNum)];
		before_Feim_CP=ChannelOutput_delayN.*conj(ChannelOutput);
		
		for i=1:length(before_Feim_CP)-GuardNum
			Gamam_PP_temp=0;
			Feim_CP_temp=0;
			Feim_P_temp=0;
			for j=1:GuardNum
                Gamam_PP_temp=Gamam_PP_temp+before_Feim_CP(:,(i+j-1));
                Feim_CP_temp=Feim_CP_temp+(abs(ChannelOutput(:,(i+j-1))))^2+(abs(ChannelOutput_delayN(:,(i+j-1))))^2;
                Feim_P_temp=Feim_P_temp+conj(ChannelOutput(:,(i+j-1))+ChannelOutput_delayN(:,(i+j-1)))*M_sinal(1,j);
			end
			Gamam_PP(1,i)=Gamam_PP_temp;
			Feim_CP(1,i)=Feim_CP_temp;
			Feim_P(1,i)=Feim_P_temp;
		end
			
		ChannelOutput2=[zeros(1,length(M_sinal)) ChannelOutput];
			
		for ii=1:length(ChannelOutput2)-length(M_sinal);
			Gamam_P_temp=0;
			for jj=1:length(M_sinal)
                Gamam_P_temp=Gamam_P_temp+M_sinal(1,jj)*conj(ChannelOutput2(1,jj+ii-1));
			end
			Gamam_P(1,ii)=Gamam_P_temp;
		end
		
		Gamam_P=Gamam_P(1,GuardNum+1:end);
		
		argmax_before_CP=real(Gamam_PP)-(Rho*Feim_CP)/2;
		argmax_before_P=(1+Rho)*real(Gamam_P)-Rho*real(Feim_P);
		argmax_before=(1-Rho)*argmax_before_P+Rho*argmax_before_CP;
		
		error_time=0;
		error_ratio_time=0;
		for i=2:FrameNum-1
			argmax=max(argmax_before(1,(i*(SubCarrNum+GuardNum)-GuardNum-100):1:(i*(SubCarrNum+GuardNum)-GuardNum+100)));
			for jj=(i*(SubCarrNum+GuardNum)-GuardNum-50):1:(i*(SubCarrNum+GuardNum)-GuardNum+50)
                if argmax_before(1,jj)==argmax 
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
		Error_ratio_time(1,SNR+1)=error_ratio_time/(FrameNum-2);
	end
	
    if Mf_id==1
       picture_style='-o';
    elseif Mf_id==2
       picture_style='-^';
    elseif Mf_id==3
       picture_style='-v';
    elseif Mf_id==4
       picture_style='-*';
    elseif Mf_id==5
       picture_style='-x';
    end
    
	figure(1)
	semilogy(0:1:5,Error_ratio_time,picture_style);
	grid on;
	hold on;
   
end

%Output figures
figure(2)
correlation_flag=length(M_sinal); 
[R_xx,lags] = xcorr(M_sinal,correlation_flag,'coeff');
subplot(211);
plot(1:length(M_sinal),real(M_sinal));
grid on;

subplot(212);
plot(lags,real(R_xx))
grid on;

figure(3)
subplot(211);
plot(1:length(argmax_before_CP),argmax_before_CP);
grid on;

subplot(212);
for nframe = 1:FrameNum-1
    plot(0:(SubCarrNum+GuardNum-1),argmax_before_CP(1,(nframe-1)*(SubCarrNum+GuardNum)+1:nframe*(SubCarrNum+GuardNum)));  
    hold on
end
grid on;

figure(4)
subplot(211);
plot(1:length(argmax_before_P),argmax_before_P);
grid on;

subplot(212);
for nframe = 1:FrameNum-1
    plot(0:(SubCarrNum+GuardNum-1),argmax_before_P(1,(nframe-1)*(SubCarrNum+GuardNum)+1:nframe*(SubCarrNum+GuardNum)));  
    hold on
end
grid on;

figure(5)
subplot(211);
plot(1:length(argmax_before),argmax_before);
grid on;

subplot(212);
for nframe = 1:FrameNum-1
    plot(0:(SubCarrNum+GuardNum-1),argmax_before(1,(nframe-1)*(SubCarrNum+GuardNum)+1:nframe*(SubCarrNum+GuardNum)));  
    hold on
end
grid on;

figure(6)
subplot(311);
for nframe = 1:FrameNum-1
    plot(0:(SubCarrNum+GuardNum-1),argmax_before_CP(1,(nframe-1)*(SubCarrNum+GuardNum)+1:nframe*(SubCarrNum+GuardNum)));  
    hold on
end
grid on;

subplot(312);
for nframe = 1:FrameNum-1
    plot(0:(SubCarrNum+GuardNum-1),argmax_before_P(1,(nframe-1)*(SubCarrNum+GuardNum)+1:nframe*(SubCarrNum+GuardNum)));  
    hold on
end
grid on;

subplot(313);
for nframe = 1:FrameNum-1
    plot(0:(SubCarrNum+GuardNum-1),argmax_before(1,(nframe-1)*(SubCarrNum+GuardNum)+1:nframe*(SubCarrNum+GuardNum)));  
    hold on
end
grid on;

TimeElapse = toc/60;	% Simulation minutes
fprintf('****** The total time used is : %3.4f minutes ********\n', TimeElapse);