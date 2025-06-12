% Simulation of《Symbol time offset estimation in choerent OFDM systems》
close all;
clear;
clc;
tic;
fc = 2e9;			%Carrier frequency
fs = 20e+6; 		%System bandwith is 5MHz
V=0;
FrameNum=10;
SubCarrNum = 128;	%Number of carriers
GuardNum = 16;      %Number of cyclic prefix

Mf=8;
PilotCarrNum = ceil( SubCarrNum/Mf );
PilotSymbol = (1+j)/sqrt(2);
PilotIndex = 1:Mf:SubCarrNum;
TmpIndex = zeros(1,SubCarrNum);
TmpIndex(PilotIndex) = PilotIndex;
DataIndex = find((1:SubCarrNum) ~= TmpIndex);
clear('TmpIndex');
Reference_sinal=ones(1,SubCarrNum);
Reference_sinal(DataIndex) = zeros(1,length(SubCarrNum));
Reference_sinal(PilotIndex) = Reference_sinal(PilotIndex)*PilotSymbol;
% IFFT of pilot datas
M_sinal_temp = sqrt(SubCarrNum) * ifft(Reference_sinal,SubCarrNum);
%Add the cyclic prefix to the pilot symbol
M_sinal=[M_sinal_temp(:,SubCarrNum-GuardNum+1:SubCarrNum) M_sinal_temp];

for multipath_yn=0:1:1

	ModType=2;          %For QPSK
	Info_length=length(DataIndex) * ModType;
        
	Error_ratio_time=zeros(1,5);
	for SNR_number=1:1:5
        
        SNR=2*SNR_number;
		Rho=SNR/(SNR+1);
		snr=10^(SNR/10);
		E=1;                % 符号能量假设为单位能量.
		Noise_Var = E/(3*ModType*snr);
				
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
            
            %Add multipath
            if multipath_yn==1
		       ChannelInput_temp=Standard_Multipath1(ChannelInput_temp,fc,V,fs);
            end
			
			if nframe==1 
               ChannelInput=ChannelInput_temp;
			else
               ChannelInput=[ChannelInput ChannelInput_temp];    
			end
		end 
		
		
		%Pass through channel
		ChannelOutput=AGWN_channel(ChannelInput,Noise_Var);
		
		
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
		
		error_ratio_time=0;
		for i=2:FrameNum-1
			argmax=max(argmax_before(1,(i*(SubCarrNum+GuardNum)-GuardNum-100):1:(i*(SubCarrNum+GuardNum)-GuardNum+100)));
			for jj=(i*(SubCarrNum+GuardNum)-GuardNum-100):1:(i*(SubCarrNum+GuardNum)-GuardNum+100)
                if argmax_before(1,jj)==argmax 
                   time_offset(1,i-1)=jj-1-i*(SubCarrNum+GuardNum)+GuardNum;
                   if time_offset(1,i-1)==0
                   else
                      error_ratio_time=error_ratio_time+1;
                   end
                end
			end
		end
		
		Error_ratio_time(1,SNR_number)=error_ratio_time/(FrameNum-2);
        
	end
	
    if multipath_yn==1
       picture_style='--*';
    else
       picture_style='-o';
    end
    
	figure(1)
	semilogy(1:1:5,Error_ratio_time,picture_style);
	grid on;
	hold on;

end

%Output figures
figure(6)
subplot(311);
for nframe = 2:FrameNum-1
    plot(0:(SubCarrNum+GuardNum-1),argmax_before_CP(1,(nframe-1)*(SubCarrNum+GuardNum)+1:nframe*(SubCarrNum+GuardNum)),'w');  
    hold on
end
grid on;

subplot(312);
for nframe = 2:FrameNum-1
    plot(0:(SubCarrNum+GuardNum-1),argmax_before_P(1,(nframe-1)*(SubCarrNum+GuardNum)+1:nframe*(SubCarrNum+GuardNum)),'w');  
    hold on
end
grid on;

subplot(313);
for nframe = 2:FrameNum-1
    plot(0:(SubCarrNum+GuardNum-1),argmax_before(1,(nframe-1)*(SubCarrNum+GuardNum)+1:nframe*(SubCarrNum+GuardNum)),'w');  
    hold on
end
grid on;

TimeElapse = toc/60;	% Simulation minutes
fprintf('****** The total time used is : %3.4f minutes ********\n', TimeElapse);