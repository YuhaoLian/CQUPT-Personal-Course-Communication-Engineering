% Simulation of《Symbol time offset estimation in choerent OFDM systems》
close all;
clear;
clc;
tic;
fc = 2.4e9;			%Carrier frequency
fs = 5e+6; 			%System bandwith is 5MHz
V=50;
FrameNum=10;
SubCarrNum =128;	%Number of carriers
GuardNum = 16;       %Number of cyclic prefix
info_length=SubCarrNum+GuardNum;
frequency_offset=0.3;  %0.1 to 0.8

Mf = 128;
Rho_power=1+(info_length/Mf);
power_coefficient=500;

%Pilot information
PilotCarrNum = ceil( SubCarrNum/Mf );
PilotSymbol = (sqrt(Rho_power)*power_coefficient*(1+j))/sqrt(2);
PilotIndex = 1:Mf:SubCarrNum;
TmpIndex = zeros(1,SubCarrNum);
TmpIndex(PilotIndex) = PilotIndex;
DataIndex = find((1:SubCarrNum) ~= TmpIndex);
clear('TmpIndex');

TransSymbol_Pilot = ones(1,SubCarrNum);
TransSymbol_Pilot(DataIndex) = zeros(1,length(SubCarrNum));
TransSymbol_Pilot(PilotIndex) = TransSymbol_Pilot(PilotIndex)*PilotSymbol;

% IFFT of pilot datas
OfdmSymbol_Pilot = sqrt(SubCarrNum) * ifft(TransSymbol_Pilot,SubCarrNum);
        
PN_code=PN_Signal(length(find(OfdmSymbol_Pilot(1:end) ~= 0)));
%     PN_code=AGWN_Generator(length(find(OfdmSymbol_Pilot(1:end) ~= 0)));
    
jj=1;
for ii=1:length(OfdmSymbol_Pilot)
    if OfdmSymbol_Pilot(1,ii) ~= 0
       OfdmSymbol_Pilot_I(1,ii)=real(OfdmSymbol_Pilot(1,ii))*real(PN_code(1,jj));
       OfdmSymbol_Pilot_Q(1,ii)=imag(OfdmSymbol_Pilot(1,ii))*imag(PN_code(1,jj));
       OfdmSymbol_Pilot(1,ii)=OfdmSymbol_Pilot_I(1,ii)+sqrt(-1)*OfdmSymbol_Pilot_Q(1,ii);
       jj=1+jj;
    end       
end
   
%Add the cyclic prefix to the pilot symbol
ChannelInput_Pilot = [OfdmSymbol_Pilot(:,SubCarrNum-GuardNum+1:SubCarrNum) OfdmSymbol_Pilot];

ModType=2;          %For QPSK
Info_length=length(DataIndex) * ModType;

kk=0;    
for SNR=0:10:10
    kk=kk+1;
    Rho=SNR/(SNR+1);
	snr=10^(SNR/10);
	E=1*power_coefficient; % 符号能量假设为单位能量(跟据导频调整).
	Noise_Var = E/(2*ModType*snr);
	       
	for nframe = 1:FrameNum
	
		TxBits = randint(1,Info_length,2);
		
		ModData = data_modulate(TxBits, 'QPSK');
		
		% Insert pilot data
		TransSymbol_Data = zeros(1,SubCarrNum);
		TransSymbol_Data(DataIndex) = ModData;
        			
		% IFFT
		OfdmSymbol_Data = sqrt(SubCarrNum) * ifft(TransSymbol_Data,SubCarrNum);
        
		% Add the cyclic prefix to the ofdm symbol
		ChannelInput_Data = [OfdmSymbol_Data(:,SubCarrNum-GuardNum+1:SubCarrNum) OfdmSymbol_Data];
		                                 
        ChannelInput_temp=ChannelInput_Pilot + ChannelInput_Data;
        ChannelInput_temp=Frequency_offset_channel(ChannelInput_temp,frequency_offset);
        
		if nframe==1 
           ChannelInput=ChannelInput_temp;
		else
           ChannelInput=[ChannelInput ChannelInput_temp];    
		end
	end 
	
	
	%Pass through channel
	ChannelOutput=AGWN_channel(ChannelInput,Noise_Var);
    ChannelOutput=multipath1(ChannelInput,2e9,50,1/2e6);
% 	ChannelOutput=Multipath(ChannelOutput,fc,V,fs);
	
	
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
            Feim_P_temp=Feim_P_temp+conj(ChannelOutput(:,(i+j-1))+ChannelOutput_delayN(:,(i+j-1)))*ChannelInput_Pilot(1,j);
		end
		Gamam_PP(1,i)=Gamam_PP_temp;
		Feim_CP(1,i)=Feim_CP_temp;
		Feim_P(1,i)=Feim_P_temp;
	end
		
	ChannelOutput2=[zeros(1,length(ChannelInput_Pilot)) ChannelOutput];
		
	for ii=1:length(ChannelOutput2)-length(ChannelInput_Pilot);
		Gamam_P_temp=0;
		for jj=1:length(ChannelInput_Pilot)
            Gamam_P_temp=Gamam_P_temp+conj(ChannelInput_Pilot(1,jj))*ChannelOutput2(1,jj+ii-1);
		end
		Gamam_P(1,ii)=Gamam_P_temp;
	end
	
	Gamam_P=Gamam_P(1,GuardNum+1:end);
	
	argmax_before_CP=abs(Gamam_PP)-(Rho*Feim_CP)/2;
	argmax_before_P=(1+Rho)*abs(Gamam_P)-Rho*real(Feim_P);
	argmax_before=(1-Rho)*argmax_before_P+Rho*argmax_before_CP;
    
	error_time=0;
	for i=2:FrameNum-1
		argmax=max(argmax_before(1,(i*(info_length)-GuardNum-100):1:(i*(info_length)-GuardNum+100)));
        for jj=(i*(info_length)-GuardNum-100):1:(i*(info_length)-GuardNum+100)
            if argmax_before(1,jj)==argmax;
               angle_offset(1,i-1)=angle(Gamam_P(jj))/pi;
               time_offset(1,i-1)=jj-1-i*(info_length)+GuardNum;
               error_time=abs(time_offset(1,i-1))^2+error_time;
           end
		end
	end
	
	MSE_time(1,kk)=sqrt(error_time/(FrameNum-2))
    angle_offset_average=sum(angle_offset)/length(angle_offset);
    sum_angle=0;
    for ii=1:length(angle_offset)
        sum_angle=(angle_offset(1,ii)-angle_offset_average)^2+sum_angle;
    end
    MSE_angle(1,kk)=sqrt(sum_angle)/length(angle_offset)
end

%Output figures
figure(1)
semilogy(1:1:kk,MSE_time,'-o');
grid on;
hold on;

figure(2)
semilogy(1:1:kk,MSE_angle,'-o');
grid on;
hold on;

figure(3)
correlation_flag=length(PN_code); 
[R_xx,lags] = xcorr(PN_code,correlation_flag,'coeff');
subplot(211);
plot(1:length(PN_code),real(PN_code));
grid on;

subplot(212);
plot(lags,real(R_xx))
grid on;

figure(4)
correlation_flag=length(ChannelInput_Pilot); 
[R_xx,lags] = xcorr(ChannelInput_Pilot,correlation_flag,'coeff');
subplot(211);
plot(1:length(ChannelInput_Pilot),real(ChannelInput_Pilot));
grid on;

subplot(212);
plot(lags,real(R_xx))
grid on;

figure(5)
subplot(211);
plot(1:length(argmax_before_CP),argmax_before_CP);
grid on;

subplot(212);
for nframe = 2:FrameNum-1
    plot(0:(info_length-1),argmax_before_CP(1,(nframe-1)*(info_length)+1:nframe*(info_length)));  
    hold on
end
grid on;

figure(6)
subplot(211);
plot(1:length(argmax_before_P),argmax_before_P);
grid on;

subplot(212);
for nframe = 2:FrameNum-1
    plot(0:(info_length-1),argmax_before_P(1,(nframe-1)*(info_length)+1:nframe*(info_length)));  
    hold on
end
grid on;

figure(7)
subplot(211);
plot(1:length(argmax_before),argmax_before);
grid on;

subplot(212);
for nframe = 2:FrameNum-1
    plot(0:(info_length-1),argmax_before(1,(nframe-1)*(info_length)+1:nframe*(info_length)));  
    hold on
end
grid on;

figure(8)
subplot(311);
for nframe = 2:FrameNum-1
    plot(0:(info_length-1),argmax_before_CP(1,(nframe-1)*(info_length)+1:nframe*(info_length)));  
    hold on
end
grid on;

subplot(312);
for nframe = 2:FrameNum-1
    plot(0:(info_length-1),argmax_before_P(1,(nframe-1)*(info_length)+1:nframe*(info_length)));  
    hold on
end
grid on;

subplot(313);
for nframe = 2:FrameNum-1
    plot(0:(info_length-1),argmax_before(1,(nframe-1)*(info_length)+1:nframe*(info_length)));  
    hold on
end
grid on;


figure(9)
subplot(211)
ChannelInput_Pilot
ChannelInput_Pilot1=multipath1(ChannelInput_Pilot,2e9,50,1/2e6);
plot(1:length(ChannelInput_Pilot1),ChannelInput_Pilot1);
grid on;

subplot(212)
[Crosscorrelation5,lags] = xcorr(ChannelInput_Pilot1,ChannelInput_Pilot,correlation_flag,'coeff');
plot(lags,abs(Crosscorrelation5))
grid on

TimeElapse = toc/60;	% Simulation minutes
fprintf('****** The total time used is : %3.4f minutes ********\n', TimeElapse);