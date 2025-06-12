function [OutSignal] = MultiPathTV(TimeSignal,fc,v,fs)
% This program acts to generate a Tapped-Delay-Line(TDL) channel coefficients
% via the modified Jakes's model, for which the first and second statistics
% of ndividual path waveforms agree well with the theoretical expectations.
%
%	[OutSignal] = MultiPath(fc, v, fs, t0, Duration, PowerDelayProfile, NonZeroPaths, TimeSignal)
%
%
%	Parameters:
%
%	Input:
%			fc: The carrier frequency for the simulated mobile channel.
%			v: The moving velocity of the receiver, where the transmitter is assumed to be fixed in location.
%			fs: The sampling frequency of the generated multipath channel.
%			t0: Initial time for channel simulation process.
%			Duration: The sampled channel responses number.
%			PowerDelayProfile: The power delay profile in dB for the Tapped-Delay-Line Rayleigh fading channel.
%			NonZeroPaths: Locations for nonzero paths.
%			TimeSignal	: Time waveform of signal to apply the channel to.
%			NoiseVar	: AWGN noise variance
%	Output:
%			OutSignal	: Signal via multipath and AWGN channel.

PowerDelayProfile = [0,-4,-9,-10,-15,-20];
NonZeroPaths = [0,1,2,3,4,6,9,13];

j = sqrt(-1);

tapNum = length(NonZeroPaths);

Ts = 1/fs;
N0 = 24;
N  = 4*N0;
M  = tapNum;
Ank= 2*pi/(N*M)*((0:N0-1)'*ones(1,M)*M+ones(N0,1)*(0:M-1)+1/4);
Wm = (2*pi*v*fc/3.e+8)/3.6;
Phase = 2*pi*[0.8 0.28];							% Random Phase 
TimeBegin=0;

for kk = TimeBegin + 1 : TimeBegin + length(TimeSignal)
	TDLCoef(:,kk-TimeBegin)=sqrt(1/N0)*sum(cos(Wm*cos(Ank)*(kk+1e8)*Ts+Phase(1))+j*sin(Wm*sin(Ank)*(kk+1e8)*Ts+Phase(2))).';
end

TempSig = zeros(tapNum,size(TDLCoef,2));
for i = 1 : tapNum
	if NonZeroPaths(1,i) > 0
		TempSig(i,:) = [zeros(1,NonZeroPaths(1,i)) TimeSignal(1,[1:end-NonZeroPaths(1,i)])];
	else
		TempSig(i,:) = TimeSignal;
	end
end

%Mean amplitude values for different paths in decimal format:
MeanAmp = exp(-1*NonZeroPaths);%10.^(PowerDelayProfile/10);
MeanAmp = sqrt( MeanAmp/sum(MeanAmp) );	%Nomorlized the power to 1.
OutSignal = MeanAmp*(TempSig.*TDLCoef)/(length(NonZeroPaths)); %pass through the multipath channel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   Add Gausian White Noise %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noise = sqrt(NoiseVar)*1/sqrt(2)*(randn(1,length(TimeSignal)) + j*randn(1,length(TimeSignal)));
%RANDN(N) is an N-by-N matrix with random entries, chosen from a normal distribution with mean zero, variance one and standard deviation one.
%RANDN(M,N) and RANDN([M,N]) are M-by-N matrices with random entries.
% OutSignal = OutSignal + noise;
