function  [Output_signal]=AGWN_channel(Input_signal,NoiseVar)

data_length=length(Input_signal);
%Add noise
noise = sqrt(NoiseVar)*(randn(1,data_length) + j*randn(1,data_length))/sqrt(2);
Output_signal=Input_signal+ noise;
return  