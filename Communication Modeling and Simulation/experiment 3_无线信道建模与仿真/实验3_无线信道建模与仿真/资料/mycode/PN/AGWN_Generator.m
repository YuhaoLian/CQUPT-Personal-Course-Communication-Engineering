function  [Output_signal]=AGWN_Generator(Input_signal_number)

Output_signal = (randn(1,Input_signal_number) + j*randn(1,Input_signal_number))/sqrt(2);

return  