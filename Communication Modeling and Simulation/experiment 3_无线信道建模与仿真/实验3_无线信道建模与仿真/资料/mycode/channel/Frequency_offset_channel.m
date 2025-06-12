function  [Output_signal]=multipath(original_signal,fT)

N=length(original_signal);

Output_signal=original_signal.*exp(sqrt(-1)*2*pi*fT*(1:N)/N);
