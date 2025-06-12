clear all;
clc;
N=64;
for FT=1:9
    fT=0.1*FT;
    interference_power=0;
    signal_power=0;
    for m=0:1:(N-1);
        for l=0:1:(N-1)
            c1=(sin(pi*((l-m)+fT))*exp(i*pi*(N-1)*((l-m)+fT)/N))/(N*sin(pi*((l-m)+fT)/N));
            c=abs(c1)^2;
            if l~=m
               interference_power=c+interference_power;
            else
               signal_power=c+signal_power;
            end 
        end
    end
    SIN(FT)=10*log10(signal_power/interference_power);
end;
FT=0.1:0.1:0.9;
plot(FT,SIN(FT*10),'-',FT,SIN(FT*10),'*');
grid on;
xlabel('\DeltafT');
ylabel('SIN(dB)');
