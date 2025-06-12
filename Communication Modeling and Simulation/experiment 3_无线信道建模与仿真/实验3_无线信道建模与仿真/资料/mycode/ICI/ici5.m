clear all;
clc;
N=64;
Z_temp=zeros(N,N/2);
Z=zeros(N/2,N/2);
zz=zeros(N,N);
for FT=1:9
    fT=0.1*FT;
    interference_power=0;
    signal_power=0;
    
    for m=0:1:(N-1);
        for l=0:1:(N-1)
            zz(m+1,l+1)=(sin(pi*((l-m)+fT))*exp(i*pi*(N-1)*((l-m)+fT)/N))/(N*sin(pi*((l-m)+fT)/N));
        end
    end
    
    for m=1:1:N;
        for l=1:1:N/2
            Z_temp(m,l)=zz(m,(l*2-1))-zz(m,l*2);
        end
    end
    
    for m=1:1:N/2;
        for l=1:1:N/2;
            Z(m,l)=Z_temp((m*2-1),l)-Z_temp(m*2,l);
            if l==m
               signal_power=abs(Z(m,l))^2+signal_power;
            else
               interference_power=abs(Z(m,l))^2+interference_power;
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
