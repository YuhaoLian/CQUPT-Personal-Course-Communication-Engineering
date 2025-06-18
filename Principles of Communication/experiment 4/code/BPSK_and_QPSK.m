clc
clear
close all;

%% 初始化参数
data_len = 100000; % 二进制数据长度为10^5
SNR_dB = -2:10; % SNR信噪比dB值
SNR = 10.^(SNR_dB/10); % Eb/N0
Eb = 1; % 每比特能量
N0 = Eb./SNR ; % 噪声功率
iter = 500; %迭代次数
error2 = zeros(iter,length(SNR_dB));          % 错误符号个数
simu_ber_BPSK = zeros(1,length(SNR_dB));         % 仿真BPSK误比特率
theory_ber_BPSK = zeros(1,length(SNR_dB));   % 理论BPSK误比特率
demod2_signal= zeros(1,data_len);         % 解调信号

%% 基带信号的产生
data_source = round(rand(1,data_len));  % 二进制随机序列

%% BPSK 频带调制  
send_signal2 = (data_source - 1/2)*2; 

%% 高斯信道白噪声
for iter = 1:iter
    for z = 1:length(SNR_dB)
         noise2 = sqrt(N0(z)/2) * randn(1,data_len); % 白噪声
         receive_signal2 = send_signal2 + noise2;
         demod_signal2 = zeros(1,data_len);     %解调信号的产生
         % 解调与检测
            for w = 1:data_len
                    if (receive_signal2(w) > 0)
                    demod_signal2(w) = 1;              % 1 if the received signal is greater than 0
                    else
                    demod_signal2(w) = 0;              % 0 if the received signal is less than 0
                    end
            end
            % 计算错误符号的个数
           for w = 1:data_len
               if(demod_signal2(w) ~=data_source(w) )
                      error2(iter,z) = error2(iter,z) + 1;    % Number of error bits
               end
           end
             %计算误比特率
            simu_ber_BPSK(iter,z) = error2(iter,z) / data_len;         % Simulation bit error rate
            theory_ber_BPSK(z) = qfunc(sqrt(2*SNR(z)));   % Theoretical bit error rate
    end
end
%% map  
figure(1);
stem(data_source);
title("Binary random sequence");
axis([0,50,0,1]);
figure(2);
stem(send_signal2);
title("BPSK baseband modulation -- transmit signal");
axis([0,50,-1.5,1.5]);

figure(4);
stem(noise2);
title("white Gaussian noise");
axis([0,50,-0.5,0.5]);

figure(5)
stem(receive_signal2);
title("receive signals");
axis([0,50,-1.5,1.5]);

figure(7)
stem(demod_signal2);
title("demodulation signals");
axis([0,50,0,1]);

figure(8);
semilogy(SNR_dB,mean(simu_ber_BPSK,1),'M-X',SNR_dB,theory_ber_BPSK,'k-s');     
grid on;                                      
axis([0 10 10^-5 10^-1])                      
xlabel('Eb/N0 (dB)');                     
ylabel('BER');                                  
 legend('BPSK Simulation bit error rate','BPSK Theoretical bit error rate');  

scatterplot(send_signal2);
title('Transmission signal constellation diagram');
scatterplot(receive_signal2);
title('Receiving signal constellation diagram');
scatterplot(demod_signal2);
title('Decoded signal constellation');

%% QPSK
s00 = 1/sqrt(2)*(-1-1j);
s01 = 1/sqrt(2)*(-1+1j);
s10 = 1/sqrt(2)*(1-1j);
s11 = 1/sqrt(2)*(1+1j);
itermax = 500;  %迭代次数
num = 5000;
snr = 0:10;
snrn = 10.^(snr/10);
ber = zeros(itermax,length(snr));

for iter = 1:itermax;
    s = 1/sqrt(2)*(randsrc(1,num)+1j*randsrc(1,num));
    for i = 1:length(snr)
        N(i) = (10.^(-snr(i)/10));
        n = sqrt(N(i)/2).*(randn(1,num)+1j.*randn(1,num));
        r = s+n;
        
        for k = 1:num
            d00 = abs(r(k)-s00);
            d01 = abs(r(k)-s01);
            d10 = abs(r(k)-s10);
            d11 = abs(r(k)-s11);
            if((d00<d01)&(d00<d10)&(d00<d11))
                y(k) = -1/sqrt(2)-1j/sqrt(2);
            elseif((d01<d00)&(d01<d10)&(d01<d11))
                y(k) = -1/sqrt(2)+1j/sqrt(2);
            elseif((d10<d00)&(d10<d01)&(d10<d11))
                y(k) = 1/sqrt(2)-1j/sqrt(2);
            else
                y(k) = 1/sqrt(2)+1j/sqrt(2);
            end
        end
        ber(iter,i) = length(find(s~=y))/(num);
    end
end
bern = sum(ber)/itermax;
q = 2*qfunc(sqrt(snrn)).*(1-0.5*qfunc(sqrt(snrn)));

figure(12);
semilogy(snr,bern,'M-X',snr,q,'k-s');     
grid on;                                                          
xlabel('Eb/N0 (dB)');                     
ylabel('BER');                                  
legend('QPSK Simulation bit error rate','QPSK Theoretical bit error rate');  

scatterplot(s);
title('Transmission signal constellation diagram');
scatterplot(r);
title('Receiving signal constellation diagram');
scatterplot(y);
title('Decoded signal constellation');
