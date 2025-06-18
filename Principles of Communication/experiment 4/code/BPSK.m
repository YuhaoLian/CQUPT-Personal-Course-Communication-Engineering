%% <<<<<<<<<<<<<<<<<<< BPSK Modulation and Demodulation >>>>>>>>>>>>>>>>>>>
clc;
clear all;
close all;

%% ******************* 初始化 ******************** %%
fc = 8; %载波频率
Rb = fc; %码元速率
T = 1 / Rb; % 一个码元传送的时间
dt = 1/2000; %一个码元的抽样时间
fs = 1/dt; %采样频率
nfft=1024; %采样的点数
N = 100000;%10^5个二进制比特
t=0:dt:N*T-dt;% 所有码元的时间
ta=0:dt:T-dt; % 一个码元时间向量
%此时，一个码元序列被抽样T/dt次

A = 1;%振幅等于1
SNR_dB = -2:10; %信噪比——SNR dB值
SNR = 10.^(SNR_dB./10); %信噪比——SNR
Eb = A^2*T/2; % PSK频带信号每比特能量
N0 = Eb./SNR ; %PSK频带信号对应的噪声功率
theo_err_prb=[];
iter = 10; %迭代次数
error_rate = zeros(iter,length(SNR_dB));

%% 随机比特序列、以及PAM生成 %%
a = 2*(round(rand(1,N))-1/2);
gt = (ones(1,T/dt))'*a;  %生成T/dt列的二进制比特序列
g_PAM = gt(:)';
t1 = length(t);

%% 载波生成以及BPSK频带调制 %%
f = A*cos(2*pi*fc*t); %生成载波信号
psk = g_PAM.*f ; %PSK频带调制信号的生成

%% BPSK频带已调信号功率谱 %%
% window=boxcar(length(psk));%使用矩形窗
% nw=4;                                   %为时间带宽积，缺省值为4
% [Pxx_PAM,ff_PAM]=pmtm(g_PAM,nw,nfft,fs); %用多窗口法(NW=4)估计功率谱
% [Pxx_PSK,ff_PSK]=pmtm(psk,nw,nfft,fs); %用多窗口法(NW=4)估计功率谱
% P_PAM= [fliplr(Pxx_PAM) Pxx_PAM];
% ff_PAM = [-fliplr(ff_PAM) ff_PAM];
% P_PSK= [fliplr(Pxx_PSK) Pxx_PSK];

%% 时域信号绘图 %%
figure(1);
subplot(3,1,1);
plot(t,g_PAM,'LineWidth',2);
title('基带PAM信号(Digital Signal)');
xlabel('时间(Time(Sec))');
ylabel('幅度(Amplitude(Volts))');
grid on;
axis([0,3,-1.1,1.1]);
subplot(3,1,2);
plot(t,f);
title('载波(Carrier Signal)');
xlabel('时间(Time(Sec))');
ylabel('幅度(Amplitude(Volts))');
grid on;
axis([0,3,-1.1,1.1]);
subplot(3,1,3);
plot(t,psk,'r');
title('BPSK已调信号(BPSK Modulated Signal)');
xlabel('时间(Time(Sec))');
ylabel('幅度(Amplitude(Volts))');
grid on;
axis([0,3,-1.1,1.1]);

%% 功率谱绘图
% figure(2);
% subplot(2,1,1);
% plot(ff_PAM,P_PAM);
% title('基带PAM信号的功率谱');
% xlabel('频率(Hz)');
% ylabel('幅度');
% grid on;
% xlim([-200,200]);
% subplot(2,1,2);
% plot(ff_PSK,P_PSK);
% title('频带PSK信号的功率谱');
% xlabel('频率(Hz)');
% ylabel('幅度');
% grid on;
% xlim([-200,200]);

%% PSK信号转化为信号空间中 %%
PSK = psk ./ sqrt(1/Eb) ./ cos(2*pi*fc*t);
s = length(ta);
for n = s:s:length(PSK)
    psk_signal(1,n/s) = PSK(1,n);
end

%% AWGN信道传输 最佳接受检测  %%
for i = 1:iter
    for p = 1:length(SNR_dB)
        noise = sqrt(N0(p)) * randn(1,length(psk_signal)); %产生PSK频带信号对应的加性高斯白噪声
        out_psk = psk_signal +noise;
        for k =1:length(out_psk)
            if(out_psk(k) > 0)             %最小距离判决
                detect_bits(1,k) =  1./sqrt(1/Eb);
            else
                detect_bits(1,k) = -1./sqrt(1/Eb);
            end
        end
        theo_err_prb(p)=qfunc(sqrt(SNR(p)));
        %% erro pro
        simulated_err_prb(p) = 1 - sum(detect_bits == a./sqrt(1/Eb))/N;
        if p == 8
            r_bits_5dB = out_psk;
        elseif p == 13
            r_bits_10dB = out_psk;
        end
    end
    error_rate(i,:) = simulated_err_prb;
end
err_rate = mean(error_rate,1);
figure;
semilogy(SNR_dB,theo_err_prb,'M-X',SNR_dB,err_rate,'k-s');
grid on;
xlabel('Es/No(dB)');
ylabel('误码率');
title('BPSK调制误码率');
legend('理论误码率','仿真误码率');
scatterplot(r_bits_5dB);
title('Receiving signal constellation diagram under 5dB');
line([0 0],[-3 3],'Color','r');
grid on;
scatterplot(r_bits_10dB);
title('Receiving signal constellation diagram under 10dB');
line([0 0],[-3 3],'Color','r');
grid on;
scatterplot(detect_bits);
title('Decoded signal constellation');
line([0 0],[-3 3],'Color','r');
grid on;