%% <<<<<<<<<<<<<<<<<<< BPSK Modulation and Demodulation >>>>>>>>>>>>>>>>>>>
clc;
clear all;
close all;

%% ******************* 初始化 ******************** %%
M = 4; %四进制
fc = 8; %载波频率
Rb = fc; %码元速率
T = 1 / Rb; % 一个码元传送的时间
dt = 1/2000; %一个码元的抽样时间
fs = 1/dt; %采样频率
nfft=1024; %采样的点数
Nbit = 100000;%10^5个二进制比特
N = Nbit/log2(M);
t=0:dt:N*T-dt;% 所有码元的时间
ta=0:dt:T-dt; % 一个码元时间向量
%此时，一个码元序列被抽样T/dt次

A = 1;%振幅等于1
SNR_dB = -2:10; %信噪比——SNR dB值
SNR = 10.^(SNR_dB./10); %信噪比——SNR
Es= A^2*T/2; % QPSK频带信号对应的维度的每符号能量
N0 = Es./SNR ; %QPSK频带信号对应的噪声功率
theo_err_prb=[];
iter = 10; %迭代次数
error_rate = zeros(iter,length(SNR_dB));

%% 随机比特序列、以及PAM生成 %%
a = 2*(round(rand(1,Nbit))-1/2); %二进制比特随机序列
for i = 1:N
    bI(1,i) = a(1,2*i);
    bq(1,i) = a(1,2*i - 1);
end

gtl = (ones(1,T/dt))'*bI;  %生成T/dt列的二进制比特序列
g_PAM_I = gtl(:)';
gtq = (ones(1,T/dt))'*bq;  %生成T/dt列的二进制比特序列
g_PAM_q = gtq(:)';
t1 = length(t);

%% 载波生成以及QPSK频带调制 %%
f_I = A*cos(2*pi*fc*t); %生成载波信号
f_q = A*sin(2*pi*fc*t); %生成载波信号
qpsk = g_PAM_I.*f_I  - g_PAM_q.*f_q; %QPSK频带调制信号的生成

%% QPSK频带已调信号功率谱 %%
% window=boxcar(length(psk));%使用矩形窗
% nw=4;                                   %为时间带宽积，缺省值为4
% [Pxx_PAM,ff_PAM]=pmtm(g_PAM,nw,nfft,fs); %用多窗口法(NW=4)估计功率谱
% [Pxx_PSK,ff_PSK]=pmtm(psk,nw,nfft,fs); %用多窗口法(NW=4)估计功率谱
% P_PAM= [fliplr(Pxx_PAM) Pxx_PAM];
% ff_PAM = [-fliplr(ff_PAM) ff_PAM];
% P_PSK= [fliplr(Pxx_PSK) Pxx_PSK];
% ff_PSK = [-fliplr(ff_PSK) ff_PSK];

%% 时域信号绘图 %%
figure(1);
subplot(3,1,1);
plot(t,g_PAM_I,'LineWidth',2);
title('I(t)基带PAM信号(Digital Signal)');
xlabel('时间(Time(Sec))');
ylabel('幅度(Amplitude(Volts))');
grid on;
axis([0,3,-1.1,1.1]);
subplot(3,1,2);
plot(t,g_PAM_q,'LineWidth',2);
title('Q(t)基带PAM信号(Digital Signal)');
xlabel('时间(Time(Sec))');
ylabel('幅度(Amplitude(Volts))');
grid on;
axis([0,3,-1.1,1.1]);
subplot(3,1,3);
plot(t,qpsk);
title('QPSK已调信号(QPSK Modulated Signal)');
xlabel('时间(Time(Sec))');
ylabel('幅度(Amplitude(Volts))');
grid on;
axis([0,3,-1.5,1.5]);

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

%% QPSK信号转化为信号空间中 %%
QPSK_I = g_PAM_I./ sqrt(1/Es);
QPSK_q = g_PAM_q./ sqrt(1/Es);
QPSK_signal = QPSK_I + 1j * QPSK_q;
s = length(ta);
for n = s:s:length(QPSK_signal)
    qpsk_signal(1,n/s) = QPSK_signal(1,n);
end

%% AWGN信道传输 %%
s00 = 1/sqrt(1/Es)*(-1-1j);
s01 = 1/sqrt(1/Es)*(-1+1j);
s10 = 1/sqrt(1/Es)*(1-1j);
s11 = 1/sqrt(1/Es)*(1+1j);

for i = 1:iter
    for p = 1:length(SNR_dB)
        noise1 = sqrt(N0(p)) * randn(1,length(qpsk_signal)); %产生QPSK频带信号对应的加性高斯白噪声
        noise2 = sqrt(N0(p)) * randn(1,length(qpsk_signal)); %产生QPSK频带信号对应的加性高斯白噪声
        out_qpsk = qpsk_signal +(noise1 + 1j*noise2);
        for k =1:length(out_qpsk)
            d00 = abs(out_qpsk(k)-s00);
            d01 = abs(out_qpsk(k)-s01);
            d10 = abs(out_qpsk(k)-s10);
            d11 = abs(out_qpsk(k)-s11);
            if((d00<d01)&&(d00<d10)&&(d00<d11))   %最小距离判决
                detect(1,k) = -1/sqrt(1/Es)-1j/sqrt(1/Es);
            elseif((d01<d00)&&(d01<d10)&&(d01<d11))
                detect(1,k) = -1/sqrt(1/Es)+1j/sqrt(1/Es);
            elseif((d10<d00)&&(d10<d01)&&(d10<d11))
                detect(1,k) = 1/sqrt(1/Es)-1j/sqrt(1/Es);
            else
                detect(1,k) = 1/sqrt(1/Es)+1j/sqrt(1/Es);
            end
        end
        theo_err_prb(p)=2*qfunc(sqrt(SNR(p)))-qfunc(sqrt(SNR(p))).*qfunc(sqrt(SNR(p)));
        %% erro pro
        simulated_err_prb(p) = 1 - sum(detect == qpsk_signal)/N;
        if p == 8
            r_bits_5dB = out_qpsk;
        elseif p == 13
            r_bits_10dB = out_qpsk;
        end
    end
    error_rate(i,:) = simulated_err_prb;
end
err_rate = mean(error_rate,1);
figure;
semilogy(SNR_dB,theo_err_prb,'M-X',SNR_dB,err_rate,'k-s');
grid on;
xlabel('Es/No(dB)');
ylabel('BER');
title('QPSK调制BER(迭代10次)');
legend('理论BER','仿真BER');
scatterplot(r_bits_5dB);
title('Receiving signal constellation diagram under 5dB');
line([0 0],[-3 3],'Color','r');
line([-3 3],[0 0],'Color','r');
grid on;
scatterplot(r_bits_10dB);
title('Receiving signal constellation diagram under 10dB');
line([0 0],[-3 3],'Color','r');
line([-3 3],[0 0],'Color','r');
grid on;
scatterplot(detect);
title('Decoded signal constellation');
line([0 0],[-3 3],'Color','r');
line([-3 3],[0 0],'Color','r');
grid on;