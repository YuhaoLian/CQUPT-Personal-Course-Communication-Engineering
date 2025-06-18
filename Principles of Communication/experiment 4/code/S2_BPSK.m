clc;
clear;
close all;
%% 仿真参数
fc = 1000;
Rb = 100;
fs = 8000;
ts = 1/fs;

%% 升余弦滚降滤波器
Fd=1;  %输入信号的采样率
Fs=8;   %输出信号的采样率
Delay=3;  %滤波器的群时延
R=0.8;    %滚降因子
[yf,tf] = rcosine(Fd,Fs,'sqrt',R,Delay);


%% 基带信号生成
info = randi([0,1],1,10000);  % 产生待发送的二进制序列
info_map = info*2-1;   % 转化成双极码
sm_rect = rectpulse(info_map,80);           % 矩形脉冲

%% 升余弦滚降滤波的产生
sm_no = filter(yf,tf,sm_rect)/(Fs^0.5);


%% BPSK 信号生成
t1 = (0:length(sm_rect)-1)*ts;
sc1 = cos(2*pi*fc*t1);                % 高频载波
s_bpsk_rect =  sc1.*  sm_rect;                                 % BPSK信号产生
s_bpsk_no = sc1 .* sm_no(1,1:800000);


%% PAM时域绘图
 % Time-domain Waveforms
figure();
subplot(3,1,1);
plot(t1,sc1);
title('高频载波 ');
xlabel("时间 (s)");
axis([0.02 0.08 -1.5 1.5]);

subplot(3,1,2);
plot(t1,sm_rect);
title('矩形脉冲基带信号');
xlabel('时间 (s)');
axis([0.02 0.08 -1.5 1.5]);

subplot(3,1,3);
plot(t1,s_bpsk_rect);
title('BPSK调制信号');
xlabel('时间(s)');
axis([0.02 0.08 -1.5 1.5]);


%% 升余弦滚降绘图
 %Time-domain Waveforms
figure();
subplot(3,1,1);
plot(t1,sc1);
title('高频载波 ');
xlabel("时间 (s)");
axis([0.02 0.08 -1.5 1.5]);

subplot(3,1,2);
plot(t1,sm_no(1,1:800000));
title('升余弦滚降基带信号');
xlabel('时间 (s)');
axis([0.02 0.08 -1.5 1.5]);

subplot(3,1,3);
plot(t1,s_bpsk_no);
title('升余弦滚降BPSK调制信号');
xlabel('时间(s)');
axis([0.02 0.08 -1.5 1.5]);


%% 带通滤波
wc_bpf = [2*pi*(fc-Rb)/fs, 2*pi*(fc+Rb)/fs];
bpf_hn = fir1(256, wc_bpf/pi);                % 带通滤波器设计

figure();
[h,w]=freqz(bpf_hn, 1, 4096);	
subplot(3,1,1);
plot(w*fs/(2*pi),20*log10(abs(h)),'r'); 
title('BPF频率响应');
xlabel('频率 (Hz)');
axis([400 1600 -120 10]);
grid on;

s_bpsk_bpf = filter(bpf_hn, 1, s_bpsk_rect);  % 滤波
subplot(3,1,2);
plot(t1,s_bpsk_bpf);
title('经过带通滤波器的BPSK信号');
xlabel('Time (s)');
axis([0.02 0.08 -1.5 1.5]);

no_bpsk_bpf = filter(bpf_hn, 1, s_bpsk_no);  % 滤波
subplot(3,1,3);
plot(t1,no_bpsk_bpf);
title('经过带通滤波器的滚降滤波BPSK信号');
xlabel('Time (s)');
axis([0.02 0.08 -1.5 1.5]);

%% 基带与频带频率特性比较
fft_n = 2^nextpow2(length(s_bpsk_bpf));
SPSK_f = fftshift(fft(s_bpsk_rect,fft_n));
PAM_f = fftshift(fft(sm_rect,fft_n));


figure();
fHz=(-fft_n/2: fft_n/2-1)/fft_n*fs;
subplot(2,1,1);
plot(fHz,abs(PAM_f ));
grid on;
title('PAM基带信号频率特性');
xlabel('频率 (Hz)');
subplot(2,1,2);
plot(fHz,abs(SPSK_f));
grid on;
title('BPSK信号频率特性');
xlabel('频率 (Hz)');

%% 基带与频带频率特性比较
fft_n = 2^nextpow2(length(s_bpsk_bpf));
PAMNO_f = fftshift(fft(sm_no,fft_n));
SPSKNO_f = fftshift(fft(s_bpsk_no,fft_n));
PAM_f = fftshift(fft(sm_rect,fft_n));


figure();
fHz=(-fft_n/2: fft_n/2-1)/fft_n*fs;
subplot(3,1,1);
plot(fHz,abs(PAM_f ));
grid on;
title('PAM基带信号频率特性');
xlabel('频率 (Hz)');
subplot(3,1,2);
plot(fHz,abs(PAMNO_f ));
grid on;
title('经过升余弦滚降滤波的PAM基带信号频率特性');
xlabel('频率 (Hz)');
subplot(3,1,3);
plot(fHz,abs(SPSKNO_f));
grid on;
title('经过升余弦滚降滤波的BPSK信号频率特性');
xlabel('频率 (Hz)');


%% 频率特性
fft_n = 2^nextpow2(length(s_bpsk_bpf));
SPSK_f = fftshift(fft(s_bpsk_rect,fft_n));

figure();
fHz=(-fft_n/2: fft_n/2-1)/fft_n*fs;
subplot(3,1,1);
plot(fHz,abs(SPSK_f));
grid on;
title('BPSK信号频率特性');
xlabel('频率 (Hz)');

SPSK_bpf_f = fftshift(fft(s_bpsk_bpf,fft_n));
fHz=(-fft_n/2: fft_n/2-1)/fft_n*fs;
subplot(3,1,2);
plot(fHz,abs(SPSK_bpf_f));
grid on;
title('经过带通滤波器的BPSK信号频率特性');
xlabel('频率 (Hz)');

SPSK_bpf_no = fftshift(fft(no_bpsk_bpf,fft_n));
fHz=(-fft_n/2: fft_n/2-1)/fft_n*fs;
subplot(3,1,3);
plot(fHz,abs(SPSK_bpf_no));
grid on;
title('经过带通滤波器的滚降滤波BPSK信号频率特性');
xlabel('频率 (Hz)');

%% 相干解调原理
s_mul =  s_bpsk_rect.*2.*sc1;
s_mul_no = no_bpsk_bpf .* 2 .* sc1;

figure()
wc = 2*2*pi*Rb/fs;
lbf_hn = fir1(32,wc/pi);               % 低通滤波器
[h,w]=freqz(lbf_hn, 1, 4096);												
plot(w*fs/(2*pi),20*log10(abs(h)),'r'); 
title('低通滤波器频率响应');
xlabel('频率(Hz)');
axis([0 1500 -120 10]);
grid on;

s_lpf = filter(lbf_hn, 1, s_mul);      % 低通滤波
s_lpf_no = filter(lbf_hn, 1, s_mul_no);      % 低通滤波

%% 时域绘图
figure();
subplot(2,1,1);
plot(t1(1:1200),sm_rect(1:1200),'-r','LineWidth',2);
axis([0 0.15 -1.1 1.1]);
title('矩形脉冲基带信号');
xlabel('Time (s)');
grid on;

subplot(2,1,2);
plot(t1(1:1200),s_lpf(1:1200),'-b');
axis([0 0.15 -1.1 1.1]);
title('PAM矩形脉冲解调的基带信号');
xlabel('Time (s)');
grid on;
%% 频域绘图
figure();
R_f = fftshift(fft(s_lpf,fft_n));

fHz=(-fft_n/2: fft_n/2-1)/fft_n*fs;
subplot(2,1,1);
plot(fHz,abs(PAM_f ));
grid on;
title('PAM基带信号频率特性');
xlabel('频率 (Hz)');
subplot(2,1,2);
plot(fHz,abs(R_f));
grid on;
title('BPSK相干解调信号频率特性');
xlabel('频率 (Hz)');

%% 时域绘图
figure();
subplot(3,1,1);
plot(t1(1:1200),sm_rect(1:1200),'-r','LineWidth',2);
axis([0 0.15 -1.1 1.1]);
title('矩形脉冲基带信号');
xlabel('Time (s)');
grid on;

subplot(3,1,2);
plot(t1(1:1200),s_lpf(1:1200),'-b');
axis([0 0.15 -1.1 1.1]);
title('PAM矩形脉冲解调的基带信号');
xlabel('Time (s)');
grid on;

subplot(3,1,3);
plot(t1(1:1200),s_lpf_no(1:1200),'-b');
axis([0 0.15 -1.3 1.3]);
title('滚降滤波解调的基带信号');
xlabel('Time (s)');
grid on;

%% 频域绘图
figure();
R_f = fftshift(fft(s_lpf,fft_n));
RNO_f = fftshift(fft(s_lpf_no,fft_n));

fHz=(-fft_n/2: fft_n/2-1)/fft_n*fs;
subplot(3,1,1);
plot(fHz,abs(PAM_f ));
grid on;
title('PAM基带信号频率特性');
xlabel('频率 (Hz)');
subplot(3,1,2);
plot(fHz,abs(R_f));
grid on;
title('BPSK相干解调信号频率特性');
xlabel('频率 (Hz)');
subplot(3,1,3);
plot(fHz,abs(RNO_f));
grid on;
title('经过升余弦滚降滤波的BPSK相干解调信号频率特性');
xlabel('频率 (Hz)');
