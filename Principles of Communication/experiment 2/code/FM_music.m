clc;
clear;
close all;

% 确定仿真参数
% 调制信号幅度：5
% 调制信号频率：5Hz
% 载频：25Hz
% 调频指数：10
% 小信噪比：10dB
% 大信噪比：40dB


%% 初始化
% fs = 99999;
[xn, fs] = audioread('motherland.wav');
fm = 5;
fc = 25;
dt=1/fs;   %设定步长
t=0:dt:10.4 - dt;
NFFT = 2^16; %傅里叶变换点数
f = linspace(-0.5/dt, 0.5/dt, NFFT); %生成-fs/2到 fs/2的NFFT个点
int_mt = zeros(1,length(t)) ;
kf = 10; %调频指数


%% 需要同学生成调制信号、载波信号以及已调信号
mt = xn'; %m(t)调制信号
ct = cos(2*pi*fc*t); %载波
int_mt(1) = 0;
%对mt进行积分
for i = 1:length(t)-1
    int_mt(i+1) = int_mt(i) + mt(i)*dt;
end
sft = cos(2*pi*fc*t + kf.*int_mt);

%% 绘制三类信号的时域图
figure(1) %绘制图形
subplot(311);
plot(t, mt);
axis([0 10.4 -0.5 0.5])
grid on;
xlabel('t');
title('调制信号时域图');
subplot(312);
plot(t, ct);
axis([0 1 -1.1 1.1])
xlabel('t');
grid on;
title('载波时域图');
subplot(313);
plot(t,sft) ;
axis([0 1 -1.1 1.1])
xlabel('t')
title('已调信号时域图');
grid on;

%% 绘制调制信号与已调信号的频域图
[f1,fmm]=T2F(t,mt);
[f2,Fm]=T2F(t,sft);%傅里叶变换得到频谱
figure(2) %绘制已调信号时频域图
subplot(211);
plot(f1,abs(fmm)) ;
axis([-fs/2 fs/2 0 max(abs(fmm))]);
xlabel('f(Hz)')
title('调制信号频谱');
subplot(212)
plot(f2,abs(Fm));
axis([-fs/2 fs/2 0 max(abs(Fm))]);
title('已调信号频谱');
xlabel('f(Hz)')


%**************figure(3)无噪声解调信号的时域图******************
figure(3)
subplot(2,1,1);plot(t,mt);		           %绘制调制信号的时域图
xlabel('时间t');
title('调制信号的时域图');
subplot(2,1,2);plot(t,sft);		       %绘制已调信号的时域图
xlabel('时间t');
title('无噪声条件下已调信号的时域图');

%% 绘制已调信号时域图和信噪比条件下含高斯白噪声的已调信号时域图
A = 1;
SNRdB = 10:5:40;
SNRdB_str = ['10' ,'15' ,'20' ,'25' ,'30' ,'35' ,'40'];
SNR = 10.^(SNRdB/10);
sigma = A*A./SNR ;
out_sft = zeros(length(SNRdB),length(sft));
for i = 1:length(SNRdB)
    noise = sqrt(sigma(i))*randn([1 length(t)]);
    out_sft(i,:) = noise + sft;
    subplot(length(SNRdB),1,i)
    plot(t,out_sft(i,:)) ;
    axis([0 1 -1.1 1.1])
    xlabel('t')
    title(['已调信号时域图信噪比为',SNRdB_str(1,2*i-1:2*i),'dB']);
end

sound(sft,fs);
