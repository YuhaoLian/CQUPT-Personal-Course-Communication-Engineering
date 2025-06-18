clear all;
close all;
clc;

%% DSB的调制及非同相位解调

%% 初始化
% fs = 128;  %抽样频率有点小
fs = 99999;
fm = 1000;
fc = 10000;
dt = 1/fs;
t = 0: dt : 50;
NFFT = 2^16; %傅里叶变换点数
f = linspace(-fs/2, fs/2, NFFT); %生成-fs/2到 fs/2的NFFT个点
iter = 4;
phase = [pi/8 pi/4 pi/3 pi/2];
phase_str = ['pi/8', 'pi/4' ,'pi/3' ,'pi/2'];

% 生成调制信号、载波信号以及已调信号
s_m = cos(2*pi*fm*t); %m(t)调制信号
s_c = cos(2*pi*fc*t); %载波
s_DSB = s_m .*s_c; %DSB调制信号

%% 低通滤波
wc=3*2*pi*fm/fs;                        %解调器低通滤波器截止频率，3fm
lbf_hn=fir1(32,wc/pi);                     %低通滤波器系统函数，FIR滤波器
freqz(lbf_hn, 1, 1024, fs);

figure;
%% 生成调制信号、载波信号以及已调信号
for i =1:iter
    % DSB
    s_cc = cos(2*pi*fc*t + phase(1,i));
    s_cor_temp_DSB = s_DSB.* 2 .* s_cc; %通过乘法器
    dem_cor_out_DSB = filter(lbf_hn,1,s_cor_temp_DSB); %通过低通滤波器
    subplot(2,2,i)
    plot(t,dem_cor_out_DSB); grid on;
    xlabel('t(s)')
    title(['解调DSB信号相位相差',phase_str(4*i - 3:i*4)]);
    % xlim([0, 5]); %这个不平滑
    xlim([0, 0.01]); %这个没毛病
    ylim([-1, 1]); %这个没毛病
end
    