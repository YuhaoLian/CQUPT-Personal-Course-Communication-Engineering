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
fm = 5;
fc = 25;
dt=0.0001;   %设定步长
t=0:dt:3;
NFFT = 2^16; %傅里叶变换点数
f = linspace(-0.5/dt, 0.5/dt, NFFT); %生成-fs/2到 fs/2的NFFT个点
int_mt = zeros(1,length(t)) ;
kf = 10:10:50; %调频指数
sft = zeros(length(kf),length(t));
kf_str = ['10','20','30','40','50'];

%% 生成调制信号、载波信号以及已调信号
mt = 5 * cos(2*pi*fm*t); %m(t)调制信号
ct = cos(2*pi*fc*t); %载波
int_mt(1) = 0;
%对mt进行积分
for i = 1:length(t)-1
    int_mt(i+1) = int_mt(i) + mt(i)*dt;
end

for k = 1:length(kf)
    sft(k,:) = cos(2*pi*fc*t + kf(k).*int_mt);
    [f2,Fm]=T2F(t,sft(k,:));%傅里叶变换得到频谱
    subplot(length(kf),1,k)
    plot(f2,abs(Fm));
    axis([-100 100 0 max(abs(Fm))]);
    title(['调频指数kf =',kf_str(1,2*k-1:2*k),' 下已调信号频谱']);
    xlabel('f(Hz)')
end



