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
dt=0.005;   %设定步长
t=0:dt:3;
NFFT = 2^16; %傅里叶变换点数
f = linspace(-0.5/dt, 0.5/dt, NFFT); %生成-fs/2到 fs/2的NFFT个点
int_mt = zeros(1,length(t)) ;
kf = 10:2:20; %调频指数
% kd = 1000; %调频灵敏度
recieve_signal_0dB = zeros(length(kf),length(t)-1);
recieve_signal_30dB = zeros(length(kf),length(t)-1);
Pn = zeros(1,length(kf)); %噪声功率
Ps = zeros(1,length(kf)); %信号功率
rec_SNR_dB = zeros(1,length(kf)); 
rec_SNR = zeros(1,length(kf)); 

for ii = 1:length(kf)
    %% 生成调制信号、载波信号以及已调信号
    mt = 5 * cos(2*pi*fm*t); %m(t)调制信号
    ct = cos(2*pi*fc*t); %载波
    int_mt(1) = 0;
    %对mt进行积分
    for i = 1:length(t)-1
        int_mt(i+1) = int_mt(i) + mt(i)*dt;
    end
    sft = cos(2*pi*fc*t + kf(ii).*int_mt);
    
    %% 生成无噪声情况下解调信号
    sdt = [diff(sft)];        % 先对已调信号做微分
    rt = abs(hilbert(sdt));     % 再构造解析信号求包络,得到非相干解调信号;
    rt = rt - mean(rt);
    recieve_signal_0dB(ii,:)  = rt;     

    
    %% 绘制30dB下已调信号时域图和信噪比条件下含高斯白噪声的已调信号时域图
    A = 1;
    SNRdB = 30;
    SNR = 10.^(SNRdB/10);
    sigma = A*A./SNR ;
    out_sft = zeros(length(SNRdB),length(sft));
    noise = sqrt(sigma)*randn([1 length(t)]);
    out_sft = noise + sft;
    sdt = [diff(out_sft)];        % 先对已调信号做微分
    rt = abs(hilbert(sdt));     % 再构造解析信号求包络,得到非相干解调信号;
    rt = rt - mean(rt);
    recieve_signal_30dB(ii,:) = rt;
    Ps(ii) = trapz(0:dt:3-dt,recieve_signal_0dB(ii,:).^2);
    Pn(ii) = trapz(0:dt:3-dt,(recieve_signal_30dB(ii,:)-recieve_signal_0dB(ii,:)).^2);

end

rec_SNR = Ps./Pn;
rec_SNR_dB = 10*log10(rec_SNR );
kf_str = ['10' , '12'  ,'14'  ,'16' , '18' , '20'];
figure
for i = 1:length(kf)
    subplot(length(kf),1,i)
    plot(0:dt:3-dt,recieve_signal_30dB(i,:));
    title(['解调信号时域图——调频灵敏度为为',kf_str(1,2*i-1:2*i)]);
end

