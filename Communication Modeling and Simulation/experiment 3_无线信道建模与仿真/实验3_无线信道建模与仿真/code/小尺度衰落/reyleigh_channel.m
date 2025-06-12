clear;
close all;
%T.S.Rappaport著 《无线通信原理与应用（第二版）》 电子工业出版社 P53 图5.23
% 参数初始化
N = 1024*1024; % 频域样点数
v = 33.33; % 移动终端的速度（米/秒）
c = 3e8; % 光速（米/秒）
fc = 900e6; % 载波频率
lamda = c / fc; % 载波波长
fm = v / lamda; % 最大多普勒频移
deltaf = 2 * fm / (N-1); % 频域采样间隔
T = 1 / deltaf; % 时域采样间隔
g = linspace(-fm, fm, N); % 频率点的值

% 产生有多普勒频移引起的多普勒功率普密度函数 sqrt(SEZ(f))
S = zeros(N, 1);
for i = 2:N-1
f = fc + g(i);
S(i) =1.5 / (pi * fm * (1 - ((f - fc) / fm) ^2) ^0.5);
end
figure(1)
plot(S)
S = sqrt(S);

% 产生同相噪声源
Ia = randn(N/2, 1);
Ib = randn(N/2, 1);
%Ig_positive和Ig_negative是以0点为中心对称的
Ig_positive = Ia + Ib*j; % 得到Ig的正数部分
Ig_negative = flipud(conj(Ig_positive)); % 得到Ig的负数部分
Ig = [sqrt(real(Ig_positive).^2 + imag(Ig_positive).^2); sqrt(real(Ig_negative).^2 + imag(Ig_negative).^2)];
Ifad = ifft(Ig .* S); %反傅立叶变换
Ifad = abs(Ifad).^2; % 求模的平方

% 产生正交噪声源
Qa = randn(N/2, 1);
Qb = randn(N/2, 1);
%Qg_positive和Qg_negative是以0点为中心对称的
Qg_positive = Qa + Qb*j; % 得到Qg的正数部分
Qg_negative = flipud(conj(Qg_positive));% 得到Qg的负数部分
Qg = [sqrt(real(Qg_positive).^2 + imag(Qg_positive).^2); sqrt(real(Qg_negative).^2 + imag(Qg_negative).^2)];
Qfad = ifft(Qg .* S);%反傅立叶变换
Qfad = abs(Qfad) .* exp(j*(angle(Qfad)-pi/2)); % 相位旋转90度
Qfad = abs(Qfad).^2;% 求模的平方

% 合并两路噪声源
rt = sqrt(Ifad + Qfad);
mean = cumsum(rt);
mean = mean(N)/N;
rt=rt/mean;
rt=rt';

figure(2)
hist(rt(1,10000:1:N/2),100);
grid on;
title('包络幅度统计量');

figure(3)
[f,xi]=ksdensity(rt(1,10000:1:N/2));
plot(xi,f)
grid on;
title('包络的概率密度函数');

figure(4)
plot(0:1:256,10*log10(rt(1,10000:1:(10000+256))))
axis tight

grid on;
xlabel('采样点');
ylabel('包络增益');
title('某一段时间内的包络曲线 (256个采样点)')

