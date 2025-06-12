%姓名：连宇昊；学号：2020210140；题号：31；
clc;
clear;


fs = 8000;  %模拟抽样频率
t = 0:1/fs:3;
xt = cos(4*pi*t) + sin(10*pi*t);

fs = 20;   %因为fm=5Hz,根据奈奎斯特采样定理，fs>=2fm=10Hz，此处选择20Hz

% 此时fs =20hz ,x(n)=cos(n*pi/5) + sin(n*pi/2)
% 所以此时序列周期为20，即序列取40个采样点即可
n = 0:39; %取40个采样点即可，此处多取一点便于观察
tt = n / fs;   % tt = nTs带入
xn = cos(4*pi*tt) + sin(10*pi*tt);

% 因为序列为周期序列，对其谱分析时，N = kT 即整数倍的序列周期，此处N取40
N = 40;
Xn = fft(xn,N);
Xm = abs(Xn);
Xp = angle(Xn);


figure(1)
subplot(2,2,1)
plot(t,xt);
xlabel('t(s)');
ylabel('x(t)');
title('x(t)在0~3秒内时域连续波形图');
subplot(2,2,2)
stem(n,xn,'fill');
xlabel('n');
ylabel('x(n)');
title('x(n)序列2个周期的时域波形图');
subplot(2,2,3)
stem(0:N-1,Xm,'fill');
xlabel('k');
ylabel('X(k)');
title('X(k)离散幅度谱（40点FFT）');
subplot(2,2,4)
stem(0:N-1,Xp,'fill');
xlabel('k');
ylabel('phase');
title('X(k)离散相位谱（40点FFT）');
