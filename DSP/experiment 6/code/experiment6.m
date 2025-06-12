clc;
clear;
close all;

%% 用脉冲响应不变法设计 IIR 数字低通滤波器
omega_p = 0.25 * pi;    %滤波器的通带截止频率
omega_s = 0.4 * pi; %滤波器的通带截止频率
Rp = 1;             %滤波器的通带衰减指标
As = 15;            %滤波器的阻带衰减指标
ripple = 10^(-Rp/20);   %滤波器的通带衰减对应的幅度值
Attn = 10^(-As/20);     %滤波器的阻带衰减对应的幅度值

Fs = 2000;
T = 1 / Fs;
Omgp = omega_p * Fs;
Omgs = omega_s * Fs;
[n,Omgc] = buttord(Omgp,Omgs,Rp,As,'s');    %计算阶数n和截止频率
[z0,p0,k0] = buttap(n);         %设计归一化的巴特沃思模拟滤波器原型
ba1 = k0 * real(poly(z0));      %求原型滤波器的系数b
aa1 = real(poly(p0));       %求原型滤波器的系数a
[ba,aa] = lp2lp(ba1,aa1,Omgc);      %将模拟低通原型滤波器转换为实际的低通滤波器：ba为分子多项式系数、aa为分母多项式系数
[bd,ad] = impinvar(ba,aa,Fs);       %利用脉冲响应不变法将模拟滤波器转换为数字滤波器
[H,w] = freqz(bd,ad);       %求数字系统的频率特性 
dbH = 20 * log10((abs(H) + eps)/max(abs(H)));
figure(1);
subplot(2,2,1);
plot(w / pi,abs(H));
ylabel('|H|');
title("幅频响应");
axis([0,1,0,1.1]);
set(gca,'XTickMode','manual','XTick',[0,0.25,0.4,1]);
set(gca,'YTickMode','manual','YTick',[0,Attn,ripple,1]);
grid;
subplot(2,2,2);
plot(w / pi,angle(H)/pi);
ylabel('\phi');
title("相频响应");
axis([0,1,-1,1]);
set(gca,'XTickMode','manual','XTick',[0,0.25,0.4,1]);
set(gca,'YTickMode','manual','YTick',[-1,0,1]);
grid;
subplot(2,2,3);
plot(w / pi,dbH);
ylabel('dB');
xlabel('频率（\pi）');
title("幅频响应(dB)");
axis([0,1,-40,5]);
set(gca,'XTickMode','manual','XTick',[0,0.25,0.4,1]);
set(gca,'YTickMode','manual','YTick',[-50,-15,-1,0]);
grid;
subplot(2,2,4);
zplane(bd,ad);
axis([-1.1,1.1,-1.1,1.1]); 
title('零极点图');

%% 用双线性变换法设计 IIR 数字低通滤波器
omega_p = 0.25 * pi;    %滤波器的通带截止频率
omega_s = 0.4 * pi; %滤波器的通带截止频率
Rp = 1;             %滤波器的通带衰减指标
As = 15;            %滤波器的阻带衰减指标
ripple = 10^(-Rp/20);   %滤波器的通带衰减对应的幅度值
Attn = 10^(-As/20);     %滤波器的阻带衰减对应的幅度值

Fs = 100;
T = 1 / Fs;
Omgp = tan(omega_p / 2) * Fs * 2;
Omgs = tan(omega_s / 2) * Fs * 2;
[n,Omgc] = buttord(Omgp,Omgs,Rp,As,'s');    %计算阶数n和截止频率
[z0,p0,k0] = buttap(n);         %设计归一化的巴特沃思模拟滤波器原型
ba1 = k0 * real(poly(z0));      %求原型滤波器的系数b
aa1 = real(poly(p0));       %求原型滤波器的系数a
[ba,aa] = lp2lp(ba1,aa1,Omgc);      %将模拟低通原型滤波器转换为实际的低通滤波器：ba为分子多项式系数、aa为分母多项式系数
[bd,ad] = bilinear(ba,aa,Fs);       %利用双线性变换法将模拟滤波器转换为数字滤波器
[H,w] = freqz(bd,ad);       %求数字系统的频率特性
dbH = 20 * log10((abs(H) + eps)/max(abs(H)));
figure(2);
subplot(2,2,1);
plot(w / pi,abs(H));
ylabel('|H|');
title("幅频响应");
axis([0,1,0,1.1]);
set(gca,'XTickMode','manual','XTick',[0,0.25,0.4,1]);
set(gca,'YTickMode','manual','YTick',[0,Attn,ripple,1]);
grid;
subplot(2,2,2);
plot(w / pi,angle(H)/pi);
ylabel('\phi');
title("相频响应");
axis([0,1,-1,1]);
set(gca,'XTickMode','manual','XTick',[0,0.25,0.4,1]);
set(gca,'YTickMode','manual','YTick',[-1,0,1]);
grid;
subplot(2,2,3);
plot(w / pi,dbH);
ylabel('dB');
xlabel('频率（\pi）');
title("幅频响应(dB)");
axis([0,1,-40,5]);
set(gca,'XTickMode','manual','XTick',[0,0.25,0.4,1]);
set(gca,'YTickMode','manual','YTick',[-50,-15,-1,0]);
grid;
subplot(2,2,4);
zplane(bd,ad);
axis([-1.1,1.1,-1.1,1.1]); 
title('零极点图');

%% 实验内容
% 1.设计滤波器
omega_p = 0.2 * pi;    %滤波器的通带截止频率
omega_s = 0.35 * pi; %滤波器的通带截止频率
Rp = 1;             %滤波器的通带衰减指标
As = 15;            %滤波器的阻带衰减指标
ripple = 10^(-Rp/20);   %滤波器的通带衰减对应的幅度值
Attn = 10^(-As/20);     %滤波器的阻带衰减对应的幅度值

Fs = 10;
T = 1 / Fs;
Omgp = tan(omega_p / 2) * Fs * 2;
Omgs = tan(omega_s / 2) * Fs * 2;
[n,Omgc] = buttord(Omgp,Omgs,Rp,As,'s');    %计算阶数n和截止频率
[z0,p0,k0] = buttap(n);         %设计归一化的巴特沃思模拟滤波器原型
ba1 = k0 * real(poly(z0));      %求原型滤波器的系数b
aa1 = real(poly(p0));       %求原型滤波器的系数a
[ba,aa] = lp2lp(ba1,aa1,Omgc);      %将模拟低通原型滤波器转换为实际的低通滤波器：ba为分子多项式系数、aa为分母多项式系数
[bd,ad] = bilinear(ba,aa,Fs);       %利用双线性变换法将模拟滤波器转换为数字滤波器
[H,w] = freqz(bd,ad);       %求数字系统的频率特性
dbH = 20 * log10((abs(H) + eps)/max(abs(H)));
figure(3);
subplot(2,2,1);
plot(w / pi,abs(H));
ylabel('|H|');
title("幅频响应");
axis([0,1,0,1.1]);
set(gca,'XTickMode','manual','XTick',[0,0.2,0.35,1]);
set(gca,'YTickMode','manual','YTick',[0,Attn,ripple,1]);
grid;
subplot(2,2,2);
plot(w / pi,angle(H)/pi);
ylabel('\phi');
title("相频响应");
axis([0,1,-1,1]);
set(gca,'XTickMode','manual','XTick',[0,0.2,0.35,1]);
set(gca,'YTickMode','manual','YTick',[-1,0,1]);
grid;
subplot(2,2,3);
plot(w / pi,dbH);
ylabel('dB');
xlabel('频率（\pi）');
title("幅频响应(dB)");
axis([0,1,-40,5]);
set(gca,'XTickMode','manual','XTick',[0,0.2,0.35,1]);
set(gca,'YTickMode','manual','YTick',[-50,-15,-1,0]);
grid;
subplot(2,2,4);
zplane(bd,ad);
axis([-1.1,1.1,-1.1,1.1]); 
title('零极点图');

% 滤波
xn = [-4,-2,0,-4,-6,-4,-2,-4,-6,-6,-4,-4,-6,-2,6,12,8,...
        0,-16,-38,-60,-84,-90,-66,-32,-4,-2,-4,8,12,12,10,6,6,6,...
        4,0,0,0,0,0,-2,-4,0,0,0,-2,-2,0,0,-2,-2,-2,-2,0];
figure(4);
subplot(2,2,1);
stem(0:length(xn)-1,xn,'.');
grid;
title("滤波前波形");
xk = fft(xn,1024);
k = 0:1023;
k = k / 512;
subplot(2,2,3);
X = abs(xk);
plot(k,X);
ylabel('|X|');
title("幅频响应");
xlabel('\omega/\pi');
grid;
yn = filter(bd,ad,xn);
subplot(2,2,2);
stem(0:length(yn)-1,yn,'.');
grid;
title("滤波后波形");
yk = fft(yn,1024);
k = 0:1023;
k = k / 512;
subplot(2,2,4);
Y = abs(yk);
plot(k,Y);
ylabel('|Y|');
title("幅频响应");
xlabel('\omega/\pi');
grid;

%% 选做题
% 读取音频
D = 2;
[xn, fs] = audioread('motherland.wav');
yn1 = xn(1:D:length(xn));
% sound(yn1,fs/D);   %播放损坏音质版本

t = (0:83199)/8000;



%设计数字滤波器
omega_p = 0.45 * pi;    %滤波器的通带截止频率
omega_s = 0.55 * pi; %滤波器的阻带截止频率
Rp = 1;             %滤波器的通带衰减指标
As = 30;            %滤波器的阻带衰减指标
ripple = 10^(-Rp/20);   %滤波器的通带衰减对应的幅度值
Attn = 10^(-As/20);     %滤波器的阻带衰减对应的幅度值

Fs = fs;
T = 1 / Fs;
Omgp = omega_p * Fs;
Omgs = omega_s * Fs;
[n,Omgc] = buttord(Omgp,Omgs,Rp,As,'s');    %计算阶数n和截止频率
[z0,p0,k0] = buttap(n);         %设计归一化的巴特沃思模拟滤波器原型
ba1 = k0 * real(poly(z0));      %求原型滤波器的系数b
aa1 = real(poly(p0));       %求原型滤波器的系数a
[ba,aa] = lp2lp(ba1,aa1,Omgc);      %将模拟低通原型滤波器转换为实际的低通滤波器：ba为分子多项式系数、aa为分母多项式系数
[bd,ad] = impinvar(ba,aa,Fs);       %利用脉冲响应不变法将模拟滤波器转换为数字滤波器
[H,w] = freqz(bd,ad);       %求数字系统的频率特性 
dbH = 20 * log10((abs(H) + eps)/max(abs(H)));
figure(5);
subplot(2,2,1);
plot(w / pi,abs(H));
ylabel('|H|');
title("幅频响应");
axis([0,1,0,1.1]);
set(gca,'XTickMode','manual','XTick',[0,0.45,0.55,1]);
set(gca,'YTickMode','manual','YTick',[0,Attn,ripple,1]);
grid;
subplot(2,2,2);
plot(w / pi,angle(H)/pi);
ylabel('\phi');
title("相频响应");
axis([0,1,-1,1]);
set(gca,'XTickMode','manual','XTick',[0,0.45,0.55,1]);
set(gca,'YTickMode','manual','YTick',[-1,0,1]);
grid;
subplot(2,2,3);
plot(w / pi,dbH);
ylabel('dB');
xlabel('频率（\pi）');
title("幅频响应(dB)");
axis([0,1,-40,5]);
set(gca,'XTickMode','manual','XTick',[0,0.45,0.55,1]);
set(gca,'YTickMode','manual','YTick',[-50,-30,-1,0]);
axis([0,1,-60,5]);
grid;
subplot(2,2,4);
zplane(bd,ad);
axis([-1.1,1.1,-1.1,1.1]); 
title('零极点图');

x = filter(bd,ad,xn);
x = x(1:D:length(x));
% yn2 = x(1:D:length(x));
sound(x,fs/D);    %播放抗混叠版
% audiowrite("鬼哭狼嚎版一条大河.wav",yn2,fs/D);

figure(10);
subplot(3,1,1);
plot(t,xn);
grid on;
title("2-4s原序列");
xlabel('t（s）');
axis([2,4,-0.6,0.6]);
subplot(3,1,2);
plot(t(1:D:length(xn)),yn1);
grid on;
xlabel('t（s）');
title("2-4sD=2抽取序列");
axis([2,4,-0.6,0.6]);
subplot(3,1,3);
plot(t(1:D:length(xn)),x);
axis([2,4,2*min(x)+min(x)*0.1,2*max(x)+max(x)*0.1]);
xlabel('t（s）');
title("2-4s滤波后序列");
grid on;

N=2048;
figure(11);
Xn=1/fs*fft(xn(8000:8199),N); % 从xn中取200点，N可取2018
subplot(3,1,1);
plot((0:N/2-1)*fs/N,abs(Xn(1:N/2))); % 模拟域幅度谱
grid on;
title("原音频幅度频谱");
xlabel('f（Hz）');
Yn1=D/fs*fft(yn1(4000:4099),N); % 2点取1点后，200点长变成了100点长
subplot(3,1,2);
plot((0:N/2-1)*fs/(N*D),abs(Yn1(1:N/2)));% 模拟域幅度谱
grid on;
title("D = 2原音频幅度频谱");
xlabel('f（Hz）');
subplot(3,1,3);
X1=D/fs*fft(x(4000:4099),N);
plot((0:N/2-1)*fs/(N*D),abs(X1(1:N/2)));% 模拟域幅度谱
grid on;
title("D = 2抗混叠音频幅度频谱");
xlabel('f（Hz）');

%% 思考题1
omega_p = 0.2 * pi;    %滤波器的通带截止频率
omega_s = 0.3 * pi; %滤波器的通带截止频率
Rp = 1;             %滤波器的通带衰减指标
As = 20;            %滤波器的阻带衰减指标
ripple = 10^(-Rp/20);   %滤波器的通带衰减对应的幅度值
Attn = 10^(-As/20);     %滤波器的阻带衰减对应的幅度值

Fs = 2000;
T = 1 / Fs;
Omgp = omega_p * Fs;
Omgs = omega_s * Fs;
[n,Omgc] = buttord(Omgp,Omgs,Rp,As,'s');    %计算阶数n和截止频率
[z0,p0,k0] = buttap(n);         %设计归一化的巴特沃思模拟滤波器原型
ba1 = k0 * real(poly(z0));      %求原型滤波器的系数b
aa1 = real(poly(p0));       %求原型滤波器的系数a
[ba,aa] = lp2lp(ba1,aa1,Omgc);      %将模拟低通原型滤波器转换为实际的低通滤波器：ba为分子多项式系数、aa为分母多项式系数
[bd,ad] = impinvar(ba,aa,Fs);       %利用脉冲响应不变法将模拟滤波器转换为数字滤波器
[H,w] = freqz(bd,ad);       %求数字系统的频率特性 
dbH = 20 * log10((abs(H) + eps)/max(abs(H)));
figure(6);
subplot(2,2,1);
plot(w / pi,abs(H));
ylabel('|H|');
title("幅频响应");
axis([0,1,0,1.1]);
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickMode','manual','YTick',[0,Attn,ripple,1]);
grid;
subplot(2,2,2);
plot(w / pi,angle(H)/pi);
ylabel('\phi');
title("相频响应");
axis([0,1,-1,1]);
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickMode','manual','YTick',[-1,0,1]);
grid;
subplot(2,2,3);
plot(w / pi,dbH);
ylabel('dB');
xlabel('频率（\pi）');
title("幅频响应(dB)");
axis([0,1,-40,5]);
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickMode','manual','YTick',[-50,-20,-1,0]);
grid;
subplot(2,2,4);
zplane(bd,ad);
axis([-1.1,1.1,-1.1,1.1]); 
title('零极点图');


%% 思考题2
omega_p = 0.6 * pi;    %滤波器的通带截止频率
omega_s = 0.4 * pi; %滤波器的通带截止频率
Rp = 2;             %滤波器的通带衰减指标
As = 30;            %滤波器的阻带衰减指标
ripple = 10^(-Rp/20);   %滤波器的通带衰减对应的幅度值
Attn = 10^(-As/20);     %滤波器的阻带衰减对应的幅度值

Fs = 100;
T = 1 / Fs;
Omgp = tan(omega_p / 2) * Fs * 2;
Omgs = tan(omega_s / 2) * Fs * 2;
[n,Omgc] = buttord(Omgp,Omgs,Rp,As,'s');    %计算阶数n和截止频率
[ba,aa] = butter(n,Omgc,'high','s');
[bd,ad] = bilinear(ba,aa,Fs);       %利用双线性变换法将模拟滤波器转换为数字滤波器
[H,w] = freqz(bd,ad);       %求数字系统的频率特性
dbH = 20 * log10((abs(H) + eps)/max(abs(H)));
figure(7);
subplot(2,2,1);
plot(w / pi,abs(H));
ylabel('|H|');
title("幅频响应");
axis([0,1,0,1.1]);
set(gca,'XTickMode','manual','XTick',[0,0.4,0.6,1]);
set(gca,'YTickMode','manual','YTick',[0,Attn,ripple,1]);
grid;
subplot(2,2,2);
plot(w / pi,angle(H)/pi);
ylabel('\phi');
title("相频响应");
axis([0,1,-1,1]);
set(gca,'XTickMode','manual','XTick',[0,0.4,0.6,1]);
set(gca,'YTickMode','manual','YTick',[-1,0,1]);
grid;
subplot(2,2,3);
plot(w / pi,dbH);
ylabel('dB');
xlabel('频率（\pi）');
title("幅频响应(dB)");
axis([0,1,-40,5]);
set(gca,'XTickMode','manual','XTick',[0,0.4,0.6,1]);
set(gca,'YTickMode','manual','YTick',[-50,-30,-2,0]);
grid;
subplot(2,2,4);
zplane(bd,ad);
axis([-1.1,1.1,-1.1,1.1]); 
title('零极点图');

%% 思考题3
wp=[0.2,0.6];
ws=[0.15,0.65];
rp=1;
as=45;
ripple = 10^(-Rp/20);   %滤波器的通带衰减对应的幅度值
Attn = 10^(-As/20);     %滤波器的阻带衰减对应的幅度值
[n1,wc1]=buttord(wp,ws,rp,as);
[b1,a1]=butter(n1,wc1);
figure(8);
[h1,w1]=freqz(b1,a1);
db1=20*log10(abs(h1)/max(abs(h1)));
subplot(221);
plot(w1/pi,abs(h1));
xlabel('\omega');
ylabel('|H|');
title("幅频响应");
axis([0,1,0,1.1]);
set(gca,'XTickMode','manual','XTick',[0,0.15,0.2,0.6,0.65,1]);
set(gca,'YTickMode','manual','YTick',[0,Attn,ripple,1]);
grid;
subplot(222);
plot(w1/pi,angle(h1)/pi);
ylabel('\phi');
title("相频响应");
axis([0,1,-1,1]);
grid;
set(gca,'XTickMode','manual','XTick',[0,0.15,0.2,0.6,0.65,1]);
set(gca,'YTickMode','manual','YTick',[-1,0,1]);
subplot(223);
plot(w1/pi,db1);
ylabel('dB');
xlabel('频率（\pi）');
axis([0,1,-40,5]);
set(gca,'XTickMode','manual','XTick',[0,0.15,0.2,0.6,0.65,1]);
set(gca,'YTickMode','manual','YTick',[-70,-45,-1,0]);
axis([0,1.1,-80,5]); 
grid;
subplot(224);
subplot(2,2,4);
zplane(b1,a1);
axis([-1.1,1.1,-1.1,1.1]); 
title('零极点图');


%% 思考题4
omega_p = [0.15 * pi,0.65 * pi];    %滤波器的通带截止频率
omega_s = [0.2 * pi,0.6 * pi]; %滤波器的通带截止频率
Rp = 1;             %滤波器的通带衰减指标
As = 45;            %滤波器的阻带衰减指标
ripple = 10^(-Rp/20);   %滤波器的通带衰减对应的幅度值
Attn = 10^(-As/20);     %滤波器的阻带衰减对应的幅度值

Fs = 100;
T = 1 / Fs;
Omgp = tan(omega_p / 2) * Fs * 2;
Omgs = tan(omega_s / 2) * Fs * 2;
[n,Omgc] = buttord(Omgp,Omgs,Rp,As,'s');    %计算阶数n和截止频率
[ba,aa] = butter(n,Omgc,'stop','s');
[bd,ad] = bilinear(ba,aa,Fs);       %利用双线性变换法将模拟滤波器转换为数字滤波器
[H,w] = freqz(bd,ad);       %求数字系统的频率特性
dbH = 20 * log10((abs(H) + eps)/max(abs(H)));
figure(9);
subplot(2,2,1);
plot(w / pi,abs(H));
ylabel('|H|');
title("幅频响应");
axis([0,1,0,1.1]);
set(gca,'XTickMode','manual','XTick',[0,0.15,0.2,0.6,0.65,1]);
set(gca,'YTickMode','manual','YTick',[0,Attn,ripple,1]);
grid;
subplot(2,2,2);
plot(w / pi,angle(H)/pi);
ylabel('\phi');
title("相频响应");
axis([0,1,-1,1]);
set(gca,'XTickMode','manual','XTick',[0,0.15,0.2,0.6,0.65,1]);
set(gca,'YTickMode','manual','YTick',[-1,0,1]);
grid;
subplot(2,2,3);
plot(w / pi,dbH);
ylabel('dB');
xlabel('频率（\pi）');
title("幅频响应(dB)");
axis([0,1,-40,5]);
set(gca,'XTickMode','manual','XTick',[0,0.15,0.2,0.6,0.65,1]);
set(gca,'YTickMode','manual','YTick',[-70,-As,-Rp,0]);
axis([0,1,-70,5]); 
grid;
subplot(2,2,4);
zplane(bd,ad);
axis([-1.1,1.1,-1.1,1.1]); 
title('零极点图');


%% 额外
omega_p = [0.15 * pi,0.65 * pi];    %滤波器的通带截止频率
omega_s = [0.2 * pi,0.6 * pi]; %滤波器的通带截止频率
Rp = 1;             %滤波器的通带衰减指标
As = 45;            %滤波器的阻带衰减指标
ripple = 10^(-Rp/20);   %滤波器的通带衰减对应的幅度值
Attn = 10^(-As/20);     %滤波器的阻带衰减对应的幅度值

Fs = 100;
T = 1 / Fs;
Omgp = tan(omega_p / 2) * Fs * 2;
Omgs = tan(omega_s / 2) * Fs * 2;
[n,Omgc] = buttord(Omgp,Omgs,Rp,As,'s');    %计算阶数n和截止频率
[ba,aa] = butter(n,Omgc,'stop','s');
[bd,ad] = bilinear(ba,aa,Fs);       %利用双线性变换法将模拟滤波器转换为数字滤波器
[H,w] = freqz(bd,ad);       %求数字系统的频率特性
dbH = 20 * log10((abs(H) + eps)/max(abs(H)));
figure(12);
subplot(2,2,1);
plot(w / pi,abs(H));
ylabel('|H|');
title("幅频响应");
axis([0,1,0,1.1]);
set(gca,'XTickMode','manual','XTick',[0,0.15,0.2,0.6,0.65,1]);
set(gca,'YTickMode','manual','YTick',[0,Attn,ripple,1]);
grid;
subplot(2,2,2);
plot(w / pi,angle(H)/pi);
ylabel('\phi');
title("相频响应");
axis([0,1,-1,1]);
set(gca,'XTickMode','manual','XTick',[0,0.15,0.2,0.6,0.65,1]);
set(gca,'YTickMode','manual','YTick',[-1,0,1]);
grid;
subplot(2,2,3);
plot(w / pi,dbH);
ylabel('dB');
xlabel('频率（\pi）');
title("幅频响应(dB)");
axis([0,1,-40,5]);
set(gca,'XTickMode','manual','XTick',[0,0.15,0.2,0.6,0.65,1]);
set(gca,'YTickMode','manual','YTick',[-70,-As,-Rp,0]);
axis([0,1,-70,5]); 
grid;
subplot(2,2,4);
zplane(bd,ad);
axis([-1.1,1.1,-1.1,1.1]); 
title('零极点图');
