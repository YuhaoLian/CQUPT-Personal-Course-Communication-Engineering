clc;
clear;
close all;
%% 窗函数特性的比较
N = 64;
beta = 7.865;
n = 1:N;
wbo = boxcar(N);    %矩形窗
wtr = triang(N);    %三角窗
whn = hanning(N);   %汉宁窗
whm = hamming(N);   %哈明窗
wbl = blackman(N);  %布莱克曼窗
wka = kaiser(N,beta);   %凯塞窗
plot(n,wbo,'-',n,wtr,'*',n,whn,'+',n,whm,'.',n,wbl,'o',n,wka,'d');
axis([0,N,0,1.1]);
legend('矩形','三角形','汉宁','哈明','布莱克曼','凯塞');

%% 用窗函数法设计 FIR 数字滤波器
wp = 0.3 * pi;
ws = 0.45 * pi;
delta_w = ws - wp;
N0 = (6.6 * pi / delta_w);
N = N0 + mod(N0 + 1,2);     %实现FIR类型偶对称滤波器，应确保N为奇数
windows = (hamming(N))'; 
wc = (ws + wp)/2;
hd = ideal_lp(wc,N);
b = hd.*windows;
[H,w] = freqz(b,1,1000,'whole');
H = (H(1:501))';
w = (w(1:501))';
mag = abs(H);
db = 20 * log10((mag + eps) / max(mag));
pha = angle(H)
n = 0:N-1;
dw = 2 * pi / 1000;
Rp = -(min(db(1:wp/dw+1))); %检验通带波动
As = -round(max(db(ws/dw+1:501)));  %检验最小阻带衰减
figure(2);
subplot(2,2,1);
stem(n,b,'.');
title('实际脉冲响应');
xlabel('n');
ylabel('h(n)');
subplot(2,2,2);
stem(n,windows,'.');
axis([0,N,0,1.1]);
title('窗函数特性');
xlabel('n');
ylabel('wd(n)');
subplot(2,2,3);
plot(w/pi,db);
axis([0,1,-80,10]);
title('幅度频率响应');
xlabel('频率（单位：\pi）');
ylabel('H(e^{j\omega})');
set(gca,'XTickMode','manual','XTick',[0,wp/pi,ws/pi,1]);
set(gca,'YTickMode','manual','YTick',[-50,-20,-3,0]);
grid;
subplot(2,2,4);
plot(w/pi,pha);
axis([0,1,-4,4]);
title('相位频率响应');
xlabel('频率（单位：\pi）');
ylabel('\phi(\omega)');
set(gca,'XTickMode','manual','XTick',[0,wp/pi,ws/pi,1]);
set(gca,'YTickMode','manual','YTick',[-3.1416,0.3,3.1416,4]);
grid;



%% 利用工具箱设计
clc;
clear;
wp = 0.3 * pi;
ws = 0.45 * pi;
delta_w = ws - wp;
N0 = ceil(6.6 * pi / delta_w);
N = N0 + mod(N0 + 1,2);     %实现FIR类型偶对称滤波器，应确保N为奇数
windows = (hamming(N))'; 
wc = (ws + wp)/2/pi;
b = fir1(N-1,wc,windows);
[H,w] = freqz(b,1,1000,'whole');
H = (H(1:501))';
w = (w(1:501))';
mag = abs(H);
db = 20 * log10((mag + eps) / max(mag));
pha = angle(H);
n = 0:N-1;
dw = 2 * pi / 1000;
Rp = -(min(db(1:wp/dw+1))); %检验通带波动
As = -round(max(db(ws/dw+1:501)));  %检验最小阻带衰减

%% 设计fir低通滤波器
clc;
clear;
wp = 0.24 * pi;
ws = 0.3 * pi;
delta_w = ws - wp;
N0 = ceil(11 * pi / delta_w);
N = N0 + mod(N0 + 1,2);     %实现FIR类型偶对称滤波器，应确保N为奇数
windows = (blackman(N))'; 
wc = (ws + wp)/2;
hd = ideal_lp(wc,N);
b = hd.*windows;
[H,w] = freqz(b,1,1000,'whole');
H = (H(1:501))';
w = (w(1:501))';
mag = abs(H);
db = 20 * log10((mag + eps) / max(mag));
pha = angle(H);
n = 0:N-1;
dw = 2 * pi / 1000;
Rp = -(min(db(1:int32(wp/dw+1)))); %检验通带波动
As = -round(max(db(ws/dw+1:501)));  %检验最小阻带衰减
figure(3);
subplot(2,2,1);
stem(n,b,'.');
title('实际脉冲响应');
xlabel('n');
ylabel('h(n)');
subplot(2,2,2);
stem(n,windows,'.');
axis([0,N,0,1.1]);
title('窗函数特性');
xlabel('n');
ylabel('wd(n)');
subplot(2,2,3);
plot(w/pi,db);
axis([0,1,-100,10]);
title('幅度频率响应');
xlabel('频率（单位：\pi）');
ylabel('H(e^{j\omega})');
set(gca,'XTickMode','manual','XTick',[0,wp/pi,ws/pi,1]);
set(gca,'YTickMode','manual','YTick',[-60,-20,-3,0]);
grid on;
subplot(2,2,4);
plot(w/pi,pha);
axis([0,1,-4,4]);
title('相位频率响应');
xlabel('频率（单位：\pi）');
ylabel('\phi(\omega)');
set(gca,'XTickMode','manual','XTick',[0,wp/pi,ws/pi,1]);
set(gca,'YTickMode','manual','YTick',[-3.1416,0.3,3.1416,4]);
grid;

%% 产生加性噪声
figure(4);
[xt,t] = xtg(1000);

%% 对噪声滤波
yt = fftfilt(b,xt);
Hyk=abs(fft(yt));
figure(5);
subplot(2,1,1)
plot(t,yt)
title('时域波形图')
subplot(2,1,2)
stem(Hyk,'.')
axis([80,120,min(Hyk),max(Hyk)])

%% 选做
% 内插0
clc;
clear;
[xn, fs] = audioread('motherland.wav');
I = 2;  % 实现I=2的整数倍0值内插
yn1 = [];   %预分配内存
for i=1:length(xn)
    yn1(I*i-1)=xn(i);
    yn1(I*i)= 0;
end
% sound(yn1,I*fs);

%滤波器设计
wp = 0.45 * pi;
ws = 0.55 * pi;
delta_w = ws - wp;
N0 = ceil(11 * pi / delta_w);
N = N0 + mod(N0 + 1,2);     %实现FIR类型偶对称滤波器，应确保N为奇数
windows = (blackman(N))'; 
wc = (ws + wp)/2;
hd = ideal_lp(wc,N);
b = hd.*windows;
[H,w] = freqz(b,1,1000,'whole');
H = (H(1:501))';
w = (w(1:501))';
mag = abs(H);
db = 20 * log10((mag + eps) / max(mag));
pha = angle(H);
n = 0:N-1;
dw = 2 * pi / 1000;
Rp = -(min(db(1:int32(wp/dw+1)))); %检验通带波动
As = -round(max(db(ws/dw+1:501)));  %检验最小阻带衰减
figure(6);
subplot(2,2,1);
stem(n,b,'.');
title('实际脉冲响应');
xlabel('n');
ylabel('h(n)');
subplot(2,2,2);
stem(n,windows,'.');
axis([0,N,0,1.1]);
title('窗函数特性');
xlabel('n');
ylabel('wd(n)');
subplot(2,2,3);
plot(w/pi,db);
axis([0,1,-100,10]);
title('幅度频率响应');
xlabel('频率（单位：\pi）');
ylabel('H(e^{j\omega})');
set(gca,'XTickMode','manual','XTick',[0,wp/pi,ws/pi,1]);
set(gca,'YTickMode','manual','YTick',[-60,-20,-0.1,0]);
grid on;
subplot(2,2,4);
plot(w/pi,pha);
axis([0,1,-4,4]);
title('相位频率响应');
xlabel('频率（单位：\pi）');
ylabel('\phi(\omega)');
set(gca,'XTickMode','manual','XTick',[0,wp/pi,ws/pi,1]);
set(gca,'YTickMode','manual','YTick',[-3.1416,0.3,3.1416,4]);
grid;

%滤波并画图
yn2 = fftfilt(b,yn1);
sound(yn2,I*fs);
figure(7);
subplot(4,1,1)
N=2048;
Xn=1/fs*fft(xn(8000:8199),N); % 从xn中取200点做谱分析，N可取2048
plot((0:N/2-1)*fs/N,abs(Xn(1:N/2))); % 模拟域幅度谱
xlabel('f(Hz)')
title('原始信号模拟域幅度谱')
subplot(4,1,2)
Yn1=1/(I*fs)*fft(yn1(16000:16399),N); % 内插后，200点长变成了400点长
plot((0:N/2-1)*I*fs/N,abs(Yn1(1:N/2)));
xlabel('f(Hz)')
title('I=2内插处理后的信号模拟域幅度谱')
subplot(4,1,3)
yn2=filter(b,1,yn1); % 对yn1进行滤波，b为所设计的镜像滤波器
Yn2=1/(I*fs)*fft(yn2(16000:16399),N); % 内插后，200点长变成了400点长
plot((0:N/2-1)*I*fs/N,abs(Yn2(1:N/2)));
xlabel('f(Hz)')
title('信号内插后再经过镜像滤波后的模拟域幅度谱')
subplot(4,1,4)
plot((0:N/4-1)*I*fs/N,abs(Yn2(1:N/4)));
xlabel('f(Hz)')
title('信号内插后再经过镜像滤波后的模拟域幅度谱')

%% 思考题
%低通
clc;
clear;
wp = 0.2 * pi;
ws = 0.3 * pi;
delta_w = ws - wp;
N0 = ceil(6.1 * pi / delta_w);
N = N0 + mod(N0 + 1,2);     %实现FIR类型偶对称滤波器，应确保N为奇数
windows = (triang(N))'; 
wc = (ws + wp)/2/pi;
b = fir1(N-1,wc,windows);
[H,w] = freqz(b,1,1000,'whole');
H = (H(1:501))';
w = (w(1:501))';
mag = abs(H);
db = 20 * log10((mag + eps) / max(mag));
pha = angle(H);
n = 0:N-1;
dw = 2 * pi / 1000;
Rp = -(min(db(1:wp/dw+1))); %检验通带波动
As = -round(max(db(ws/dw+1:501)));  %检验最小阻带衰减
ripple = 10^(-Rp/20);   %滤波器的通带衰减对应的幅度值
Attn = 10^(-As/20);     %滤波器的阻带衰减对应的幅度值
figure(8);
subplot(2,2,1);
stem(n,b,'.');
title('实际脉冲响应');
xlabel('n');
ylabel('h(n)');
% subplot(2,2,2);
% stem(n,windows,'.');
% axis([0,N,0,1.1]);
% title('窗函数特性');
% xlabel('n');
% ylabel('wd(n)');
subplot(2,2,2);
plot(w/pi,mag);
title('幅度频率响应');
xlabel('频率（单位：\pi）');
ylabel('H(e^{j\omega})');
set(gca,'XTickMode','manual','XTick',[0,wp/pi,ws/pi,1]);
set(gca,'YTickMode','manual','YTick',[0,Attn,ripple,1]);
grid;
subplot(2,2,3);
plot(w/pi,db);
axis([0,1,-100,10]);
title('幅度频率响应');
xlabel('频率（单位：\pi）');
ylabel('db H(e^{j\omega})');
set(gca,'XTickMode','manual','XTick',[0,wp/pi,ws/pi,1]);
set(gca,'YTickMode','manual','YTick',[-40,-As,-1,0]);
axis([0,1,-50,5]);
grid on;
subplot(2,2,4);
plot(w/pi,pha);
axis([0,1,-4,4]);
title('相位频率响应');
xlabel('频率（单位：\pi）');
ylabel('\phi(\omega)');
set(gca,'XTickMode','manual','XTick',[0,wp/pi,ws/pi,1]);
set(gca,'YTickMode','manual','YTick',[-3.1416,0.3,3.1416,4]);
grid;

%高通
clc;
clear;
wp = 0.4 * pi;
ws = 0.6 * pi;
delta_w = ws - wp;
N0 = ceil(6.2 * pi / delta_w);
N = N0 + mod(N0 + 1,2);     %实现FIR类型偶对称滤波器，应确保N为奇数
windows = (hanning(N))'; 
wc = (ws + wp)/2/pi;
b = fir1(N-1,wc,'high',windows);
[H,w] = freqz(b,1,1000,'whole');
H = (H(1:501))';
w = (w(1:501))';
mag = abs(H);
db = 20 * log10((mag + eps) / max(mag));
pha = angle(H);
n = 0:N-1;
dw = 2 * pi / 1000;
Rp = -(min(db(1:wp/dw+1))); %检验通带波动
As = -round(max(db(ws/dw+1:501)));  %检验最小阻带衰减
ripple = 10^(-Rp/20);   %滤波器的通带衰减对应的幅度值
Attn = 10^(-As/20);     %滤波器的阻带衰减对应的幅度值
figure(9);
subplot(2,2,1);
stem(n,b,'.');
title('实际脉冲响应');
xlabel('n');
ylabel('h(n)');
subplot(2,2,2);
% stem(n,windows,'.');
% axis([0,N,0,1.1]);
% title('窗函数特性');
% xlabel('n');
% ylabel('wd(n)');
plot(w/pi,mag);
title('幅度频率响应');
xlabel('频率（单位：\pi）');
ylabel('H(e^{j\omega})');
set(gca,'XTickMode','manual','XTick',[0,wp/pi,ws/pi,1]);
set(gca,'YTickMode','manual','YTick',[0,ripple,Attn]);
grid;
subplot(2,2,3);
plot(w/pi,db);
axis([0,1,-100,10]);
title('幅度频率响应');
xlabel('频率（单位：\pi）');
ylabel('db H(e^{j\omega})');
set(gca,'XTickMode','manual','XTick',[0,wp/pi,ws/pi,1]);
set(gca,'YTickMode','manual','YTick',[-70,-43,-0.2,0]);
grid on;
subplot(2,2,4);
plot(w/pi,pha);
axis([0,1,-4,4]);
title('相位频率响应');
xlabel('频率（单位：\pi）');
ylabel('\phi(\omega)');
set(gca,'XTickMode','manual','XTick',[0,wp/pi,ws/pi,1]);
set(gca,'YTickMode','manual','YTick',[-3.1416,0.3,3.1416,4]);
grid;

%带通
clc;
clear;
wpl = 0.15 * pi;
wph = 0.65 * pi;
wsl = 0.2 * pi;
wsh = 0.6 * pi;
wd1 = (wpl + wsl)/2;
wd2 = (wph + wsh)/2;
delta_w = min(abs(wsl-wpl),abs(wsh-wph));
N0 = ceil(6.6 * pi / delta_w);%%8好一点
N = N0 + mod(N0 + 1,2);     %实现FIR类型偶对称滤波器，应确保N为奇数
wn = [wd1/pi,wd2/pi];
b = fir1(N-1,wn,'bandpass', hamming(N));

[H,w] = freqz(b,1,1000,'whole');
H = (H(1:501))';
w = (w(1:501))';
mag = abs(H);
db = 20 * log10((mag + eps) / max(mag));
pha = angle(H);
n = 0:N-1;
dw = 2 * pi / 1000;
Rp = -(min(db(1:(wph-wpl)/dw+1))); %检验通带波动
As = -round(max(db((wsh-wsl)/dw+1:501)));  %检验最小阻带衰减
ripple = 10^(-Rp/20);   %滤波器的通带衰减对应的幅度值
Attn = 10^(-As/20);     %滤波器的阻带衰减对应的幅度值
figure(10);
subplot(2,2,1);
stem(n,b,'.');
title('实际脉冲响应');
xlabel('n');
ylabel('h(n)');
subplot(2,2,2);
% stem(n,windows,'.');
% axis([0,N,0,1.1]);
% title('窗函数特性');
% xlabel('n');
% ylabel('wd(n)');
plot(w/pi,mag);
title('幅度频率响应');
xlabel('频率（单位：\pi）');
ylabel('H(e^{j\omega})');
set(gca,'XTickMode','manual','XTick',[0,wpl/pi,wsl/pi,wsh/pi,wph/pi,1]);
set(gca,'YTickMode','manual','YTick',[0,ripple,Attn]);
grid;
subplot(2,2,3);
plot(w/pi,db);
axis([0,1,-100,10]);
title('幅度频率响应');
xlabel('频率（单位：\pi）');
ylabel('db H(e^{j\omega})');
set(gca,'XTickMode','manual','XTick',[0,wpl/pi,wsl/pi,wsh/pi,wph/pi,1]);
set(gca,'YTickMode','manual','YTick',[-70,-50,-1,0]);
grid on;
subplot(2,2,4);
plot(w/pi,pha);
axis([0,1,-4,4]);
title('相位频率响应');
xlabel('频率（单位：\pi）');
ylabel('\phi(\omega)');
set(gca,'XTickMode','manual','XTick',[0,wpl/pi,wsl/pi,wsh/pi,wph/pi,1]);
set(gca,'YTickMode','manual','YTick',[-3.1416,0.3,3.1416,4]);
grid;


%带阻
clc;
clear;
wpl = 0.15 * pi;
wph = 0.65 * pi;
wsl = 0.2 * pi;
wsh = 0.6 * pi;
wd1 = (wpl + wsl)/2;
wd2 = (wph + wsh)/2;
delta_w = min(abs(wsl-wpl),abs(wsh-wph));
N0 = ceil(6.2 * pi / delta_w);
N = N0 + mod(N0 + 1,2);     %实现FIR类型偶对称滤波器，应确保N为奇数
wn = [wd1/pi,wd2/pi];
b = fir1(N-1,wn,'stop', hanning(N));

[H,w] = freqz(b,1,1000,'whole');
H = (H(1:501))';
w = (w(1:501))';
mag = abs(H);
db = 20 * log10((mag + eps) / max(mag));
pha = angle(H);
n = 0:N-1;
dw = 2 * pi / 1000;
Rp = -(min(db(1:(wph-wpl)/dw+1))); %检验通带波动
As = -round(max(db((wsh-wsl)/dw+1:501)));  %检验最小阻带衰减
ripple = 10^(-Rp/20);   %滤波器的通带衰减对应的幅度值
Attn = 10^(-As/20);     %滤波器的阻带衰减对应的幅度值
figure(11);
subplot(2,2,1);
stem(n,b,'.');
title('实际脉冲响应');
xlabel('n');
ylabel('h(n)');
subplot(2,2,2);
% stem(n,windows,'.');
% axis([0,N,0,1.1]);
% title('窗函数特性');
% xlabel('n');
% ylabel('wd(n)');
plot(w/pi,mag);
title('幅度频率响应');
xlabel('频率（单位：\pi）');
ylabel('H(e^{j\omega})');
set(gca,'XTickMode','manual','XTick',[0,wpl/pi,wsl/pi,wsh/pi,wph/pi,1]);
set(gca,'YTickMode','manual','YTick',[0,ripple,Attn]);
grid;
subplot(2,2,3);
plot(w/pi,db);
axis([0,1,-100,10]);
title('幅度频率响应');
xlabel('频率（单位：\pi）');
ylabel('db H(e^{j\omega})');
set(gca,'XTickMode','manual','XTick',[0,wpl/pi,wsl/pi,wsh/pi,wph/pi,1]);
set(gca,'YTickMode','manual','YTick',[-70,-45,-1,0]);
grid on;
subplot(2,2,4);
plot(w/pi,pha);
axis([0,1,-4,4]);
title('相位频率响应');
xlabel('频率（单位：\pi）');
ylabel('\phi(\omega)');
set(gca,'XTickMode','manual','XTick',[0,wpl/pi,wsl/pi,wsh/pi,wph/pi,1]);
set(gca,'YTickMode','manual','YTick',[-3.1416,0.3,3.1416,4]);
grid;