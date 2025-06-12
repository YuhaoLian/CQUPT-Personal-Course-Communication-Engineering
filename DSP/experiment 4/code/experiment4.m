clc;
clear;
%% 任务一
N=8;
for i=0:N-1
    x1(i+1)=X1(i);
end
xk1=fft(x1,N);
figure(1);
subplot(211); 
stem(0:length(x1)-1,x1,'.');
title('x1 的波形');
subplot(212); 
stem(0:N-1,abs(xk1),'.');
title('x1 的 8 点离散幅度谱');


N=16;
xk1=fft(x1,N);
figure(2);
subplot(211); 
stem(0:length(x1)-1,x1,'.');
title('x1 的波形');
subplot(212); 
stem(0:N-1,abs(xk1),'.');
title('x1 的 16 点离散幅度谱');


%% 任务二
N=8;
for i=0:N-1
    x2(i+1)=X2(i);
end
xk2=fft(x2,N);
figure(3);
subplot(211); 
stem(0:length(x2)-1,x2,'.');
title('x2 的波形');
subplot(212); 
stem(0:N-1,abs(xk2),'.');
title('x2 的 8 点离散幅度谱');


N=16;
xk2=fft(x2,N);
figure(4);
subplot(211); 
stem(0:length(x2)-1,x2,'.');
title('x2 的波形');
subplot(212); 
stem(0:N-1,abs(xk2),'.');
title('x2 的 16 点离散幅度谱');

%% 任务三
N=8;
for i=0:N-1
    x4(i+1)=X4(i);
end
xk4=fft(x4,N)
figure(5);
subplot(211); 
stem(0:length(x4)-1,x4,'.');
title('x4 的波形');
subplot(212); 
stem(0:N-1,abs(xk4),'.');
title('x4 的 8 点离散幅度谱');


N=16;
for i=0:N-1
    x4(i+1)=X4(i);
end
xk4=fft(x4,N)
figure(6);
subplot(211); 
stem(0:length(x4)-1,x4,'.');
title('x4 的波形');
subplot(212); 
stem(0:N-1,abs(xk4),'.');
title('x4 的 16 点离散幅度谱');

%% 任务四
Fs = 64; %采样频率64hz
N = 16;
t = 0:1/Fs:(N-1)/Fs;
for i =1:length(t)
    x6(i)=X6(t(i));
end
xk6=fft(x6,N);
figure(7);
subplot(211); 
stem(0:length(x6)-1,x6,'.');
title('x6 的波形');
subplot(212); 
stem(0:N-1,abs(xk6),'.');
title('x6 的 16 点离散幅度谱');

N = 32;
t = 0:1/Fs:(N-1)/Fs;
for i =1:length(t)
    x6(i)=X6(t(i));
end
xk6=fft(x6,N);
figure(8);
subplot(211); 
stem(0:length(x6)-1,x6,'.');
title('x6 的波形');
subplot(212); 
stem(0:N-1,abs(xk6),'.');
title('x6 的 32 点离散幅度谱');

N = 64;
t = 0:1/Fs:(N-1)/Fs;
for i =1:length(t)
    x6(i)=X6(t(i));
end
xk6=fft(x6,N);
figure(9);
subplot(211); 
stem(0:length(x6)-1,x6,'.');
title('x6 的波形');
subplot(212); 
stem(0:N-1,abs(xk6),'.');
title('x6 的 64 点离散幅度谱');


%% 任务五
N=512;
[xn,fs] = audioread('F:\学习\数字信号处理实验\实验1\motherland.wav');
% n = 0:N-1;
% t = n/fs;  
% f1 = n*fs/N-fs/2;
w1 = 0:2*pi/(N-1):2*pi;
xn = xn(8000:8199);
xnk = fft(xn,N)./fs;
Hm = abs(xnk);
Hp = angle(xnk);
figure(10);
subplot(311)
plot(1:1/8000:1+199/8000,xn);
grid on;
xlabel("t");
ylabel("x(n)");
title("原信号采样序列");
subplot(312); 
plot(w1/pi,Hm);
grid on;
xlabel("ω/Π(rad/s)");
ylabel("Magnitude");
title("幅频特性曲线");
subplot(313); 
plot(w1/pi,Hp);
grid on;
xlabel("ω/Π(rad/s)");
ylabel("Phase");
title("相频特性曲线");
%% 图像处理
% lena = imread('lena_colour.bmp'); % 读原图
lena = imread('lena.bmp'); % 读原图
figure(13); 
imshow(lena);
title('原图') 
fftI = fft2(lena); % 二维离散傅里叶变换
A1 = abs(fftI); % 取模值
% 把幅度限定在[0,255]
B1=(A1-min(min(A1)))./(max(max(A1))-min(min(A1)))*255;
figure(14); 
imshow(B1); title('二维幅度谱图');
B= fftshift(B1); 
figure(15);
imshow(B); title('移到中心位置的二维频谱图')

%% 思考题1
N=8;
for i=0:N-1
    x2(i+1)=X2(i);
end
xk2=fft(x2,N);
figure(11);
subplot(321); 
stem(0:length(x2)-1,x2,'.');
title('x2 的波形');
subplot(323); 
stem(0:N-1,abs(xk2),'.');
title('x2 的 8 点离散幅度谱');

N=16;
xk2=fft(x2,N);
subplot(325); 
stem(0:N-1,abs(xk2),'.');
title('x2 的 16 点离散幅度谱');

N=8;
for i=0:N-1
    x3(i+1)=X3(i);
end
xk3=fft(x3,N);
subplot(322); 
stem(0:length(x3)-1,x3,'.');
title('x3 的波形');
subplot(324); 
stem(0:N-1,abs(xk3),'.');
title('x3 的 8 点离散幅度谱');
N=16;
xk3=fft(x3,N);
subplot(326); 
stem(0:N-1,abs(xk3),'.');
title('x3 的 16 点离散幅度谱');

%% 思考题3
xx = [1 2 2 3 3 2 2 1];
xx1 = [1 2 3 2 1];
xx2 = [1 0 2 0 2 0 3 0 3 0 2 0 2 0 1];
wx = 0:2*pi/511:2*pi;
xxk = fft(xx,512);
HHm = abs(xxk);
HHp = angle(xxk);
xxk1 = fft(xx1,512);
HHm1 = abs(xxk1);
HHp1 = angle(xxk1);
xxk2 = fft(xx2,512);
HHm2 = abs(xxk2);
HHp2 = angle(xxk2);
figure(12);
subplot(331); 
[n,m]=size(xx);
stem(0:m-1,xx,'.');
xlabel("n");
title("x(n)");
subplot(334); 
plot(wx/pi,HHm);
xlabel("\omega/\pi");
title("幅频特性曲线");
subplot(337); 
plot(wx/pi,HHp);
xlabel("\omega/\pi");
title("相频特性曲线");
subplot(332); 
[n,m]=size(xx1);
stem(0:m-1,xx1,'.');
xlabel("n");
title("x1(n)");
subplot(335); 
plot(wx/pi,HHm1);
xlabel("\omega/\pi");
title("辐频特性曲线");
subplot(338); 
plot(wx/pi,HHp1);
xlabel("\omega/\pi");
title("相频特性曲线");
subplot(333); 
[n,m]=size(xx2);
stem(0:m-1,xx2,'.');
xlabel("n");
title("x2(n)");
subplot(336); 
plot(wx/pi,HHm2);
xlabel("\omega/\pi");
title("辐频特性曲线");
subplot(339); 
plot(wx/pi,HHp2);
xlabel("\omega");
title("相频特性曲线");
