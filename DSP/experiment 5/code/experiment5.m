clc;
clear;
%% 时域采样定理的验证example
x = [1 1 1 1];
N = 64;
xk = fft(x,N);
Xm = abs(xk);
Xp = angle(xk);
figure(1);
subplot(211);
stem(0:3,x,'.');
subplot(212);
k = 0:N-1;
plot(2*k/N,Xm);
title("|FT[x(n)]|");
xlabel("\omega/\pi")

%% 时域采样定理的验证
% 1khz采样
t = 0:63;
n = t/1000;
f = 444.128*exp(-50*2^0.5.*n*pi).*sin(50*2^0.5.*n*pi).*uDT(n);
figure(2);
subplot(321);
stem(t,f,'.');
title("Fs=1khz");
fk = fft(f)/1000;
fm = abs(fk);
fp = angle(fk);
subplot(322);
plot(0:1000/63:1000,fm);
title("|T*FT[x1(n)]|");
xlabel("f(hz)");

%300hz采样
n = 0:1/300:0.063;
f = 444.128*exp(-50*2^0.5.*n*pi).*sin(50*2^0.5.*n*pi).*uDT(n);
subplot(323);
stem(0:18,f,'.');
title("Fs=300hz");
fk = fft(f)/300;
fm = abs(fk);
fp = angle(fk);
subplot(324);
plot(0:300/18:300,fm);
title("|T*FT[x2(n)]|");
xlabel("f(hz)");

%200hz采样
n = 0:1/200:0.063;
f = 444.128*exp(-50*2^0.5.*n*pi).*sin(50*2^0.5.*n*pi).*uDT(n);
subplot(325);
stem(0:12,f,'.');
title("Fs=200hz");
fk = fft(f)/200;
fm = abs(fk);
fp = angle(fk);
subplot(326);
plot(0:300/12:300,fm);
title("|T*FT[x3(n)]|");
xlabel("f(hz)");

%% 频域采样定理的验证
clc;
clear;
n = 0:26;
for i = 0:26
    y(i+1) = x(i);
end
yk = fft(y,512);
ym = abs(yk);
yp = angle(yk);
figure(3);
subplot(321);
plot(0:2*pi/511:2*pi,ym);
title("连续频谱|FT[x(n)]|");
xlabel("\omega");
subplot(322);
stem(n,y,'.');
title("三角序列");
xlabel("n");
% 16点采样
yk = fft(y,16);
ym = abs(yk);
yp = angle(yk);
subplot(323);
stem(0:15,ym,'.');
title("16点频谱|FT[x(n)]|");
xlabel("k");
subplot(324);
yy = ifft(yk,16);
stem(0:15,yy,'.');
title("16点三角序列");
xlabel("n");
% 32点采样
yk = fft(y,32);
ym = abs(yk);
yp = angle(yk);
subplot(325);
stem(0:31,ym,'.');
title("32点频谱|FT[x(n)]|");
xlabel("k");
subplot(326);
yy = ifft(yk,32);
stem(0:31,yy,'.');
title("32点三角序列");
xlabel("n");


%%
M=27;N=32;n=0:M;
%产生 M 长三角波序列 x(n)
xa=0:floor(M/2); xb= ceil(M/2)-1:-1:0; xn=[xa,xb];
for i = 1:length(xn)
    if i<=13
        xn(i) = xn(i) +1;
    else
        xn(i) = xn(i) -1;
    end
end
x(14)=x(14)+1;
Xk=fft(xn,1024); %1024 点 FFT[x(n)], 用于近似序列 x(n)的 TF
X32k=fft(xn,32) ;%32 点 FFT[x(n)]
x32n=ifft(X32k); %32 点 IFFT[X32(k)]得到 x32(n)
X16k=X32k(1:2:N); %隔点抽取 X32k 得到 X16(K)
x16n=ifft(X16k,N/2); %16 点 IFFT[X16(k)]得到 x16(n)
subplot(3,2,2);stem(n,xn,'.');box on
title('(b) 三角波序列 x(n)');xlabel('n');ylabel('x(n)');axis([0,32,0,20])
k=0:1023;wk=2*k/1024; %
subplot(3,2,1);plot(wk,abs(Xk));title('(a)FT[x(n)]');
xlabel('\omega/\pi');ylabel('|X(e^j^\omega)|');axis([0,1,0,200])
k=0:N/2-1;
subplot(3,2,3);stem(k,abs(X16k),'.');box on
title('(c) 16 点频域采样');xlabel('k');ylabel('|X_1_6(k)|');axis([0,8,0,200])
n1=0:N/2-1;
subplot(3,2,4);stem(n1,x16n,'.');box on
title('(d) 16 点 IDFT[X_1_6(k)]');xlabel('n');ylabel('x_1_6(n)');axis([0,32,0,20])
k=0:N-1;
subplot(3,2,5);stem(k,abs(X32k),'.');box on
title('(e) 32 点频域采样');xlabel('k');ylabel('|X_3_2(k)|');axis([0,16,0,200])
n1=0:N-1;
subplot(3,2,6);stem(n1,x32n,'.');box on
title('(f) 32 点 IDFT[X_3_2(k)]');xlabel('n');ylabel('x_3_2(n)');axis([0,32,0,20])

%% 音频信号
[xn,fs] = audioread('F:\学习\数字信号处理实验\实验1\motherland.wav');
x = xn(1000:2999);
wx = 0:2*pi/511:2*pi;
xxk = fft(x,512);
HHm = abs(xxk);
HHp = angle(xxk);
figure(4);
subplot(311); 
[n,m]=size(x);
stem(0:n-1,x,'.');
title("2000个采样点原序列");
subplot(312); 
plot(wx,HHm);
xlabel("\omega");
title("幅频特性曲线");
subplot(313); 
plot(wx,HHp);
xlabel("\omega");
title("相频特性曲线");
for i = 1: 2000
    if mod(i,2)
    x1((i+1)/2) = x(i);
    else
    x2(i/2) = x(i);
    end
end
xxk = fft(x1,512);
HHm = abs(xxk);
HHp = angle(xxk);
figure(5);
subplot(311); 
[n,m]=size(x1);
stem(0:m-1,x1,'.');
title("奇数抽样序列");
subplot(312); 
plot(wx,HHm);
xlabel("\omega");
title("幅频特性曲线");
subplot(313); 
plot(wx,HHp);
xlabel("\omega");
title("相频特性曲线");


xxk = fft(x2,512);
HHm = abs(xxk);
HHp = angle(xxk);
figure(6);
subplot(311); 
[n,m]=size(x2);
stem(0:m-1,x2,'.');
title("偶数抽样序列");
subplot(312); 
plot(wx,HHm);
xlabel("\omega");
title("幅频特性曲线");
subplot(313); 
plot(wx,HHp);
xlabel("\omega");
title("相频特性曲线");

xk = fft(xn,1024);
figure
Hm = abs(xk);
Hp = angle(xk);
subplot(311); 
[n,m]=size(xn);
plot(0:n-1,xn);
title("原信号");
subplot(312); 
plot(0:2*pi/1023:2*pi,Hm);
xlabel("\omega");
title("幅频特性曲线");
subplot(313); 
plot(0:2*pi/1023:2*pi,Hp);
xlabel("\omega");
title("相频特性曲线");