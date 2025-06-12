clc;
clear;
%%
%第一问
A1 = [2 16 44 56 32];
B1 = [3 3 -15 18 -12];
[R1 P1 K1] = residuez(A1,B1);

%%
%第二问
A2 = [0 2 -1.6 -0.9];
B2 = [1 -2.5 1.96 -0.48];
A3 = [0 0 0 1 -1];
B3 = [1 -0.9 -0.65 0.873 0];
[R2 P2 K2] = tf2zp(A2,B2);
[R3 P3 K3] = tf2zp(A3,B3);
figure(1);
zplane(A2,B2);
figure(2);
zplane(A3,B3);

%%
%第三问
A4 = [1];
B4 = [1 -3/4 1/8];
[H w] = freqz(A4,B4,400,'whole');
Hm = abs(H);
Hp = angle(H);
figure(3);
subplot(2,1,1);
plot(w/pi,Hm);
grid on;
xlabel("ω/Π(rad/s)");
ylabel("Magnitude");
title("离散系统的幅频特性曲线");
subplot(2,1,2);
plot(w/pi,Hp);
grid on;
xlabel("ω/Π(rad/s)");
ylabel("Phase");
title("离散系统的相频特性曲线");

%% 思考题1
A5 = [1 0 0 0 0 0 0 0 -1];
B5 = [1 0 0 0 0 0 0 0 -0.9];
[R5 P5 K5] = tf2zp(A5,B5);
figure(4);
zplane(A5,B5);
[H1 w1] = freqz(A5,B5,400,'whole');
Hm1 = abs(H1);
Hp1 = angle(H1);
figure(5);
subplot(2,1,1);
plot(w1/pi,Hm1);
grid on;
xlabel("ω/Π(rad/s)");
ylabel("Magnitude");
title("离散系统的幅频特性曲线");
subplot(2,1,2);
plot(w1/pi,Hp1);
grid on;
xlabel("ω/Π(rad/s)");
ylabel("Phase");
title("离散系统的相频特性曲线");

%% 思考题2
[xn,fs] = audioread('F:\学习\数字信号处理实验\实验3\motherland.wav');
xn = xn(8000:16000);

A6 = [1 0];
B6 = [1 0.8];
[h1,t1] = impz(A6,B6);
[y1,zf1] = filter(A6,B6,xn);

A7 = [1 0];
B7 = [1 -1];
[h2,t2] = impz(A7,B7);
[y2,zf2] = filter(A7,B7,xn);

A8 = [1 0];
B8 = [1 1.2];
[h3,t3] = impz(A8,B8);
[y3,zf3] = filter(A8,B8,xn);

figure(6);
subplot(4,1,1);
stem(t1,h1);
xlabel("t");
ylabel("h1(t)");
title("原信号");
subplot(4,1,2);
[HH1 ww1] = freqz(A6,B6,400,'whole');
plot(w1/pi,abs(HH1));
grid on;
xlabel("ω/Π(rad/s)");
ylabel("Magnitude");
title("离散系统的幅频特性曲线");
subplot(4,1,4);
plot([1:1/8000:2],y1);
xlabel("t");
ylabel("f1(t)");
title("滤波后信号");
subplot(4,1,3);
plot([1:1/8000:2],xn);
xlabel("t");
ylabel("f1(t)");
title("原信号");

figure(7);
subplot(4,1,1);
stem(t2,h2);
xlabel("t");
ylabel("h2(t)");
title("原信号");
subplot(4,1,2);
[HH1 ww1] = freqz(A7,B7,400,'whole');
plot(w1/pi,abs(HH1));
grid on;
xlabel("ω/Π(rad/s)");
ylabel("Magnitude");
title("离散系统的幅频特性曲线");
subplot(4,1,4);
plot([1:1/8000:2],y2);
xlabel("t");
ylabel("f2(t)");
title("滤波后信号");
subplot(4,1,3);
plot([1:1/8000:2],xn);
xlabel("t");
ylabel("f1(t)");
title("原信号");

figure(8);
subplot(4,1,1);
stem(t3,h3);
xlabel("t");
ylabel("h3(t)");
title("原信号");
subplot(4,1,2);
[HH1 ww1] = freqz(A8,B8,400,'whole');
plot(w1/pi,abs(HH1));
grid on;
xlabel("ω/Π(rad/s)");
ylabel("Magnitude");
title("离散系统的幅频特性曲线");
subplot(4,1,4);
plot([1:1/8000:2],y3);
xlabel("t");
ylabel("f3(t)");
title("滤波后信号");
subplot(4,1,3);
plot([1:1/8000:2],xn);
xlabel("t");
ylabel("f1(t)");
title("原信号");

%% 思考题2
A = [1 -2 2];
B = [2 -2 1];
[R P K] = residuez(A,B);
figure(9);
subplot(3,1,1);
zplane(A,B);
[H w] = freqz(A,B,400,'whole');
Hm = abs(H);
Hp = angle(H);
subplot(3,1,2);
plot(w/pi,Hm);
grid on;
xlabel("ω/Π(rad/s)");
ylabel("Magnitude");
title("幅频特性曲线");
subplot(3,1,3);
plot(w/pi,Hp);
grid on;
xlabel("ω/Π(rad/s)");
ylabel("Phase");
title("相频特性曲线");


