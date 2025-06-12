clc;
clear;
%% 连续信号
t1 = -10:1/100:10;
t2 = -3:1/100:3;
t3 = 0:1/100:10;
Sa = sin(t1)./t1;
R = heaviside(t2+1)-heaviside(t2-1);
f = 5 * exp(0.5.* t3).* sin(2*pi.* t3);

figure(1)
subplot(3,1,1);
plot(t1,Sa);
title("Sa函数");
xlabel("t");
ylabel("Sa(t)");

subplot(3,1,2);
plot(t2,R);
title("门函数");
xlabel("t");
ylabel("R(t)");
axis([-3 3 -0.2 1.2]);

subplot(3,1,3);
plot(t3,f);
title("调制函数");
xlabel("t");
ylabel("f(t)");

%% 离散信号
k1 = -3:1:3;
k2 = -6:1:6;
k3 = 0:1:60;
delta = [0 0 0 1 0 0 0];
gk = heaviside(k2) - heaviside(k2-4);
for k = 1:1:13
    if gk(k)>0.2
        gk(k) = 1;
    end
end
F = (1.1.^k3).*sin(0.05*pi.*k3);

figure(2);
subplot(3,1,1);
stem(k1,delta);
title("冲激序列");
xlabel("k");
ylabel("δ(k)");

subplot(3,1,2);
stem(k2,gk);
title("门序列");
xlabel("k");
ylabel("g4(k)");

subplot(3,1,3);
stem(k3,F);
title("调制序列");
xlabel("k");
ylabel("F(k)");

%% 音频文件
[xn,fs] = audioread('F:\学习\数字信号处理实验\实验1\motherland.wav');
% sound(xn,fs);
K=400:1:500;

[x_n,f_s] = audioread('F:\学习\数字信号处理实验\实验1\motherland.wav',[fs,2*fs]);
T=1:1/8000:2;

figure(3);
subplot(2,1,1);
stem(K,xn(400:500));
xlabel("k");
ylabel("音频函数");
title("400-500音频采样");
subplot(2,1,2);
% stem(T,x_n);
plot(T,x_n,'r-');
xlabel("t");
ylabel("1-2s音频函数");
title("1-2s音频采样");


%% 图片文件
lena = imread('lena_colour.bmp'); 
lena_j = rgb2gray(lena);                
figure(4);                    
imshow(lena_j);
imwrite(lena_j,'lena.bmp'); 
