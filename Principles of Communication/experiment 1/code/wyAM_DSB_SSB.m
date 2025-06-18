clear all;
close all;
clc;

%% AM、DSB、SSB的调制及其解调

% fs = 128;  %抽样频率有点小
fs = 99999;
fm = 1000;
fc = 10000;
dt = 1/fs;
t = 0: dt : 50;
NFFT = 2^16; %傅里叶变换点数
f = linspace(-fs/2, fs/2, NFFT); %生成-fs/2到 fs/2的NFFT个点
iter = 1:4;
phase = [pi/8 pi/4 pi/3 pi/2];


%% 生成调制信号、载波信号以及已调信号
s_m = cos(2*pi*fm*t); %m(t)调制信号
s_x = sin(2*pi*fm*t); %m(t)调制信号希尔伯特变换
s_c = cos(2*pi*fc*t); %载波
c_q = sin(2*pi*fc*t); %载波希尔波特变换
s_AM = (s_m + 1).*s_c; %AM调制信号
s_AM_envelope = s_m+ 1; %AM调制信号包络
s_DSB = s_m .*s_c; %DSB调制信号
s_DSB_envelope = abs(s_m); %DSB调制信号包络
s_SSB_up = 0.5 * s_m .*s_c - 0.5 .*s_x .*c_q; %SSB_up调制信号
s_SSB_down = 0.5 * s_m .*s_c + 0.5 .*s_x .*c_q; %SSB_down调制信号

S_m = t2f(s_m,NFFT); %m(t)调制信号频谱
S_c = t2f(s_c,NFFT); %载波信号频谱
S_AM = t2f(s_AM,NFFT); %AM调制信号频谱
S_DSB = t2f(s_DSB,NFFT); %DSB调制信号频谱
S_SSB_up = t2f(s_SSB_up,NFFT); %SSB_up调制信号频谱
S_SSB_down = t2f(s_SSB_down,NFFT); %SSB_down调制信号频谱

%% 低通滤波
wc=3*2*pi*fm/fs;                        %解调器低通滤波器截止频率，3fm
lbf_hn=fir1(32,wc/pi);                     %低通滤波器系统函数，FIR滤波器
freqz(lbf_hn, 1, 1024, fs);

%% 相干解调
% AM
s_cor_temp_AM = s_AM .* 2 .* s_c; %通过乘法器
dem_cor_out_AM = filter(lbf_hn,1,s_cor_temp_AM) - 1; %通过低通滤波器
S_cor_temp_AM = t2f(s_cor_temp_AM,NFFT); 
S_cor_out_AM = t2f(dem_cor_out_AM,NFFT); 

% DSB
s_cor_temp_DSB = s_DSB.* 2 .* s_c; %通过乘法器
dem_cor_out_DSB = filter(lbf_hn,1,s_cor_temp_DSB); %通过低通滤波器
S_cor_temp_DSB = t2f(s_cor_temp_DSB,NFFT); 
S_cor_out_DSB = t2f(dem_cor_out_DSB,NFFT); 

% SSB_up
s_cor_temp_SSB_up = s_SSB_up .* 2 .* s_c; %通过乘法器
dem_cor_out_SSB_up= filter(lbf_hn,1,s_cor_temp_SSB_up); %通过低通滤波器
S_cor_temp_SSB_up = t2f(s_cor_temp_SSB_up,NFFT); 
S_cor_out_SSB_up = t2f(dem_cor_out_SSB_up,NFFT); 

% SSB_down
s_cor_temp_SSB_down = s_SSB_down .* 2 .* s_c; %通过乘法器
dem_cor_out_SSB_down= filter(lbf_hn,1,s_cor_temp_SSB_down); %通过低通滤波器
S_cor_temp_SSB_down = t2f(s_cor_temp_SSB_down,NFFT); 
S_cor_out_SSB_down = t2f(dem_cor_out_SSB_down,NFFT); 


%% 需要同学自己编写非相干解调
%包络检波
fbe=[0 0.05 0.1 1];
damps=[1 1 0 0];
b = firpm(100,fbe,damps);
[H_filter,W_filter]=freqz(b,1,512);
s_uncor_temp_AM = s_AM .* s_c;
dem_uncor_out_AM = (pi/2)*filter(b,1,abs(s_AM )) - 1;
S_uncor_temp_AM = t2f(s_uncor_temp_AM,NFFT); 
S_uncor_out_AM = t2f( dem_uncor_out_AM,NFFT); 

%% 画图，至少需要呈现这些图，并能与DSB调制方式、SSB_up调制方式进行比较
figure;
subplot(6,2,1)
plot(t,s_m); grid on;
xlabel('t(s)')
title('调制信号')
% xlim([0, 5]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
subplot(6,2,3)
plot(t,s_c); grid on;
xlabel('t(s)')
title('载波信号')
% xlim([0, 5]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
subplot(6,2,5)
plot(t,s_AM);hold on;
plot(t,s_AM_envelope,'r');grid on;
xlabel('t(s)')
title('AM信号')
% xlim([0, 5]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
subplot(6,2,7)
plot(t,s_DSB);hold on;
plot(t,s_DSB_envelope,'r');grid on;
xlabel('t(s)')
title('DSB信号')
% xlim([0, 5]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
subplot(6,2,9)
plot(t,s_SSB_up);grid on;
xlabel('t(s)')
title('SSB上边带信号')
% xlim([0, 5]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
subplot(6,2,11)
plot(t,s_SSB_up);grid on;
xlabel('t(s)')
title('SSB下边带信号')
% xlim([0, 5]); %这个不平滑
xlim([0, 0.01]); %这个没毛病

subplot(6,2,2)
plot(f,abs(S_m));grid on;
xlabel('f(Hz)')
title('调制信号的幅度-频率特性')
% xlim([-4 4]);
xlim([0 fs/2]);
subplot(6,2,4)
plot(f,abs(S_c));grid on;
xlabel('f(Hz)')
title('载波信号的幅度-频率特性')
% xlim([-15 15]);
xlim([0 fs/2]);
subplot(6,2,6)
plot(f,abs(S_AM));grid on;
xlabel('f(Hz)')
title('AM信号的幅度-频率特性')
% xlim([-15 15]);
xlim([0 fs/2]);
subplot(6,2,8)
plot(f,abs(S_DSB));grid on;
xlabel('f(Hz)')
title('DSB信号的幅度-频率特性')
% xlim([-15 15]);
xlim([0 fs/2]);
subplot(6,2,10)
plot(f,abs(S_SSB_up));grid on;
xlabel('f(Hz)')
title('SSB上边带信号的幅度-频率特性')
% xlim([-15 15]);
xlim([0 fs/2]);
subplot(6,2,12)
plot(f,abs(S_SSB_down));grid on;
xlabel('f(Hz)')
title('SSB下边带信号的幅度-频率特性')
% xlim([-15 15]);
xlim([0 fs/2]);


figure;
subplot(3,2,1)
plot(t,s_m); grid on;
xlabel('t(s)')
title('调制信号')
% xlim([0, 5]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
subplot(3,2,2)
plot(f,abs(2*S_m/fs));grid on;
xlabel('f(Hz)')
title('调制信号的幅度-频率特性')
% xlim([0 fs/2]);
subplot(3,2,3)
plot(t,s_c); grid on;
xlabel('t(s)')
title('载波信号')
% xlim([0, 5]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
subplot(3,2,4)
plot(f,abs(2*S_c/fs));grid on;
xlabel('f(Hz)')
title('载波信号的幅度-频率特性')
subplot(3,2,5)
plot(t,s_AM);hold on;
plot(t,s_AM_envelope,'r');grid on;
xlabel('t(s)')
title('AM信号')
% xlim([0, 5]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
subplot(3,2,6)
plot(f,abs(2*S_AM/fs));grid on;
xlabel('f(Hz)')
title('AM信号的幅度-频率特性')
% xlim([-15 15]);
% xlim([0 fs/2]);

figure
subplot(3,2,1)
plot(t,s_DSB);hold on;
plot(t,s_DSB_envelope,'r');grid on;
xlabel('t(s)')
title('DSB信号')
% xlim([0, 5]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
subplot(3,2,3)
plot(t,s_SSB_up);grid on;
xlabel('t(s)')
title('SSB上边带信号')
% xlim([0, 5]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
subplot(3,2,5)
plot(t,s_SSB_down);grid on;
xlabel('t(s)')
title('SSB下边带信号')
% xlim([0, 5]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
subplot(3,2,2)
plot(f,abs(2*S_DSB/fs));grid on;
xlabel('f(Hz)')
title('DSB信号的幅度-频率特性')
% xlim([-15 15]);
% xlim([0 fs/2]);
subplot(3,2,4)
plot(f,abs(2*S_SSB_up/fs));grid on;
xlabel('f(Hz)')
title('SSB上边带信号的幅度-频率特性')
% xlim([-15 15]);
% xlim([0 fs/2]);
subplot(3,2,6)
plot(f,abs(2*S_SSB_down/fs));grid on;
xlabel('f(Hz)')
title('SSB下边带信号的幅度-频率特性')
% xlim([-15 15]);
% xlim([0 fs/2]);



%相干解调信号过程
figure;
subplot(8,2,1);
plot(t, s_cor_temp_AM);grid on;
xlabel('t(s)');
% xlim([0, 10]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
title('AM相干解调通过乘法器后时域');
subplot(8,2,3);
plot(t, dem_cor_out_AM);grid on;
xlabel('t(s)');
% xlim([0, 10]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
title('AM相干解调输出');
subplot(8,2,5);
plot(t, s_cor_temp_DSB);grid on;
xlabel('t(s)');
% xlim([0, 10]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
title('DSB相干解调通过乘法器后时域');
subplot(8,2,7);
plot(t, dem_cor_out_DSB);grid on;
xlabel('t(s)');
% xlim([0, 10]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
title('DSB相干解调输出');
subplot(8,2,9);
plot(t, s_cor_temp_SSB_up);grid on;
xlabel('t(s)');
% xlim([0, 10]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
title('SSB上边带相干解调通过乘法器后时域');
subplot(8,2,11);
plot(t, dem_cor_out_SSB_up);grid on;
xlabel('t(s)');
% xlim([0, 10]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
title('SSB上边带相干解调输出');
subplot(8,2,13);
plot(t, s_cor_temp_SSB_down);grid on;
xlabel('t(s)');
% xlim([0, 10]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
title('SSB下边带相干解调通过乘法器后时域');
subplot(8,2,15);
plot(t, dem_cor_out_SSB_down);grid on;
xlabel('t(s)');
% xlim([0, 10]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
title('SSB下边带相干解调输出');


subplot(8,2,2);
plot(f, abs(S_cor_temp_AM));grid on;
xlabel('f(HZ)');
% xlim([-25, 25]);
xlim([0 fs/2]);
title('AM相干解调通过乘法器后频域');
subplot(8,2,4);
plot(f, abs(S_cor_out_AM));grid on;
xlabel('f(HZ)');
% xlim([-10, 10]);
xlim([0 fs/2]);
title('AM相干解调输出频域');
subplot(8,2,6);
plot(f, abs(S_cor_temp_DSB));grid on;
xlabel('f(HZ)');
% xlim([-25, 25]);
xlim([0 fs/2]);
title('DSB相干解调通过乘法器后频域');
subplot(8,2,8);
plot(f, abs(S_cor_out_DSB));grid on;
xlabel('f(HZ)');
% xlim([-10, 10]);
xlim([0 fs/2]);
title('DSB相干解调输出频域');
subplot(8,2,10);
plot(f, abs(S_cor_temp_SSB_up));grid on;
xlabel('f(HZ)');
% xlim([-25, 25]);
xlim([0 fs/2]);
title('SSB上边带相干解调通过乘法器后频域');
subplot(8,2,12);
plot(f, abs(S_cor_out_SSB_up));grid on;
xlabel('f(HZ)');
% xlim([-10, 10]);
xlim([0 fs/2]);
title('SSB上边带相干解调输出频域');
subplot(8,2,14);
plot(f, abs(S_cor_temp_SSB_down));grid on;
xlabel('f(HZ)');
% xlim([-25, 25]);
xlim([0 fs/2]);
title('SSB下边带相干解调通过乘法器后频域');
subplot(8,2,16);
plot(f, abs(S_cor_out_SSB_down));grid on;
xlabel('f(HZ)');
% xlim([-10, 10]);
xlim([0 fs/2]);
title('SSB下边带相干解调输出频域');

figure;
subplot(2,2,1);
plot(t, s_cor_temp_DSB);grid on;
xlabel('t(s)');
% xlim([0, 10]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
title('DSB相干解调通过乘法器后时域');
subplot(2,2,3);
plot(t, dem_cor_out_DSB);grid on;
xlabel('t(s)');
% xlim([0, 10]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
title('DSB相干解调输出');
subplot(2,2,2);
plot(f, abs(2*S_cor_temp_AM/fs));grid on;
xlabel('f(HZ)');
% xlim([-25, 25]);
title('AM相干解调通过乘法器后频域');
subplot(2,2,4);
plot(f, abs(2*S_cor_out_AM/fs));grid on;
xlabel('f(HZ)');
% xlim([-10, 10]);
title('AM相干解调输出频域');


figure;
subplot(2,2,1);
plot(t, s_cor_temp_DSB);grid on;
xlabel('t(s)');
% xlim([0, 10]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
title('DSB相干解调通过乘法器后时域');
subplot(2,2,3);
plot(t, dem_cor_out_DSB);grid on;
xlabel('t(s)');
% xlim([0, 10]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
title('DSB相干解调输出');
subplot(2,2,2);
plot(f, abs(2*S_cor_temp_DSB/fs));grid on;
xlabel('f(HZ)');
% xlim([-25, 25]);
% xlim([0 fs/2]);
title('DSB相干解调通过乘法器后频域');
subplot(2,2,4);
plot(f, abs(2*S_cor_out_DSB/fs));grid on;
xlabel('f(HZ)');
% xlim([-10, 10]);
% xlim([0 fs/2]);
title('DSB相干解调输出频域');

figure;
subplot(4,2,1);
plot(t, s_cor_temp_SSB_up);grid on;
xlabel('t(s)');
% xlim([0, 10]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
title('SSB上边带相干解调通过乘法器后时域');
subplot(4,2,3);
plot(t, dem_cor_out_SSB_up);grid on;
xlabel('t(s)');
% xlim([0, 10]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
title('SSB上边带相干解调输出');
subplot(4,2,5);
plot(t, s_cor_temp_SSB_down);grid on;
xlabel('t(s)');
% xlim([0, 10]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
title('SSB下边带相干解调通过乘法器后时域');
subplot(4,2,7);
plot(t, dem_cor_out_SSB_down);grid on;
xlabel('t(s)');
% xlim([0, 10]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
title('SSB下边带相干解调输出');
subplot(4,2,2);
plot(f, abs(2*S_cor_temp_SSB_up/fs));grid on;
xlabel('f(HZ)');
% xlim([-25, 25]);
% xlim([0 fs/2]);
title('SSB上边带相干解调通过乘法器后频域');
subplot(4,2,4);
plot(f, abs(2*S_cor_out_SSB_up/fs));grid on;
xlabel('f(HZ)');
% xlim([-10, 10]);
% xlim([0 fs/2]);
title('SSB上边带相干解调输出频域');
subplot(4,2,6);
plot(f, abs(2*S_cor_temp_SSB_down/fs));grid on;
xlabel('f(HZ)');
% xlim([-25, 25]);
% xlim([0 fs/2]);
title('SSB下边带相干解调通过乘法器后频域');
subplot(4,2,8);
plot(f, abs(2*S_cor_out_SSB_down/fs));grid on;
xlabel('f(HZ)');
% xlim([-10, 10]);
% xlim([0 fs/2]);
title('SSB下边带相干解调输出频域');




%非相干解调信号过程
figure;
subplot(2,1,1);
plot(t, dem_uncor_out_AM);grid on;
xlabel('t(s)');
% xlim([0, 10]); %这个不平滑
xlim([0, 0.01]); %这个没毛病
title('非相干解调输出');

subplot(2,1,2);
plot(f, abs(2*S_uncor_out_AM/fs));grid on;
xlabel('f(HZ)');
% xlim([-10 10]);
% xlim([0 fs/2]);
title('非相干解调输出频域');
