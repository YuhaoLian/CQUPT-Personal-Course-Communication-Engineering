clear;
close all;
%T.S.Rappaport著 《无线通信原理与应用（第二版）》 电子工业出版社 P53 图5.23、5.24

% 参数初始化
N = 256; % 频域样点数
v = 33.33; % 移动终端的速度（米/秒）
c = 3e8; % 光速（米/秒）
fc = 900e6; % 载波频率
lamda = c / fc; % 载波波长
fm = v / lamda; % 最大多普勒频移
deltaf = 2 * fm / (N-1); % 频域采样间隔
T = 1 / deltaf; % 时域采样间隔
g = linspace(-fm, fm, N); % 频率点的值
Path_delay = [0 5 7 12 17 25];%相对时延
Path_Gain = exp(-Path_delay);%每径平均功率，按照指数规递
PathNumber=length(Path_delay);%多径数

% 产生有多普勒频移引起的多普勒功率普密度函数 sqrt(SEZ(f))
S = zeros(N, 1);
for i = 2:N-1
f = fc + g(i);
S(i) = 1.5 / (pi * fm * (1 - ((f - fc) / fm) ^2) ^0.5);
end
S = sqrt(S);

fadingcoeff=zeros(PathNumber,N/2);

%合成径数
for ii=1:PathNumber
    
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
    rt=rt(1,1+N/2:end);
    %增加时延
    if ii==1 
        fadingcoeff(ii,:) = rt;
    else
        fadingcoeff(ii,:) = [ zeros(1,Path_delay(1,ii)) rt(1,[1:end-Path_delay(1,ii)])];
    end

end

t=T:T:T*N/2;
figure(1);
a_3D=abs(fadingcoeff);
tt=t(1:N/2);
pp=1:PathNumber;
y1=zeros(size(tt));
y2=ones(size(tt));
y3=2*ones(size(tt));
y4=3*ones(size(tt));
y5=4*ones(size(tt));
y6=5*ones(size(tt));
plot3(tt,y1 ,10*log10(a_3D(1,1:N/2)),tt,y2 ,10*log10(a_3D(2,1:N/2)),tt,y3 ,10*log10(a_3D(3,1:N/2)),tt,y4 ,10*log10(a_3D(4,1:N/2)),tt,y5 ,10*log10(a_3D(5,1:N/2)),tt,y6 ,10*log10(a_3D(6,1:N/2)))
grid on
title('不考虑平均功率递减时各径包络的比较');

correlation_flag=N/2;

figure(2);
subplot(3,1,1)
[Crosscorrelation1,lags] = xcov(fadingcoeff(1,:),fadingcoeff(2,:),correlation_flag,'coeff');
plot(lags,Crosscorrelation1)
grid
title('第1径和第3径的相关系数');

subplot(3,1,2)
[Crosscorrelation2,lags] = xcov(fadingcoeff(1,:),fadingcoeff(3,:),correlation_flag,'coeff');
plot(lags,Crosscorrelation2)
grid
title('第1径和第2径的相关系数');

subplot(3,1,3)
[Crosscorrelation3,lags] = xcov(fadingcoeff(1,:),fadingcoeff(4,:),correlation_flag,'coeff');
plot(lags,Crosscorrelation3)
grid
title('第1径和第4径的相关系数');

%考虑平均功率递减，各径叠后的最终结果
OutSignal=Path_Gain*fadingcoeff;
