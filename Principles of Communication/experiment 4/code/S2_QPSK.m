close all;
clear;

% --------QPSK Constellations
% clear all;
% close all;

x_qpsk = randi([0,1],1,10000)*2-1 + 1i*(randi([0,1],1,10000)*2-1); % case 1

scatterplot(x_qpsk);
xlabel('I');
ylabel('Q');
title('无噪声星座');
line([0 0],[-3 3],'Color','r');
line([-3 3],[0 0],'Color','r');
axis([-2 2 -2 2]);
grid on;



%% 请分别实现信噪比为30dB、10dB、0dB星座图
% A = 1;%振幅等于1
SNR_dB = 0:10:30; %信噪比——SNR dB值
SNR = 10.^(SNR_dB./10); %信噪比——SNR
Es = 1; % QPSK频带信号每个维度每符号能量
N0 = Es./SNR ; %QPSK频带信号每个维度对应的噪声功率
SNR_str = ['0' '10' '20' '30'];

for i = 1:length(SNR_dB)
    noise = sqrt(N0(i)) * (randn(1,length(x_qpsk)) + 1j*randn(1,length(x_qpsk))); %产生QPSK频带信号对应的加性高斯白噪声
    out_qpsk = noise + x_qpsk;
    scatterplot(out_qpsk);
    xlabel('I');
    ylabel('Q');
    if i == 1
        title([SNR_str(1,i) 'dB噪声星座']);
    else
        title([SNR_str(1,2*(i -1):2*(i -1) + 1) 'dB噪声星座']);
    end
    line([0 0],[-3 3],'Color','r');
    line([-3 3],[0 0],'Color','r');
    axis([-2 2 -2 2]);
    grid on;
end





