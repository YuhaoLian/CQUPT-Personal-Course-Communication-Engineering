clc
clear;
close all;

d0=1;%天线远场参考距离1m
gamma=3.71;%路径损耗因子
c=3*10^8;%光速
f=900*10^6;%频率900MHz
P0=43;%单位dBm
x=1:0.1:8;
d=2.^x;%单位千米，2Km到256Km
sigma=0.3;%阴影分布标准差

for i=1:length(d)
    L(i)=32.4+20*log10(f/1e6)+20*log10(d(i));
    shadow(i)=10*log10(lognrnd(0,sigma));
    L2(i)=L(i)+shadow(i);
    P1(i)=P0-L(i);
    P2(i)=P0-L2(i);
end;

figure(1)
plot(log10(d),P1,'r','linewidth',2)
hold on;
grid on;
plot(log10(d),P2,'b','linewidth',2)
xlabel('log10(距离(Km))')
ylabel('接收功率(dBm)')
legend('只考虑自由空间损耗的接收功率','考虑自由空间损耗+阴影的接收功率')

figure(2)
plot(d,P1,'r','linewidth',2)
hold on;
grid on;
plot(d,P2,'b','linewidth',2)
xlabel('距离(Km)')
ylabel('接收功率(dBm)')
legend('只考虑自由空间损耗的接收功率','考虑自由空间损耗+阴影的接收功率')

figure(3)
plot(d,L,'r','linewidth',2)
hold on;
grid on;
plot(d,L2,'b','linewidth',2)
xlabel('距离(Km)')
ylabel('损耗(dB)')
legend('自由空间损耗','自由空间损耗+阴影','Location','southeast')
%小区用户均匀分布
%范围x[-100,100],y[-100,100]
X=100;
Y=100;
x=[1:1:X];
y=[1:1:Y];

for i=1:X
for j=1:Y
      shadow1=10*log10(lognrnd(0,sigma));
      d1=sqrt(abs(x(i)-X/2)^2+abs(y(j)-Y/2)^2)/1000;
      L1(i,j)=32.45+20*log10(f/1e6)+20*log10(d1);
      P1(i,j)=P0-L1(i,j);
      P2(i,j)=shadow1;
      P3(i,j)=P1(i,j)+P2(i,j);
end
end;

figure(4)
mesh(x,y,P1)
title('自由空间损耗3D图');
xlabel('x(Km)')
ylabel('y(Km)')
zlabel('接收功率（dBm）')

figure(5)
mesh(x,y,P2)
title('阴影衰落3D图');
xlabel('x(Km)')
ylabel('y(Km)')
zlabel('接收功率（dBm）')

figure(6)
mesh(x,y,P3)
title('自由空间损耗+阴影3D图');
xlabel('x(Km)')
ylabel('y(Km)')
zlabel('接收功率（dBm）')
