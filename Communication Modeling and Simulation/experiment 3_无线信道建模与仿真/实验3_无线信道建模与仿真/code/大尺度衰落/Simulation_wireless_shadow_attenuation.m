clc
clear;
close all;

d0=1;%����Զ���ο�����1m
gamma=3.71;%·���������
c=3*10^8;%����
f=900*10^6;%Ƶ��900MHz
P0=43;%��λdBm
x=1:0.1:8;
d=2.^x;%��λǧ�ף�2Km��256Km
sigma=0.3;%��Ӱ�ֲ���׼��

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
xlabel('log10(����(Km))')
ylabel('���չ���(dBm)')
legend('ֻ�������ɿռ���ĵĽ��չ���','�������ɿռ����+��Ӱ�Ľ��չ���')

figure(2)
plot(d,P1,'r','linewidth',2)
hold on;
grid on;
plot(d,P2,'b','linewidth',2)
xlabel('����(Km)')
ylabel('���չ���(dBm)')
legend('ֻ�������ɿռ���ĵĽ��չ���','�������ɿռ����+��Ӱ�Ľ��չ���')

figure(3)
plot(d,L,'r','linewidth',2)
hold on;
grid on;
plot(d,L2,'b','linewidth',2)
xlabel('����(Km)')
ylabel('���(dB)')
legend('���ɿռ����','���ɿռ����+��Ӱ','Location','southeast')
%С���û����ȷֲ�
%��Χx[-100,100],y[-100,100]
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
title('���ɿռ����3Dͼ');
xlabel('x(Km)')
ylabel('y(Km)')
zlabel('���չ��ʣ�dBm��')

figure(5)
mesh(x,y,P2)
title('��Ӱ˥��3Dͼ');
xlabel('x(Km)')
ylabel('y(Km)')
zlabel('���չ��ʣ�dBm��')

figure(6)
mesh(x,y,P3)
title('���ɿռ����+��Ӱ3Dͼ');
xlabel('x(Km)')
ylabel('y(Km)')
zlabel('���չ��ʣ�dBm��')
