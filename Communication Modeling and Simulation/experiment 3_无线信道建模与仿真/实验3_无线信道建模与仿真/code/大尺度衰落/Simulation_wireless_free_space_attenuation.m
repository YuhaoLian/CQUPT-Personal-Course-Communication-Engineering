clc;
clear all
close all
d=0.1:0.1:100;
y1=wireless_free_space_attenuation(d,900);
y2=wireless_free_space_attenuation(d,1800);
plot(d,y1,d,y2);
xlabel('����(Km)');
ylabel('��ģ�dB)');
title('���ɿռ����');
legend('���ɿռ�f=900MHz','���ɿռ�f=1800MHz','Location','southeast')