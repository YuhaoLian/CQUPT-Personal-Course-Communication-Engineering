clc;
clear all;
close all;
f=1800;
d=0.02:0.01:5;
y0=wireless_free_space_attenuation(d,f);
y1=wireless_Walfish_Ikegami_LOS_attenuation(f,d);
Model=1;
Hm=1.5;
Hb=17;
w=20;
b=40;
Phi=90;
Hroof=15;
y2=wireless_Walfish_Ikegami_NLOS_attenuation(Model,f,d,Hm,Hb,Hroof,w,b,Phi);
plot(d,y0,d,y1,d,y2);
xlabel('����(Km)')
grid on
ylabel('·�����(dB)')
title('COST231-WIģ��·�����')
legend('���ɿռ�f=1800MHz','�Ӿ�·�����f=1800MHz','���Ӿ�·�����f=1800MHz','Location','southeast')