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
xlabel('距离(Km)')
grid on
ylabel('路径损耗(dB)')
title('COST231-WI模型路径损耗')
legend('自由空间f=1800MHz','视距路径损耗f=1800MHz','非视距路径损耗f=1800MHz','Location','southeast')