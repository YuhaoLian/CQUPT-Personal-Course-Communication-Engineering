clc;
clear all
close all
d=0.1:0.1:100;
y1=wireless_free_space_attenuation(d,900);
y2=wireless_free_space_attenuation(d,1800);
plot(d,y1,d,y2);
xlabel('距离(Km)');
ylabel('损耗（dB)');
title('自由空间损耗');
legend('自由空间f=900MHz','自由空间f=1800MHz','Location','southeast')