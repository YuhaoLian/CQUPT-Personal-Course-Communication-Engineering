clc;
close all
clear all
f=900;
d=0.1:0.1:100;
y0=wireless_free_space_attenuation(d,f);
d=0.1:0.1:100;
y1=wireless_hata_attenuation(1,f,1.5,50,d);
y2=wireless_hata_attenuation(2,f,1.5,50,d);
y3=wireless_hata_attenuation(3,f,1.5,100,d);
y4=wireless_hata_attenuation(4,f,1.5,100,d);
plot(d,y0,d,y1,d,y2,d,y3,d,y4);
xlabel('����(Km)');
ylabel('��ģ�dB)');
grid on;
title('Okumura-Hataģ�����');
legend('���ɿռ�f=900MHz','��С����f=900MHz,Hm=1.5m,Hb=50m','�����f=900MHz,Hm=1.5m,Hb=50m','����f=900MHz,Hm=1.5m,Hb=100m','ũ��f=900MHz,Hm=1.5m,Hb=100m','Location','southeast','Location','southeast')