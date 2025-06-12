%Compare and analysis of Jakes model of a mobile rayleigh fading channel.
% Julius CQUPT Spring 2004.
clc;
clear;
close all;

fprintf ('Program of Compare and Analysis of Jakes Model of a Mobile Rayleigh Fading Channel\n\n');
fprintf ('To select "Jakes, Microwave Mobile Communications, 1973 IEEE Press." ,press "1"\n\n');
fprintf ('To select "Improved Models for the Generation of Multiple Uncorrelated Rayleigh Fading Waveforms, 2002 IEEE Communications Letters." ,press "2"\n\n');
fprintf ('To select "Modified Jakes` Model for Simulating Multiple Uncorrelated Fading Waveforms Model1,2000 IEEE Communications Letters." ,press "3"\n\n');
fprintf ('To select "Modified Jakes` Model for Simulating Multiple Uncorrelated Fading Waveforms Model2,2000 IEEE Communications Letters." ,press "4"\n\n');
fprintf ('To select "The modified Jakes` model, for which the first and second statistics of each path waveforms agree well with the theoretical expectations." ,press "5"\n\n');

defaults = input('Please press your seleciong(1~5):');
fprintf ('Please waite........');

Fc =  2e9;                  % carrier frequency (Hz)
V=100;                      % mobile speed  3 Km/hr
C=3e8;                      %light speed m/s
Fm=(V*Fc*1000)/(C*3600)    % Doppler shift
Fs = Fm;%3e6;%10000*Fm%20e6;    % sampling frequency(system band)
Ts=1/Fs;
infrom_length=1024*100;     % number of samples
t=Ts:Ts:Ts*infrom_length;
path_num=6;

if defaults==1
   fadingcoeff_3D=Normal_Jakes_Multipath_Est(infrom_length,Fc,V,Fs,path_num);
elseif defaults==2
   fadingcoeff_3D=Improved_Jakes_Multipath_Est(infrom_length,Fc,V,Fs,path_num); 
elseif defaults==3
   fadingcoeff_3D=Modified_Jakes_Multipath_Est1(infrom_length,Fc,V,Fs,path_num);
elseif defaults==4
   fadingcoeff_3D=Modified_Jakes_Multipath_Est2(infrom_length,Fc,V,Fs,path_num); 
elseif defaults==5
   fadingcoeff_3D=Standard_Jakes_Multipath_Est1(infrom_length,Fc,V,Fs,path_num); 
else
   fprintf ('Press a wrong number!');
   exit;
end

correlation_flag=128;
fadingcoeff=fadingcoeff_3D(1,:);
a = abs(fadingcoeff);
xi = real(fadingcoeff);
xq = imag(fadingcoeff);
[R_xx,lags] = xcov(fadingcoeff,correlation_flag,'coeff');
[Sx,w] = periodogram(R_xx,[],'twosided',correlation_flag*2,Fs);

figure(1);
subplot(2,2,1)
hist(a,100);
title('Envelope');

subplot(2,2,2)
hist(xi,100);
title('In-phase');

subplot(2,2,3)
hist(xq,100);
title('Quadrature');

figure(2);
subplot(2,1,1)
plot(t(1:1024), 10*log10(a(1:1024)));
axis([0 1024/Fs min(10*log10(a(1:1024)))-1 max(10*log10(a(1:1024)))+1]);
xlabel('Time (sec)');
ylabel('Gain(dB)');
title('Fading Envelope (1024 samples)')

subplot(2,1,2)
psdplot(fftshift(Sx),w-Fs/2,'Hz','linear','PSD');
grid on
axis([-Fm-5 Fm+5 0 max(Sx)*1.2]);

figure(3);
subplot(2,1,1)
plot(lags,R_xx)
grid
title('Auto-correlation');

subplot(2,1,2)
[Crosscorrelation,lags] = xcov(xq,xi,correlation_flag,'coeff');
plot(lags,Crosscorrelation)
grid
title('Cross-correlation of In-phase & Quadrature');

figure(4);
a_3D=abs(fadingcoeff_3D);
tt=t(1:1024);
pp=1:path_num;
y1=zeros(size(tt));
y2=ones(size(tt));
y3=2*ones(size(tt));
y4=3*ones(size(tt));
y5=4*ones(size(tt));
y6=5*ones(size(tt));
plot3(tt,y1 ,10*log10(a_3D(1,1:1024)),tt,y2 ,10*log10(a_3D(2,1:1024)),tt,y3 ,10*log10(a_3D(3,1:1024)),tt,y4 ,10*log10(a_3D(4,1:1024)),tt,y5 ,10*log10(a_3D(5,1:1024)),tt,y6 ,10*log10(a_3D(6,1:1024)))
grid on

figure(5);
subplot(3,1,1)
[Crosscorrelation1,lags] = xcov(fadingcoeff_3D(1,:),fadingcoeff_3D(2,:),correlation_flag,'coeff');
plot(lags,Crosscorrelation1)
grid
title('Cross-correlation of path1 & path2');

S = corrcoef(fadingcoeff_3D(1,:),fadingcoeff_3D(2,:))

subplot(3,1,2)
[Crosscorrelation2,lags] = xcov(fadingcoeff_3D(1,:),fadingcoeff_3D(3,:),correlation_flag,'coeff');
plot(lags,Crosscorrelation2)
grid
title('Cross-correlation of path1 & path3');

S = corrcoef(fadingcoeff_3D(1,:),fadingcoeff_3D(3,:))

subplot(3,1,3)
[Crosscorrelation3,lags] = xcov(fadingcoeff_3D(1,:),fadingcoeff_3D(4,:),correlation_flag,'coeff');
plot(lags,Crosscorrelation3)
grid
title('Cross-correlation of path1 & path4');

S = corrcoef(fadingcoeff_3D(1,:),fadingcoeff_3D(4,:))

% subplot(5,1,4)
% [Crosscorrelation4,lags] = xcov(fadingcoeff_3D(1,:),fadingcoeff_3D(5,:),correlation_flag,'coeff');
% plot(lags,Crosscorrelation4)
% grid
% title('Cross-correlation of path1 & path5');
% 
% subplot(5,1,5)
% [Crosscorrelation5,lags] = xcov(fadingcoeff_3D(1,:),fadingcoeff_3D(6,:),correlation_flag,'coeff');
% plot(lags,Crosscorrelation5)
% grid
% title('Cross-correlation of path1 & path6');

% clc;