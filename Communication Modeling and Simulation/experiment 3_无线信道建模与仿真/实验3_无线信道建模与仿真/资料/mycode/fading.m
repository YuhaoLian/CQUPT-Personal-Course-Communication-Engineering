% Jakes model of a mobile Rayleigh fading channel.
% From Jakes, Microwave Mobile Communications, IEEE Press.
% EE 252. San Jose State University. Spring 2003.
clc;
clear;
close all;
Fc = 900e+06;        % carrier frequency (Hz) --> ISM band
wc = 2*pi*Fc;        % carrier frequency (rad/sec)
V = 3000/3600;       % mobile speed (= 3 Km/hr) in m/sec
N0 = 8;              % The number of paths

lambda = 3e+08/Fc;   % wavelength at Fc
Fm = V/lambda;       % Doppler shift
Fs = 10*Fm;           % sampling frequency
M = 100000;            % number of samples

paths = 16;
Ts = 1/Fs;           % sampling period
t = [0:Ts:(M-1)*Ts]; % time array

fading = jakes(Fm,Fs,M,N0,paths);

a = abs(fading);
xi = real(fading);
xq = imag(fading);
[R_xx,lags] = xcorr(fading,64,'coeff');
[Sx,w] = periodogram(R_xx,[],'twosided',128,Fs);

figure
subplot(3,1,1)
hist(a,100);
title('Envelope');
subplot(3,1,2)
hist(xi,100);
title('In-phase');
subplot(3,1,3)
hist(xq,100);
title('Quadrature');

figure
plot(t(1:128), 10*log10(a(1:128)));
axis([0 128/Fs min(10*log10(a(1:128)))-1 max(10*log10(a(1:128)))+1]);
xlabel('Time (sec)');
ylabel('Gain');
title('Fading Envelope (128 samples)')

figure;
subplot(2,1,1)
plot(lags,R_xx)
title('Autocorrelation');
grid
subplot(2,1,2)
psdplot(fftshift(Sx),w-Fs/2,'Hz','linear','PSD');
axis([-Fm-5 Fm+5 0 max(Sx)*1.2]);
