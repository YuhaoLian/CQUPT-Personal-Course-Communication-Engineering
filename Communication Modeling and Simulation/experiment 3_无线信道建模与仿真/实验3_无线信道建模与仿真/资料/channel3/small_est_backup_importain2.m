clear
close all
clc

bandNum=2;
userNum=16;
sampleNum=128;
subcarrierNum=1024*bandNum;
subchannelNum=32*bandNum;
fs=6e6*bandNum; % sampling frequency or system bandwidth
Ts=1/fs;
Ts2=2e-2;

[Path_Delay,Path_loss,Doppler_Frequency]=multipath_detail(fs);
path_number=length(Path_Delay);

[fading_factor]=small_term_fading(userNum, sampleNum, fs, Ts2);

delay_max=max(Path_Delay)+1;
fading_factor_tmp=zeros(1,sampleNum);
fading_factor_tmp1=zeros(path_number,sampleNum);
fading_factor_tmp2=zeros(userNum,delay_max,sampleNum);

for i=1:sampleNum
    for j=1:path_number
        fading_factor_tmp1(j,i)=fading_factor(1,j,i);
	end
end

k=0;
n=1;
for j=1:delay_max
    k=k+1;
    if k==Path_Delay(1,n)+1
        fading_factor_tmp2(:,j,:)=fading_factor(:,n,:);
        n=n+1;
    else
        fading_factor_tmp2(:,j,:)=0;
    end
end

fading_factor_tmp=sum(fading_factor_tmp1);

interval_num=sampleNum;
SubcarrierNumPerSubchannle=subcarrierNum/subchannelNum;
frequency_subchannel_time=zeros(userNum,interval_num,subchannelNum);
for n=1:interval_num
	for i=1:userNum
        user=fading_factor_tmp2(i,:,n);
        frequency_subcarrier_temp=fft(user,subcarrierNum);
        frequency_subcarrier_time(i,n,:)=10*log10(abs(frequency_subcarrier_temp).^2);
        for ii=1:subchannelNum
            for iii=1:SubcarrierNumPerSubchannle 
                frequency_subchannel_time(i,n,ii)=frequency_subchannel_time(i,n,ii)+abs(frequency_subcarrier_temp(1,iii+(ii-1)*SubcarrierNumPerSubchannle)).^2;
            end
        end        
    end
end

frequency_subchannel_time = 10*log10((frequency_subchannel_time/SubcarrierNumPerSubchannle));

figure(1)
for i=1:userNum
    for j=1:subcarrierNum
        frequency_subcarrier1(i,j)=frequency_subcarrier_time(i,1,j);
    end
end
mesh(1:subcarrierNum,1:userNum,frequency_subcarrier1)

figure(2)
for i=1:userNum
    for j=1:subchannelNum
        frequency_subchannel2(i,j)=frequency_subchannel_time(i,1,j);
    end
end
mesh(1:subchannelNum,1:userNum,frequency_subchannel2)

figure(3)
plot(1:1:subcarrierNum, frequency_subcarrier1(1,:), '-r');
hold on;
plot(1+(SubcarrierNumPerSubchannle):(SubcarrierNumPerSubchannle):subcarrierNum+(SubcarrierNumPerSubchannle), frequency_subchannel2(1,:), '-g');

figure(4)
for i=1:interval_num
    for j=1:subchannelNum
        frequency_subchannel1(i,j)=frequency_subchannel_time(1,i,j);
    end
end
mesh(1:subchannelNum,1:interval_num,frequency_subchannel1)

figure(5)
plot([1:sampleNum]*Ts2, 20*log10(abs(fading_factor_tmp)))
title('Rayleigh Fading')
xlabel('Time (sec)')
ylabel('Gain(dB)')
title('Fading Envelope')

save ('cellptest.bat','frequency_subchannel_time')

% figure(6);
% a_3D=abs(fading_factor_tmp1);
% t=Ts:Ts:Ts*sampleNum;
% tt=t(1:sampleNum);
% pp=1:path_number;
% y1=zeros(size(tt));
% y2=ones(size(tt));
% y3=2*ones(size(tt));
% y4=3*ones(size(tt));
% y5=4*ones(size(tt));
% y6=5*ones(size(tt));
% plot3(tt,y1 ,10*log10(a_3D(1,1:sampleNum)),tt,y2 ,20*log10(a_3D(2,1:sampleNum)),tt,y3 ,10*log10(a_3D(3,1:sampleNum)),tt,y4 ,10*log10(a_3D(4,1:sampleNum)),tt,y5 ,10*log10(a_3D(5,1:sampleNum)),tt,y6 ,10*log10(a_3D(6,1:sampleNum)))
% grid on
% 
% figure(7);
% a  = abs(fading_factor_tmp);
% xi = real(fading_factor_tmp);
% xq = imag(fading_factor_tmp);
% subplot(2,2,1)
% hist(a,100);
% title('Envelope');
% 
% subplot(2,2,2)
% hist(xi,100);
% title('In-phase');
% 
% subplot(2,2,3)
% hist(xq,100);
% title('Quadrature');
% 
% figure(8);
% correlation_flag=128;
% subplot(3,1,1)
% [Crosscorrelation1,lags] = xcov(fading_factor_tmp1(1,:),fading_factor_tmp1(2,:),correlation_flag,'coeff');
% plot(lags,Crosscorrelation1)
% grid
% title('Cross-correlation of path1 & path2');
% 
% subplot(3,1,2)
% [Crosscorrelation2,lags] = xcov(fading_factor_tmp1(1,:),fading_factor_tmp1(3,:),correlation_flag,'coeff');
% plot(lags,Crosscorrelation2)
% grid
% title('Cross-correlation of path1 & path3');
% 
% subplot(3,1,3)
% [Crosscorrelation3,lags] = xcov(fading_factor_tmp1(1,:),fading_factor_tmp1(4,:),correlation_flag,'coeff');
% plot(lags,Crosscorrelation3)
% grid
% title('Cross-correlation of path1 & path4');
