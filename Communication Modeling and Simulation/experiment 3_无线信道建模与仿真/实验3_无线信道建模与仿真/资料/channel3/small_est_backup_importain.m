clear
close all
clc

bandNum=1;
userNum=16;
sampleNum=128;
subcarrierNum1=1024*bandNum;
subcarrierNum2=1024*2*bandNum;
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
subchannel=32;
group_number=subcarrierNum1/subchannel;
frequency_channel_group=zeros(userNum,interval_num,subchannel);
for n=1:interval_num
	for i=1:userNum
        user=fading_factor_tmp2(i,:,n);
        frequency_channel_temp=fft(user,subcarrierNum1);
        frequency_channel(i,n,:)=10*log10(abs(frequency_channel_temp).^2);
        for ii=1:subchannel
            for iii=1:group_number 
                frequency_channel_group(i,n,ii)=frequency_channel_group(i,n,ii)+abs(frequency_channel_temp(1,iii+(ii-1)*group_number)).^2;
            end
        end        
    end
end

frequency_channel_group = 10*log10((frequency_channel_group/group_number));

figure(1)
for i=1:userNum
    for j=1:subcarrierNum1
        frequency_channel1(i,j)=frequency_channel(i,1,j);
    end
end
mesh(1:subcarrierNum1,1:userNum,frequency_channel1)

figure(2)
for i=1:userNum
    for j=1:subchannel
        frequency_channel3(i,j)=frequency_channel_group(i,1,j);
    end
end
mesh(1:subchannel,1:userNum,frequency_channel3)

figure(3)
for i=1:userNum
    user=fading_factor_tmp2(i,:,1);
    frequency_channel_compare(i,:)=10*log10(abs(fft(user,subcarrierNum2)).^2);
end
mesh(1:subcarrierNum2,1:userNum,frequency_channel_compare)

figure(4)
plot(1:1:subcarrierNum2, frequency_channel_compare(1,:));
hold on;
plot(1:(subcarrierNum2/subcarrierNum1):subcarrierNum2, frequency_channel1(1,:), '-r');
hold on;
plot(1:(subcarrierNum2/subchannel):subcarrierNum2, frequency_channel3(1,:), '-g');

figure(5)
for i=1:interval_num
    for j=1:subcarrierNum1
        frequency_channel2(i,j)=frequency_channel(1,i,j);
    end
end
mesh(1:subcarrierNum1,1:interval_num,frequency_channel2)

figure(6)
plot([1:sampleNum]*Ts2, 20*log10(abs(fading_factor_tmp)))
title('Rayleigh Fading')
xlabel('Time (sec)')
ylabel('Gain(dB)')
title('Fading Envelope (1024 samples)')

figure(7);
a_3D=abs(fading_factor_tmp1);
t=Ts:Ts:Ts*sampleNum;
tt=t(1:sampleNum);
pp=1:path_number;
y1=zeros(size(tt));
y2=ones(size(tt));
y3=2*ones(size(tt));
y4=3*ones(size(tt));
y5=4*ones(size(tt));
y6=5*ones(size(tt));
plot3(tt,y1 ,10*log10(a_3D(1,1:sampleNum)),tt,y2 ,20*log10(a_3D(2,1:sampleNum)),tt,y3 ,10*log10(a_3D(3,1:sampleNum)),tt,y4 ,10*log10(a_3D(4,1:sampleNum)),tt,y5 ,10*log10(a_3D(5,1:sampleNum)),tt,y6 ,10*log10(a_3D(6,1:sampleNum)))
grid on

figure(8);
a  = abs(fading_factor_tmp);
xi = real(fading_factor_tmp);
xq = imag(fading_factor_tmp);
subplot(2,2,1)
hist(a,100);
title('Envelope');

subplot(2,2,2)
hist(xi,100);
title('In-phase');

subplot(2,2,3)
hist(xq,100);
title('Quadrature');

figure(9);
correlation_flag=128;
subplot(3,1,1)
[Crosscorrelation1,lags] = xcov(fading_factor_tmp1(1,:),fading_factor_tmp1(2,:),correlation_flag,'coeff');
plot(lags,Crosscorrelation1)
grid
title('Cross-correlation of path1 & path2');

subplot(3,1,2)
[Crosscorrelation2,lags] = xcov(fading_factor_tmp1(1,:),fading_factor_tmp1(3,:),correlation_flag,'coeff');
plot(lags,Crosscorrelation2)
grid
title('Cross-correlation of path1 & path3');

subplot(3,1,3)
[Crosscorrelation3,lags] = xcov(fading_factor_tmp1(1,:),fading_factor_tmp1(4,:),correlation_flag,'coeff');
plot(lags,Crosscorrelation3)
grid
title('Cross-correlation of path1 & path4');
