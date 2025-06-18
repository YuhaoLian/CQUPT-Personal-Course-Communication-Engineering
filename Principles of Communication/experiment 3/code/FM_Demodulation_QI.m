function [sft_demodulation] = FM_Demodulation_QI(sft,t,fc,fm,fs)
%% =============================================================
%%******************************    QI-Demodulation   *********************************
%%=============================================================
i_data=zeros(1,length(t));
q_data=zeros(1,length(t));
i_data = sft .* cos(2*pi*fc*t);
q_data = sft .* sin(2*pi*fc*t);

%低通滤波器除去I、Q的高频分量
wc=2*2*pi*fm/fs;                        %解调器低通滤波器截止频率，2fm
lbf_hn=fir1(32,wc/pi);                     %低通滤波器系统函数，FIR滤波器
[Hm,w] = freqz(lbf_hn, 1, 1024, fs);
adc_i = filter(lbf_hn,1,i_data);
adc_q = filter(lbf_hn,1,q_data);


%FM 解调(1)
c_len = length(adc_i);
cr = zeros(1,c_len);
cj = zeros(1,c_len);

for i=2:c_len
    cr(i) = adc_i(i)*adc_i(i) + adc_q(i)*adc_q(i); %I(n)^2+Q(n)^2
    cj(i) = adc_i(i-1)*adc_q(i) - adc_i(i)*adc_q(i-1);  %I(n-1)*Q(n) -I(n)*Q(n-1)
end

sft_demodulation= zeros(1,c_len);

for i=1:c_len
    if cr(i) == 0
        sft_demodulation(i) =0;
    else
        sft_demodulation(i) =(cj(i)/cr(i));
    end
end

end

