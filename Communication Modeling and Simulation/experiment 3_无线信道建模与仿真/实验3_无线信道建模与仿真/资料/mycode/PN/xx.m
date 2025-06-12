% function  [max_index,ti]=downlink_OFDM_blind(sub_number,SNR,frequency_space,delta,timing_offset)
%sub_number--the number of subcarrier
%SNR--Signal to Noise Ratio
%frequency_space-- the interval of the adjacent pilot bits
%delta--the frequency offset normalized
%timing_offset-- the timing_offset

clear;
clc;
sub_number=128;
SNR=1;
frequency_space=128;
delta=0.5;
timing_offset=0;

frame_number=1;
simulate_number=sub_number*frame_number;
u_number=1;
user_number=1;
%the pilot data
for  i=1:simulate_number
    if (mod((i-1),frequency_space)~=0)
     user(1,i)=0;
   else
     user(1,i)=1024;  
   end
end
%form the frame
for  i=1:frame_number
    for jj=1:u_number
        for k=1:sub_number
            trans_signal((i-1)*sub_number*u_number+(jj-1)*sub_number+k)=user(jj,(i-1)*sub_number+k);
        end
    end
end
temp_3=trans_signal(1:128);
trans_signal(129:256)=temp_3;
%trans_signal
%perform ifft+cyclic prefix
%cprefix_index=0;
%cprefix_index=0;
cprefix_index=sub_number/4;
for  i=1:frame_number*u_number
    for  jj=((i-1)*sub_number+1):i*sub_number
         temp(jj-(i-1)*sub_number)=trans_signal(jj);
    end
    temp1=ifft(temp);
    for  k=((i-1)*sub_number+1):i*sub_number
         trans_signal(k)=temp1(k-(i-1)*sub_number);                                 %trans_signal is the ifft of trans_signal3
    end
    for  t=1:cprefix_index
         temp2(t)=temp1(sub_number-cprefix_index+t);
    end
    for  m=(cprefix_index+1):(sub_number+cprefix_index)
         temp2(m)=temp1(m-cprefix_index);
     end
    for  l=((i-1)*(sub_number+cprefix_index)+1):(i*(sub_number+cprefix_index))
         trans_signal1(l)=temp2(l-(i-1)*(sub_number+cprefix_index));                %trans_signal1 is the trans_signal+cyclic prefix
     end
 end
 %In order to prevent the matrix to across the bound 
  trans_signal1((frame_number*u_number)*((sub_number+cprefix_index)+1):(frame_number*u_number+1)*(sub_number+cprefix_index))=0;
  %trans_signal1
  M=sub_number/frequency_space
  reference_signal(1:frame_number*sub_number)=0;
  reference_signal1(1:frame_number*sub_number)=0;
  reference_signal1=pn_user1(frequency_space,0);
  %reference_signal=reference_signal1;
  %reference_signal1=pn_user1(frequency_space,0);

for i=1:frequency_space
    reference_signal((i-1)*M+1)=reference_signal1(i);
end
% reference_signal
%for i=1:sub_number
%    trans_signal1(cprefix_index+i)=reference_signal(i)*trans_signal1(i+cprefix_index);
%    trans_signal1(sub_number+2*cprefix_index+i)=reference_signal(i)*trans_signal1(i+2*cprefix_index+sub_number);
%end

%trans_signal1
%In order to prevent the matrix to across the bound 
trans_signal1((frame_number)*((sub_number+cprefix_index)+1):(frame_number+1)*(sub_number+cprefix_index))=0;
M=sub_number/frequency_space;
reference_signal(1:frame_number*(sub_number+cprefix_index))=0;
reference_signal1(1:frame_number*(sub_number+cprefix_index))=0;
reference_signal1=pn_user1(frame_number*frequency_space*1.25,0);
%separate the channel signal
for i=1:frame_number*frequency_space
    reference_signal((i-1)*M+1)=reference_signal1(i);
end
%reference_signal
%trans_signal1
for i=1:frame_number*(sub_number+cprefix_index)
    trans_signal1(i)=reference_signal(i)*trans_signal1(i);
end
% trans_signal1

%Add the multipath channel
[trans_signal1]=Standard_Multipath1(trans_signal1,2e9,50,1/2e6);
%Add the timing error
trans_signal_t(1:(frame_number*u_number+1)*(sub_number+cprefix_index))=0;
for i=1:frame_number*user_number
      
     for l=((i-1)*(sub_number+cprefix_index)+1:(i*(sub_number+cprefix_index)))
         trans_signal_t(l+timing_offset)=trans_signal1(l);
     end
end
%Add the frequency offset
for  jj=1:(sub_number+cprefix_index)*(frame_number*user_number+1)
    trans_signal4(jj)=trans_signal_t(jj)*exp(j*2*pi*delta*(jj-1)/sub_number);
end

  
%Added White Noise Channel
trans_signal4=awgn(trans_signal4,SNR,'measured');
%New timing correction algorithm
sum(1:M)=0;
reference_signal
for k=1:M
    temp1(1:frame_number*(sub_number+cprefix_index))=0;
    for kk=1:frame_number*(sub_number+cprefix_index)
        temp1(kk+k-1)=reference_signal(kk);
    end
    temp1
    for i=1:frame_number*(sub_number+cprefix_index)
         sum(k)=sum(k)+conj(temp1(i))*trans_signal4(i);
    end
end

sum_data=abs(sum)/2;
[max_value,max_index]=max(abs(sum));
max_value
max_index=max_index-1
%max_index=0;
%Compansating the timing error according to the result of the above algorithm
for i=1:frame_number*user_number
     for l=((i-1)*(sub_number+cprefix_index)+1+max_index:(i*(sub_number+cprefix_index))+max_index)
         trans_signal5(l-max_index)=trans_signal4(l);
     end
end
%trans_signal5
%PN code removed
for i=1:frame_number*(sub_number+cprefix_index)
    trans_signal6(i)=conj(reference_signal(i))*trans_signal5(i)/2;
end
%trans_signal6
frac_frequency_offset=0;
for i=1:frame_number*frequency_space-1
    frac_frequency_offset=frac_frequency_offset+angle(conj(trans_signal6((i-1)*M+1))*trans_signal6(i*M+1))*sub_number/(2*pi*M);
end
ti=frac_frequency_offset/(frame_number*frequency_space-1)
%frequency offset estimation
%fpname='d:\dissertation-program\data_pn5db.mat'; 
%fp = fopen(fpname,'a+');
%for M1 =1:128
%    fprintf(fp,'%d\t%10.8f\n',M1,sum_data(M1));
%end
%fclose(fp);
