function  [fading_factor]=small_term_fading(sampleNum, fs, Ts2)

Ts=1/fs;

[Path_Delay,Path_loss,Doppler_Frequency]=multipath_detail(fs);
path_number=length(Path_Delay);

MeanAmp = 10.^(Path_loss/10);
MeanAmp = sqrt( MeanAmp/sum(MeanAmp) );	%Nomorlized the power to 1.

N=200; %number of input waves
n=[0:N-1];
cita=2*pi*n/N; %input angles


    
for k=1:path_number
	alfa=randn(1,N); %magnitudes of input waves, path k
	alfanormed(k,:)=alfa/sqrt(sum(alfa.^2)); %normalized magnitude, path k
    phaseinit=randn(1,N);
	phaseinit(k,:)=phaseinit*2*pi/max(phaseinit); %initial phases, path k
end


for j=1:path_number
    
    TempSig=zeros(1,sampleNum);
    for i=0:sampleNum-1
        t=i*Ts2;%Ts*6e6;
        is=MeanAmp(1,j)*sum(alfanormed(j,:).*cos(2*pi*Doppler_Frequency(1,path_number)*t*cos(cita)+phaseinit(j,:)));
        qs=MeanAmp(1,j)*sum(alfanormed(j,:).*sin(2*pi*Doppler_Frequency(1,path_number)*t*cos(cita)+phaseinit(j,:)));
        I(j,i+1)=is;
        Q(j,i+1)=qs;
        fading_factors=qs+sqrt(-1)*is;
        TempSig(1,i+1)=fading_factors;
    end
    
    fading_factor(j,:) =  TempSig;
    
end

