function  shadow_fading(Dradius, f, area_id,scale)

D=10;               %In urban area, the correlation at distance 10m is 0.3
correlation=0.3;
a=-log(correlation)/D;

N=5;
M=2*N^2;           %Number of frequency samples
Naverage=1;       %Number of channel realizations in calculating channel correlation
                   %Increase Naverage to get more accurate correlation of
                   %simulated channel

if (area_id==1) || (area_id==2)
    SE=5.2;
elseif (area_id==3) || (area_id==4)
    SE=6.6;
else
    SE=0;
end

taoL=0.65*((log10(f/1e6))^2)-1.3*log10(f/1e6)+SE;

Gcutoff=1/1000;
fcutoff=a/2/pi*sqrt(Gcutoff^(-2/3)-1);               %cutoff frequency: \Phi(fcutoff)=\Phi(0)*Gcutoff
IntegrateF=1-a/(2*pi)/sqrt(fcutoff^2+a^2/4/pi^2);    %engergy in cutoff frequency

temp=1/fcutoff*10;
x=-temp:temp/100:temp;
L=length(x);
y=x;
   
u=rand(1,M);
fr=a/(2*pi)*((u.^(-2)-1).^0.5);
theta=pi*(rand(1,M)-0.5);
fxsample=fr.*cos(theta);
fysample=fr.*sin(theta);

fxsample=2*pi*fxsample;
fysample=2*pi*fysample;

c=sqrt(2/M);
theta=rand(M)*2*pi;

counter=0;
for k=1:scale:Dradius
    %for j=1:Dradius
    s=zeros(1,Dradius/scale);
    %counter=counter+1;
    data_temp1=1:scale:Dradius;
    data_temp2=k*ones(1,Dradius/scale);
    k
    for i=1:M
        s=s+(cos(fxsample(i)*data_temp1+fysample(i)*data_temp2+theta(i)))*c*taoL;
    end
    if k==1
        fp_shadow = fopen('shadow_fading.bin', 'wt');  
        fprintf(fp_shadow, '%d\t', s);
        fclose(fp_shadow);
    else
        fp_shadow = fopen('shadow_fading.bin', 'at');  
        fprintf(fp_shadow, '%d\t', s);
        fclose(fp_shadow);   
    end
        %end
end


% x=1:1:Dradius;%/fcutoff*10;
% y=1:1:Dradius;
% L=length(x);
% X=[];
% Y=[];
% for i=1:L
%     X=[X; x];
%     Y=[Y, y'];
% end
% s=zeros(L,L);
% for i=1:M
%     s=s+cos(fxsample(i)*X+fysample(i)*Y+theta(i));
% end
% s=c*s*taoL;
% out_shadow=s;



