clear;
close all;
%T.S.Rappaport�� ������ͨ��ԭ����Ӧ�ã��ڶ��棩�� ���ӹ�ҵ������ P53 ͼ5.23��5.24

% ������ʼ��
N = 256; % Ƶ��������
v = 33.33; % �ƶ��ն˵��ٶȣ���/�룩
c = 3e8; % ���٣���/�룩
fc = 900e6; % �ز�Ƶ��
lamda = c / fc; % �ز�����
fm = v / lamda; % ��������Ƶ��
deltaf = 2 * fm / (N-1); % Ƶ��������
T = 1 / deltaf; % ʱ��������
g = linspace(-fm, fm, N); % Ƶ�ʵ��ֵ
Path_delay = [0 5 7 12 17 25];%���ʱ��
Path_Gain = exp(-Path_delay);%ÿ��ƽ�����ʣ�����ָ�����
PathNumber=length(Path_delay);%�ྶ��

% �����ж�����Ƶ������Ķ����չ������ܶȺ��� sqrt(SEZ(f))
S = zeros(N, 1);
for i = 2:N-1
f = fc + g(i);
S(i) = 1.5 / (pi * fm * (1 - ((f - fc) / fm) ^2) ^0.5);
end
S = sqrt(S);

fadingcoeff=zeros(PathNumber,N/2);

%�ϳɾ���
for ii=1:PathNumber
    
	% ����ͬ������Դ
	Ia = randn(N/2, 1);
	Ib = randn(N/2, 1);
	%Ig_positive��Ig_negative����0��Ϊ���ĶԳƵ�
	Ig_positive = Ia + Ib*j; % �õ�Ig����������
	Ig_negative = flipud(conj(Ig_positive)); % �õ�Ig�ĸ�������
	Ig = [sqrt(real(Ig_positive).^2 + imag(Ig_positive).^2); sqrt(real(Ig_negative).^2 + imag(Ig_negative).^2)];
	Ifad = ifft(Ig .* S); %������Ҷ�任
	Ifad = abs(Ifad).^2; % ��ģ��ƽ��
	
	% ������������Դ
	Qa = randn(N/2, 1);
	Qb = randn(N/2, 1);
	%Qg_positive��Qg_negative����0��Ϊ���ĶԳƵ�
	Qg_positive = Qa + Qb*j; % �õ�Qg����������
	Qg_negative = flipud(conj(Qg_positive));% �õ�Qg�ĸ�������
	Qg = [sqrt(real(Qg_positive).^2 + imag(Qg_positive).^2); sqrt(real(Qg_negative).^2 + imag(Qg_negative).^2)];
	Qfad = ifft(Qg .* S);%������Ҷ�任
	Qfad = abs(Qfad) .* exp(j*(angle(Qfad)-pi/2)); % ��λ��ת90��
	Qfad = abs(Qfad).^2;% ��ģ��ƽ��
	
	% �ϲ���·����Դ
	rt = sqrt(Ifad + Qfad);
	mean = cumsum(rt);
	mean = mean(N)/N;
	rt=rt/mean;
	rt=rt';
    rt=rt(1,1+N/2:end);
    %����ʱ��
    if ii==1 
        fadingcoeff(ii,:) = rt;
    else
        fadingcoeff(ii,:) = [ zeros(1,Path_delay(1,ii)) rt(1,[1:end-Path_delay(1,ii)])];
    end

end

t=T:T:T*N/2;
figure(1);
a_3D=abs(fadingcoeff);
tt=t(1:N/2);
pp=1:PathNumber;
y1=zeros(size(tt));
y2=ones(size(tt));
y3=2*ones(size(tt));
y4=3*ones(size(tt));
y5=4*ones(size(tt));
y6=5*ones(size(tt));
plot3(tt,y1 ,10*log10(a_3D(1,1:N/2)),tt,y2 ,10*log10(a_3D(2,1:N/2)),tt,y3 ,10*log10(a_3D(3,1:N/2)),tt,y4 ,10*log10(a_3D(4,1:N/2)),tt,y5 ,10*log10(a_3D(5,1:N/2)),tt,y6 ,10*log10(a_3D(6,1:N/2)))
grid on
title('������ƽ�����ʵݼ�ʱ��������ıȽ�');

correlation_flag=N/2;

figure(2);
subplot(3,1,1)
[Crosscorrelation1,lags] = xcov(fadingcoeff(1,:),fadingcoeff(2,:),correlation_flag,'coeff');
plot(lags,Crosscorrelation1)
grid
title('��1���͵�3�������ϵ��');

subplot(3,1,2)
[Crosscorrelation2,lags] = xcov(fadingcoeff(1,:),fadingcoeff(3,:),correlation_flag,'coeff');
plot(lags,Crosscorrelation2)
grid
title('��1���͵�2�������ϵ��');

subplot(3,1,3)
[Crosscorrelation3,lags] = xcov(fadingcoeff(1,:),fadingcoeff(4,:),correlation_flag,'coeff');
plot(lags,Crosscorrelation3)
grid
title('��1���͵�4�������ϵ��');

%����ƽ�����ʵݼ���������������ս��
OutSignal=Path_Gain*fadingcoeff;
