clear;
close all;
%T.S.Rappaport�� ������ͨ��ԭ����Ӧ�ã��ڶ��棩�� ���ӹ�ҵ������ P53 ͼ5.23
% ������ʼ��
N = 1024*1024; % Ƶ��������
v = 33.33; % �ƶ��ն˵��ٶȣ���/�룩
c = 3e8; % ���٣���/�룩
fc = 900e6; % �ز�Ƶ��
lamda = c / fc; % �ز�����
fm = v / lamda; % ��������Ƶ��
deltaf = 2 * fm / (N-1); % Ƶ��������
T = 1 / deltaf; % ʱ��������
g = linspace(-fm, fm, N); % Ƶ�ʵ��ֵ

% �����ж�����Ƶ������Ķ����չ������ܶȺ��� sqrt(SEZ(f))
S = zeros(N, 1);
for i = 2:N-1
f = fc + g(i);
S(i) =1.5 / (pi * fm * (1 - ((f - fc) / fm) ^2) ^0.5);
end
figure(1)
plot(S)
S = sqrt(S);

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

figure(2)
hist(rt(1,10000:1:N/2),100);
grid on;
title('�������ͳ����');

figure(3)
[f,xi]=ksdensity(rt(1,10000:1:N/2));
plot(xi,f)
grid on;
title('����ĸ����ܶȺ���');

figure(4)
plot(0:1:256,10*log10(rt(1,10000:1:(10000+256))))
axis tight

grid on;
xlabel('������');
ylabel('��������');
title('ĳһ��ʱ���ڵİ������� (256��������)')

