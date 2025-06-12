clear all;
clc
N=16;
fT=0.1;
mm=0:1:(N-1);
m=mm-N/2;
c=(sin(pi*(m+fT)).*exp(i*pi*(N-1)*(m+fT)/N))./(N*sin(pi*(m+fT)/N));
plot(mm,real(c),'s',mm,imag(c),'d',mm,abs(c),'o');
legend('实部','虚部','模值')
xlabel('N');
ylabel('子载波偏差');
