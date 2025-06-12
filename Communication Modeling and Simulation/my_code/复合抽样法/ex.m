clc;
clear;


figure;
xR = Composite_sampling(100);
subplot(2,2,1);
cdfplot(xR); 
hold on
ezplot('0.6*(1-exp(-0.5*x))+0.4*(1-exp(-3*x))',[0,10])
hold off
xlabel('x');
ylabel('f(x)');
title('100个随机数')

xR = Composite_sampling(1000);
subplot(2,2,2);
cdfplot(xR); 
hold on
ezplot('0.6*(1-exp(-0.5*x))+0.4*(1-exp(-3*x))',[0,10])
hold off
xlabel('x');
ylabel('f(x)');
title('1000个随机数')

xR = Composite_sampling(10000);
subplot(2,2,3);
cdfplot(xR); 
hold on
ezplot('0.6*(1-exp(-0.5*x))+0.4*(1-exp(-3*x))',[0,10])
hold off
xlabel('x');
ylabel('f(x)');
title('10000个随机数')

xR = Composite_sampling(100000);
subplot(2,2,4);
cdfplot(xR); 
hold on
ezplot('0.6*(1-exp(-0.5*x))+0.4*(1-exp(-3*x))',[0,10])
hold off
xlabel('x');
ylabel('f(x)');
title('100000个随机数')