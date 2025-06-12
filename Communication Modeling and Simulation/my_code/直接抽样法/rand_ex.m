
clear all;
clc;


%% ---------------- Parameters ---------------%%
Lambda = 1/3;

%% ----------------   指数分布  ---------------%%
figure;

% 100个随机数
Randnum=-log(unifrnd(0,1,1,100))/Lambda;
[y,x]=ksdensity(Randnum);
subplot(2,2,1);
plot(x,y,'o','MarkerFaceColor','r')
hold on;
ezplot('1/3*exp(-1/3*x)',[0,10])
xlabel('x');
ylabel('f(x)');
title('100个随机数')

% 1000个随机数
Randnum=-log(unifrnd(0,1,1,1000))/Lambda;
[y,x]=ksdensity(Randnum);
subplot(2,2,2);
plot(x,y,'o','MarkerFaceColor','r')
hold on;
ezplot('1/3*exp(-1/3*x)',[0,10])
xlabel('x');
ylabel('f(x)');
title('1000个随机数')

% 10000个随机数
Randnum=-log(unifrnd(0,1,1,10000))/Lambda;
[y,x]=ksdensity(Randnum);
subplot(2,2,3);
plot(x,y,'o','MarkerFaceColor','r')
hold on;
ezplot('1/3*exp(-1/3*x)',[0,10])
xlabel('x');
ylabel('f(x)');
title('10000个随机数')

% 100000个随机数
Randnum=-log(unifrnd(0,1,1,10000000))/Lambda;
[y,x]=ksdensity(Randnum);
subplot(2,2,4);
plot(x,y,'o','MarkerFaceColor','r')
hold on;
ezplot('1/3*exp(-1/3*x)',[0,10])
xlabel('x');
ylabel('f(x)');
title('10000000个随机数')