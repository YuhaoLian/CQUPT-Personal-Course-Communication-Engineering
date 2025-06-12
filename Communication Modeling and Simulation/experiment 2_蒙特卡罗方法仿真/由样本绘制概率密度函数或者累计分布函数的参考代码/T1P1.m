close all 
clear all
clc
Randnum=(-2)*log(unifrnd(0,1,1,10000000));
[y,x]=ksdensity(Randnum);
plot(x,y,'bo')
hold on;
ezplot('0.5*exp(-0.5*x)',[0,20])