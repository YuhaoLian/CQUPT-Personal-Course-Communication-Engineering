close all;
clear all;
clc;
%3000000个标准正态分布的样本
x=randn(3000000,1);
%拟合出300个样本的概率密度
[f,xi]=ksdensity(x);
%绘制图形
subplot(211)
plot(x)
title('样本数据(Sample Data)')
subplot(212)
plot(xi,f)
title('拟合出的概率密度函数(PDF)')