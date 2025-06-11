clear; % 清除变量
clf; % 清除当前图像窗口
I=1.602e-6; % 元电流大小
R=0.3;
k=1/1e-7; % 磁场常量
S=2*pi*R*R;
x=-10:0.5:10; % 横坐标范围
y=x; % 纵坐标范围
[X,Y]=meshgrid(x,y); % 设置坐标网格点 meshgrid生成网格采样点
r=sqrt((X).^2+(Y).^2); % 第一个顺磁偶极子到场点的距离
B=(I*k*S./((r).^3))*2*Y./r; % 计算磁感应强度
Bx=(I*k*S./((X).^3)).*3.*(X.*Y)./((r).^2);
By=(I*k*S./((X).^3)).*(2-3.*(X.*X)./((r).^2));
title('磁偶极子周围的磁场'), % 显示标题
hold on % 保持图像
figure(1);
l = streamslice(X,Y,Bx,By);
set(l,'LineWidth',2)
set(l,'Color','r');
xlabel('x'); % 显示横坐标
ylabel('y'), % 显示纵坐标
plot(R,0,'gx');
plot(-R,0,'go');
figure(2);
mesh(X,Y,B)

