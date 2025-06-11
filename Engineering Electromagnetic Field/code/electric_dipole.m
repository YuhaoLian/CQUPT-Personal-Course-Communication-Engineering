clear; % 清除变量
clf; % 清除当前图像窗口
q=1.602e-19; % 元电荷电量
k=1/(4*pi*8.85e-12); % 静电力常量
a=1.6; % 正电荷的单位距离
b=1.6; % 负电荷的单位距离
x=-10:0.5:10; % 横坐标范围
y=x; % 纵坐标范围
[X,Y]=meshgrid(x,y); % 设置坐标网格点 meshgrid生成网格采样点
rp=sqrt((X-a).^2+(Y).^2); % 第一个正电荷到场点的距离
rm=sqrt((X+b).^2+(Y).^2); % 第一个负电荷到场点的距离
V=q*k*(1./rp-1./rm); % 计算电势
[Ex,Ey]=gradient(-V); % 计算场强
AE=sqrt(Ex.^2+Ey.^2);
Ex=Ex./AE;
Ey=Ey./AE; % 场强归一化，使箭头等长
cv=linspace(min(min(V)),max(max(V)),250); % 产生250个电位值
contourf(X,Y,V,cv,'r-') % 用红实线画填色等位线图
title('电偶极子的场'), % 显示标题
hold on % 保持图像
quiver(X,Y,Ex,Ey,0.7) % 第五输入宗量0.7使场强箭头长短适中
plot(a,0,'wo',a,0,'w+') % 用绿线画正电荷位置
plot(-b,0,'go',-b,0,'g-') % 用黄线画负电荷位置
xlabel('x'); % 显示横坐标
ylabel('y'), % 显示纵坐标
hold off % 保持图像
surf(X,Y,V)
xlabel('\itx'); 
zlabel('\itz');
hold off % 保持图像
mesh(X,Y,AE)
xlabel('\itx'); 
zlabel('\itz');
hold off % 保持图像

