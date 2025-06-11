% clear; % 清除变量
clf; % 清除当前图像窗口
q=1.602e-19; % 元电荷电量
k=1/(4*pi*8.85e-12); % 静电力常量
a=1.6; % 正电荷的单位距离
b=1.6; % 负电荷的单位距离
x=-5:1:5; % 横坐标范围
y=x; % y坐标范围
z=x;
xx=-5:0.25:5; % 横坐标范围
yy=xx
zz=xx
points1 = zeros(41,41,41);
points2 = zeros(3,41*41*41);
% c = zeros(1,41*41*41);
A=[a,0,0]';
B=[-b,0,0];
for i=1:41
    for j=1:41
        for k=1:41
            points2(:,k+41*(j-1)+41*41*(i-1))=[xx(i),yy(j),zz(k)]';
        end
    end
end
%z=x.*y'  % z坐标范围
[X,Y,Z]=meshgrid(xx,yy,zz); % 设置坐标网格点 meshgrid生成网格采样点
r=sqrt((X-a).^2+(Y).^2+(Z).^2); % 第一个正电荷到场点的距离
R=sqrt((X+b).^2+(Y).^2+(Z).^2); % 第一个负电荷到场点的距离
V=q*k*(1./r-1./R); % 计算电势
[Ex,Ey,Ez]=gradient(-V); % 计算场强
AE=sqrt(Ex.^2+Ey.^2+Ez.^2);
M=AE(:,:,1)
% for i=1:41
%     for j=1:41
%         for k=1:41
%             M=AE(:,:,k)
%             c(1,k+41*(j-1)+41*41*(i-1))=M(i,j)
%         end
%     end
% end
x=-5:1.5:5; % 横坐标范围
y=x; % y坐标范围
z=x;
[X,Y,Z]=meshgrid(x,y,z); % 设置坐标网格点 meshgrid生成网格采样点
r=sqrt((X-a).^2+(Y).^2+(Z).^2); % 第一个正电荷到场点的距离
R=sqrt((X+b).^2+(Y).^2+(Z).^2); % 第一个负电荷到场点的距离
V=q*k*(1./r-1./R); % 计算电势
[Ex,Ey,Ez]=gradient(-V); % 计算场强
AE=sqrt(Ex.^2+Ey.^2+Ez.^2);
scatter3(points2(1,:),points2(2,:),points2(3,:),2,c,'.')
colorbar
Ex=Ex./AE;
Ey=Ey./AE;
Ez=Ez./AE;
% 场强归一化，使箭头等长
hold on;
quiver3(X,Y,Z,Ex,Ey,Ez,0.7) % 第五输入宗量0.7使场强箭头长短适中
hold on;
title('三维电偶极子周围的电场'), % 显示标题
hold on % 保持图像
plot3(1.6,0,0,'yo',1.6,0,0,'g+') % 用绿线画正电荷位置
plot3(-1.6,0,0,'go',-1.6,0,0,'g-') % 用黄线画负电荷位置
hold on;
figure(1)
xlabel('x'); % 显示横坐标
ylabel('y'), % 显示纵坐标
zlabel('z'); % 显示横坐标
% cv=linspace(min(min(min(V))),max(max(max(V))),1000); % 产生100个电位值
% verts = stream3(X,Y,Z,V,V,V,X,Y,Z);
% l=streamline(verts);
% set(l,'LineWidth',2);
% set(l,'Color','y');
view(3)

                    
                    


                           
         
