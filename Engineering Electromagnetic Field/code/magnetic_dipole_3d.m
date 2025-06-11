% clear; % 清除变量
clf; % 清除当前图像窗口
I=1.602e-8; % 元电流大小
R=0.3;
k=1/1e-7; % 磁场常量
S=2*pi*R*R;
x=-4:1:4; % 横坐标范围
y=x; % 纵坐标范围
z=y;
xx=-4:0.25:4; % 横坐标范围
yy=xx
zz=xx
points1 = zeros(33,33,33);
points2 = zeros(3,33*33*33);
% c = zeros(1,33*33*33);
for i=1:33
    for j=1:33
        for k=1:33
            points2(:,k+33*(j-1)+33*33*(i-1))=[xx(i),yy(j),zz(k)]';
        end
    end
end
[Xx,Yy,Zz]=meshgrid(xx,yy,zz); % 设置坐标网格点 meshgrid生成网格采样点
r=sqrt((Xx).^2+(Yy).^2); % 第一个顺磁偶极子到场点的距离
rr=sqrt((Xx).^2+(Zz).^2);
Bx=(I*k*S./((Xx).^3)).*3.*(Xx.*Zz)./((rr).^2).*(Xx./r);
By=Bx;
Bz=(I*k*S./((Xx).^3)).*(2-3.*(Xx.*Xx)./((rr).^2));
B=sqrt(Bx.*Bx+By.*By+Bz.*Bz)
% for i=1:33
%     for j=1:33
%         for k=1:33
%             M=B(:,:,k)
%             c(1,k+33*(j-1)+33*33*(i-1))=M(i,j)
%         end
%     end
% end
[X,Y,Z]=meshgrid(x,y,z); % 设置坐标网格点 meshgrid生成网格采样点
r=sqrt((X).^2+(Y).^2); % 第一个顺磁偶极子到场点的距离
rr=sqrt((X).^2+(Z).^2);
Bx=(I*k*S./((X).^3)).*3.*(X.*Z)./((rr).^2).*(X./r);
By=Bx;
Bz=(I*k*S./((X).^3)).*(2-3.*(X.*X)./((rr).^2));
title('磁偶极子周围的磁场'), % 显示标题
hold on % 保持图像
scatter3(points2(1,:),points2(2,:),points2(3,:),2,c,'.')
colorbar
xlabel('x'); % 显示横坐标
ylabel('y'), % 显示纵坐标
h = 0; % 高度
pos = [0,0]; % 圆心位置
t=0:0.001:(2*pi);  % 圆滑性设置
t=[t,0];
o=plot3(pos(1)+R*sin(t),pos(2)+R*cos(t), h*ones(size(t)));
set(o,'LineWidth',3);
hold on;
verts = stream3(X,Y,Z,Bx,By,Bz,X,Y,Z);
l=streamline(verts);
view(3);
set(l,'LineWidth',2);
set(l,'Color','k');
