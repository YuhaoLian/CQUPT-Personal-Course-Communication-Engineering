function [y,ny] = Convolution(x,h,n_x,n_h)
%Convolution 求卷积和并进行画图
%     x:输入序列
%     h:系统序列
%     n_x:输入序列n范围
%     n_h:系统序列n范围   

[nx,mx] = size(x);
for i = 1:mx
    if n_x(i) >= 0
        flag = i;
        break;
    end
end
X = x(flag:mx); 
nX = 0:1:mx-flag;
[nh,mh] = size(h);
for i = 1:mh
    if n_h(i) >= 0
        flag = i;
        break;
    end
end
H = h(flag:mh);
nH = 0:1:mh-flag;
ny1 = nX(1) + nH(1);
ny2 = nX(end) + nH(end);
ny = ny1:ny2;
y = conv(X,H);


figure(1);
subplot(3,1,1);
stem(n_h,H,'fill');
grid on;
xlabel("k");
ylabel("h(k)");
title("系统序列 h(k)");

subplot(3,1,2);
stem(n_x,X,'fill');
grid on;
xlabel("k");
ylabel("x(k)");
title("输入序列 x(k)");



subplot(3,1,3);
stem(ny,y,'fill');
grid on;
xlabel("k");
ylabel("x(k)*h(k)");
title("零状态响应序列 x(k)*h(k)");
end

