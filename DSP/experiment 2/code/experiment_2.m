clc;
clear;

a = [3 4 1];
b = [1 1];
figure(1);
subplot(2,1,1);
impz(b,a,30);
grid on;
title("系统取样响应");

c = [2.5 6 10];
d = [1];
subplot(2,1,2);
impz(d,c,30);
grid on;
title("系统取样响应");

nh = -3:14;
nx = -3:9;
h = ((7/8).^nh).*(uDT(nh)-uDT(nh-10));
x = uDT(nx)-uDT(nx-5);
% 求卷积
y = conv(x,h);
ny1 = nx(1) + nh(1);
ny2 = nx(end) + nh(end);
ny = ny1:ny2;

figure(2);
subplot(3,1,1);
stem(nh,h,'fill');
grid on;
xlabel("k");
ylabel("h(k)");
title("h(k)");

subplot(3,1,2);
stem(nx,x,'fill');
grid on;
xlabel("k");
ylabel("x(k)");
title("x(k)");

subplot(3,1,3);
stem(ny,y,'fill');
grid on;
xlabel("k");
ylabel("x(k)*h(k)");
title("x(k)*h(k)");


%% 图像
lena = imread('test_10.png');
figure(3);
subplot(2,2,1);
imshow(lena);
title("原始图像");
Gx = [-1 0 1;-2 0 2;-1 0 1];
Gy = [1 2 1;0 0 0;-1 -2 -1];

lena_gx(:,:,1) = conv2(Gx,lena(:,:,1));
lena_gx(:,:,2) = conv2(Gx,lena(:,:,2));
lena_gx(:,:,3) = conv2(Gx,lena(:,:,3));
subplot(2,2,2);
imshow(lena_gx);
title("采用 Gx 进行卷积滤波后的图像");


lena_gy(:,:,1) = conv2(Gy,lena(:,:,1));
lena_gy(:,:,2) = conv2(Gy,lena(:,:,2));
lena_gy(:,:,3) = conv2(Gy,lena(:,:,3));
subplot(2,2,3);
imshow(lena_gy);
title("采用 Gy 进行卷积滤波后的图像");

lena_gxy(:,:,1) = (lena_gy(:,:,1) + lena_gx(:,:,1))/2;
lena_gxy(:,:,2) = (lena_gy(:,:,2) + lena_gx(:,:,2))/2;
lena_gxy(:,:,3) = (lena_gy(:,:,3) + lena_gx(:,:,3))/2;

% lena_gxy(:,:,1) = sqrt(lena_gy(:,:,1).* lena_gy(:,:,1) + lena_gx(:,:,1).* lena_gx(:,:,1));
% % lena_gxy(:,:,2) = sqrt(lena_gy(:,:,2).* lena_gy(:,:,2) + lena_gx(:,:,2).* lena_gx(:,:,2));
% lena_gxy(:,:,3) = sqrt(lena_gy(:,:,3).* lena_gy(:,:,3) + lena_gx(:,:,3).* lena_gx(:,:,3));


subplot(2,2,4);
imshow(lena_gxy);
title("采用 Gx 和 Gy 进行卷积滤波后的图像");

