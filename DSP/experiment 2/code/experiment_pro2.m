clc;
clear;
x = [3,11,7,0,-1,4,2];
nx = -3:3;
h = [2 3 0 -5 2 1];
nh = -1:4;
subplot(3,1,1);
stem(nh,h,'fill');
grid on;
xlabel("k");
ylabel("h(k)");
title("系统序列 h(k)");
subplot(3,1,2);
stem(nx,x,'fill');
grid on;
xlabel("k");
ylabel("x(k)");
title("输入序列 x(k)");
y = conv(x,h);
nn = nx(1) + nh(1):nx(end) + nh(end);

subplot(3,1,3);
stem(nn,y,'fill');
grid on;
xlabel("k");
ylabel("x(k)*h(k)");
title("零状态响应序列 x(k)*h(k)");
% [y ny] = Convolution(x,h,nx,nh);

%% 图像其余滤波方式

lena = imread('test_10.png');
J=imnoise(lena,'salt & pepper',0.02);
figure(2);
subplot(3,2,1);
imshow(lena);
title("原始图像");
subplot(3,2,2),
imshow(J);
title('添加盐椒噪声')

%% 中值滤波
b(:,:,1) = medfilt2(J(:,:,1),[3,3]);
b(:,:,2) = medfilt2(J(:,:,2),[3,3]);
b(:,:,3) = medfilt2(J(:,:,3),[3,3]);
subplot(3,2,3),
imshow(b);
title('中值滤波');
%% 拉氏算子锐化图像
h=[0,1,0;1,-4,0;0,1,0];%拉氏算子
c(:,:,1) = filter2(h,lena(:,:,1));
c(:,:,2) = filter2(h,lena(:,:,2));
c(:,:,3) = filter2(h,lena(:,:,3));
subplot(3,2,4),
imshow(c);
title('拉氏算子');

%% prewitt算子
H = fspecial('prewitt');%应用prewitt算子锐化图像
d(:,:,1) = filter2(H,lena(:,:,1));
d(:,:,2) = filter2(H,lena(:,:,2));
d(:,:,3) = filter2(H,lena(:,:,3));
subplot(3,2,5),
imshow(d);
title('prewitt算子');


%% log算子
H=fspecial('log');%应用log算子锐化图像
e(:,:,1) = filter2(H,lena(:,:,1));
e(:,:,2) = filter2(H,lena(:,:,2));
e(:,:,3) = filter2(H,lena(:,:,3));
subplot(3,2,6),
imshow(e);
title('log算子');