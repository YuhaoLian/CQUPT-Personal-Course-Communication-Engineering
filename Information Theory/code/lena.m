clc;
clear;

lena_img = imread('lena.bmp');             
figure();                    
imshow(lena_img);

figure();     
imhist(lena_img);
count = imhist(lena_img);

p = count/sum(count);

h = -p.*log2(p);
h(find(isnan(h)==1)) = 0;

HH = sum(h);
H = entropy(lena_img);