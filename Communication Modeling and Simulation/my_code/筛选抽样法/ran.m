clc;
clear;

R0=2;
R1=5;

[xRandnum,yRandnum,zRandnum] = Screening_sampling(R0,R1,100);
subplot(2,2,1);
scatter3(xRandnum,yRandnum,zRandnum,'.');
xlim([-5,5]);
title('100个散点');

[xRandnum,yRandnum,zRandnum] = Screening_sampling(R0,R1,1000);
subplot(2,2,2);
scatter3(xRandnum,yRandnum,zRandnum,'.');
xlim([-5,5]);
title('1000个散点');

[xRandnum,yRandnum,zRandnum] = Screening_sampling(R0,R1,10000);
subplot(2,2,3);
scatter3(xRandnum,yRandnum,zRandnum,'.');
xlim([-5,5]);
title('10000个散点');

[xRandnum,yRandnum,zRandnum] = Screening_sampling(R0,R1,100000);
subplot(2,2,4);
scatter3(xRandnum,yRandnum,zRandnum,'.');
xlim([-5,5]);
title('100000个散点');