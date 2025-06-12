clc;
close all;
clear;
tic;

area_id=1;

hB=100;

hT=10;

f=600e6;

for d =1:200000
    d
    path_lose(d)=-long_term_fading(d, f, area_id, hB, hT);
end

TimeElapse = 1000*toc;	% Simulation minutes
fprintf('****** The total time used is : %3.16f milliseconds ********\n', TimeElapse);

save('path_lose.mat','path_lose');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Estimation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(1:d, path_lose(1,:));
