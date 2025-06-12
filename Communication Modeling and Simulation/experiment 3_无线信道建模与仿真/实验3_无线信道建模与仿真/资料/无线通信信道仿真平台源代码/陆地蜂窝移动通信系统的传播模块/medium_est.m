clc;
close all;
clear;
tic;

area_id=1;

f=600e6;

Dradius=200000;

scale=100;

shadow_fading(Dradius, f, area_id, scale);

TimeElapse = 1000*toc;	% Simulation minutes
fprintf('****** The total time used is : %3.16f milliseconds ********\n', TimeElapse);

out_fading = load('shadow_fading.bin');
save('shadow_fading.mat','out_fading');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Estimation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hist(out_fading,100);
mean_value=mean(out_fading)
tao=sum((out_fading-mean_value).^2)/length(out_fading)

Dradius_tmp=Dradius/scale;
for i=1:Dradius_tmp
    out_fading_est(i,:)=out_fading(1,((i-1)*Dradius_tmp+1):1:i*Dradius_tmp);
end
figure(1)
meshc(1:Dradius_tmp,1:Dradius_tmp,out_fading_est);
figure(2)
meshc(1:Dradius,1:Dradius,out_fading1);
