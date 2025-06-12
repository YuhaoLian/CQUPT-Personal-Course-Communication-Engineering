clc
close all
clear

area_id=1;

hB=30;

hT=10;

f=600e6;

X1=101;

X2=101;
% 
% out_fading1=zeros(X1 ,X2 ); 
% 
% out_fading2=zeros(X1 ,X2 );

% for i = 1:X1
%     for j = 1:X2
%         d1=abs(i-1-(X1-1)/2);
%         d2=abs(j-1-(X2-1)/2);
%         d=sqrt(d1^2+d2^2);
%         out_fading1(j,i)=-long_term_fading(d, f, area_id, hB, hT);
%     end
% end

% out_fading2=out_fading1+shadow_fading(X1, X2, f, area_id);
out_fading2=shadow_fading(X1, X2, f, area_id);
figure(1)
meshc(1:X1,1:X2,out_fading2);
figure(2)
contourf(1:X1,1:X2,out_fading2);
% colorbar
% colormap(bone)

% plot(1:X2,out_fading1(501,:))
% hold on
% plot(1:X2,out_fading2(101,:))