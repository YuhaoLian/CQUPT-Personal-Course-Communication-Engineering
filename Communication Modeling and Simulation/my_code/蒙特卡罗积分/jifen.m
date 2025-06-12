data=rand(100000000,3);
 x=4*data(:,1);
 y=-1+3*data(:,2);
 z=16*data(:,3);
 II=find(x>=y.^2&x<=y+2&z<=x.*(y.^2));
 M=length(II);
 V=192*M/100000000