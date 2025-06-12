function xR = Composite_sampling(mm)
R=unifrnd(0,1,mm,1);R1=exprnd(1/3,mm,1);
R2=exprnd(2,mm,1);xR=zeros(mm,1);
for ii=1:mm
    if R(ii,1)<=0.4
       xR(ii,1)=R1(ii,1);
    else
       xR(ii,1)=R2(ii,1); 
   end
end   

