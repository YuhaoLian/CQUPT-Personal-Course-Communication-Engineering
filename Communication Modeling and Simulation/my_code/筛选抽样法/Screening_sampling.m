function  [xRandnum,yRandnum,zRandnum]= Screening_sampling(R0,R1,mm)
xRandnum=zeros(1,mm);
yRandnum=zeros(1,mm);
zRandnum=zeros(1,mm);
ii=1;
while ii<mm
      Randnum1=unifrnd(-5,5); 
      Randnum2=unifrnd(-5,5); 
      Randnum3=unifrnd(-5,5); 
      s=Randnum1^2+Randnum2^2+Randnum3^2;
      if  (s<=R1^2)&&(s>=R0^2)
       xRandnum(1,ii)=Randnum1;
       yRandnum(1,ii)=Randnum2;
       zRandnum(1,ii)=Randnum3; 
       ii=ii+1;
     end
end


  
