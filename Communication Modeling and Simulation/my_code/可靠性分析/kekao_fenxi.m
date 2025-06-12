function Rguji=kekao_fenxi(t,thetaa1,thetaa2,thetab1,thetab2,mm)
%t 是要求系统生存的寿命%thetaa1 是元件A1的数学期望
%thetaa2 是元件A2的数学期望%thetab1 是元件B1的数学期望 
%thetab2 是元件B2的数学期望%mm 是随机实验次数
frq=0;randnuma1 = exprnd(thetaa1,1,mm);
randnuma2 = exprnd(thetaa2,1,mm);
randnumb1 = exprnd(thetab1,1,mm);
randnumb2 = exprnd(thetab2,1,mm);
for ii=1:mm
    if (randnuma1(1,ii)>t)&(randnuma2(1,ii)>t)        pass1=1;
    else         pass1=0;
    end 
    if (randnumb1(1,ii)>t)&(randnumb2(1,ii)>t)         pass2=1;
    else        pass2=0;
    end 
    if (pass1+pass2)>=1        frq=frq+1;
    end    
end,Rguji=frq/mm