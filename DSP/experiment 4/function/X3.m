function y= X3(n)
%X3 信号3
if n>=0&&n<=3
    y = 4 - n;
elseif n>=4&&n<=7
    y = n - 3;
else
    y = 0;
end 
end