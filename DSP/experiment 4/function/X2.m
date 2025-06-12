function y= X2(n)
%X2 信号2
if n >= 0&&n <= 3
    y = n + 1;
elseif n >= 4&&n <= 7
    y = 8 - n;
else
    y = 0;
end 
end