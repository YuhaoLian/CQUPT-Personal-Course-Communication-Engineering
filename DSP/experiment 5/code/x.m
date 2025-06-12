function y = x(n)
    if (n >= 0)&&(n <= 13)
        y = n + 1;
    elseif n >= 14&&n <= 27
        y = 27 - n;
    else
        y = 0;
    end
end

