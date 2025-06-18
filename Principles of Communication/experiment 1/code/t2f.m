%∏µ¿Ô“∂±‰ªª
function S=t2f(s,NFFT)
S = fft(s, NFFT);
S = fftshift(S);
end
