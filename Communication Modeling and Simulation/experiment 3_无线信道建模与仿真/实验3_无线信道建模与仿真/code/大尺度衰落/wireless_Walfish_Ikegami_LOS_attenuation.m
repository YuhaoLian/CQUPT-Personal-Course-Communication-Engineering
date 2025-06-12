function y=wireless_Walfish_Ikegami_LOS_attenuation(f,d)
    y=42.6+26*log(d)/log(10)+20*log(f)/log(10);