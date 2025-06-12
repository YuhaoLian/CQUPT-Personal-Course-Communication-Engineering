function  [LNLOS]=long_term_fading(d, f, area_id, hB, hT)

% d: transmitter-receiver-distance in meters
% f: carrier frequency (Hz)
% area_id: 1(small/medium cities), 2(large cities), 3(suburban), 4(open areas)
% hB: height of Base antenna
% hT: height of CPE antenna

C=3e8;%light speed m/s
f0=1e9;
GB=12;%Base antenna gain(dBi) 
GT=10;%CPE antenna gain(dBi) 
uexcess=6;

LLOS=20*log10(4*pi/C)+20*log10(d)+20*log10(f)-10*log10(GB)-10*log10(GT);

[LHATA uHATA]=excess_path_loss(area_id, hB, hT);

Lexcess=LHATA+uHATA*log10(d)+uexcess*log10(f/f0);

LNLOS=LLOS+Lexcess;