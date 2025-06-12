% From Jakes, Microwave Mobile Communications, 1973 IEEE Press.

function  [fadingcoeff]=multipath(infrom_length,Fc,V,Fs,path_num)

C=3e8;                                %light speed Km/s
Fm=(V*Fc*1000)/(C*3600);              %Doppler shift
Ts = 1/Fs;                            %Chip duration according to digital bandwidth Fs
Wm=2*pi*Fm;
Ts_adjustment=7.5+Ts*infrom_length*rand(1);
t=Ts+Ts_adjustment:Ts:Ts*infrom_length+Ts_adjustment;   % Time array

M=40;
N = 2*(2*M+1);
n = 1 : M;
f = Fm*cos(2*pi*n/N);	    %Frequency vector

fadingcoeff=zeros(path_num,infrom_length);

for ii=1:path_num
    
	%Initial phases
	alpha = 0;
	beta_n = pi*(n + 2*(ii-1))/(M+1); 	%Phase vector
	
	%real part
	Xc0 = sqrt(2)*cos(alpha)*cos(Wm*t);
	Xc = Xc0 + 2*cos(beta_n)*cos(2*pi*f'*t);
	%imaginary part
	Xs0 = sqrt(2)*sin(alpha)*cos(Wm*t);
	Xs = Xs0 + 2*sin(beta_n)*cos(2*pi*f'*t);
	
	%complex fading function
	fadingcoeff(ii,:) = (1/sqrt(2*M+1))*(Xc+sqrt(-1)*Xs);

end