function ray = jakes(fm,fs,M,N_0,index)
% Jakes Model of a Rayleigh fading channel
% Written by Chandra Athaudage 12/10/2001, CUBIN
%
% jakes(fm,fs,M,N_0,index)
%
% fm Doppler frequency 
% fs sampling frequency
% M  number of samples
% N_0 number of paths
% index (between 1 & N_0) uncorrelated Rayleigh fading functions

N = 2*(2*N_0+1);
nn = 1 : N_0;
f = fm*cos(2*pi*nn/N);	    %frequency vector
Ts = 1/fs; 					%sampling interval (discrete signals)
onset = N_0*(1/fm);		    %onset before steady model output
Te = 2*onset + M*Ts; 	    %evaluation period
t = 0 : Ts : Te;			%time vector

%%%%%%%%%%%%%%%   
%Initial phases
alpha = 0;
beta_n = pi*(nn + 2*(index-1))/(N_0+1); 	%phase vector

%%%%%%%%%%%
%simulation
%%%%%%%%%%%
%real part
Xc0 = sqrt(2)*cos(alpha)*cos(2*pi*fm*t);
Xc = Xc0 + 2*cos(beta_n)*cos(2*pi*f'*t);
%imaginary part
Xs0 = sqrt(2)*sin(alpha)*cos(2*pi*fm*t);
Xs = Xs0 + 2*sin(beta_n)*cos(2*pi*f'*t);

%%%%%%%%%%%%%%%%%%%%%%%%
%complex fading function
T  = (1/sqrt(2*N_0+1))*(Xc+sqrt(-1)*Xs);
%%%%%%%%%
Tstart = onset + 0.9*rand(1)*onset; %avoid onset and random start 
Mstart = round(Tstart/Ts);
Mend = Mstart + M -1;
ray = T(Mstart:Mend);