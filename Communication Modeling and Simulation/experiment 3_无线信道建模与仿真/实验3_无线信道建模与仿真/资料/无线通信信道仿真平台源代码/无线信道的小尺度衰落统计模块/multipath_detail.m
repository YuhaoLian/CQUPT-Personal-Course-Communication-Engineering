function  [Path_Delay,Path_loss,Doppler_Frequency]=multipath_detail(fs)

Path_loss = [-9, -6, -3, 0, 0, 0, 0, 0, -3.5, -6.5, -9.5, -12.5];%dB
Path_Delay = [-25, -15, -5, -0.3, -0.1, 0, 0.12, 0.33, 5.5, 17, 28, 35];%us
Doppler_Frequency = [0.1, 0.3, 0.5, 1, 10, 0.05, 15, 2.1, 0.9, 0.6, 0.27, 0.11 ];%Hz

% Path_loss = [0,-1,-9,-10,-15,-20];
% 
% Path_Delay = [0,2,4,6,9,23]*1e6/fs;

tao=((10.^(Path_loss/10))*Path_Delay')/sum(10.^(Path_loss/10));
tao2=((10.^(Path_loss/10))*((Path_Delay).^2)')/sum(10.^(Path_loss/10));
Bc=1e6/(5*sqrt(tao2-tao^2))

Path_Delay = round((Path_Delay-min(Path_Delay))*fs/1e6);
% Path_loss = [0,-1,-9,-10,-15,-20];
% Path_Delay = [0,2,4,6,9,23];
% Doppler_Frequency = [0,0,0,0,0,0];
% Path_loss = [0, 0, 0, -3.5, -6.5, -9.5, -12.5];%dB
% Path_Delay = [0, 0.12, 0.33, 5.5, 17, 28, 35];%us
% Doppler_Frequency = [0.05, 15, 2.1, 0.9, 0.6, 0.27, 0.11 ];%Hz
% Path_Delay = round((Path_Delay-min(Path_Delay))*fs/1e6);