% File: multipath_sim.m
% Matlab Simulink IEEE 802.15.3a compliant channel model 
% Author: Tim Becker (MSEE Degree, December 2004)
% Advisor: Robert Morelos-Zaragoza. San Jose State University

function [sys,x0,str,ts] = multipath_sim(t,x,u,flag,T,nT,LAMBDA,lambda,GAMMA,gamma,sigma1_dB,sigma2_dB,std_shadow,LOSflag)
% Differs from multipath_sim2.m in that h begins at t=0 for NLOS channels.

v = 10000;
switch flag,

  % Initialization
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes(T,nT,v);

  % Outputs
  case 3,
    sys=mdlOutputs(t,x,u,T,LAMBDA,lambda,GAMMA,gamma,sigma1_dB,sigma2_dB,std_shadow,LOSflag,v);
    
  % Unhandled flags
  case {1, 2, 4, 9}
      sys = [];

  % Unexpected flags
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);

end

% end multipath_sim



%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================

function [sys,x0,str,ts]=mdlInitializeSizes(T,nT,v)

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = v;
sizes.NumInputs      = 0;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

% initialize the initial conditions
x0  = [];

% str is always an empty matrix
str = [];

% initialize the array of sample times
ts  = [T*nT*1e-9 0];

% end mdlInitializeSizes



%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================

function sys=mdlOutputs(t,x,u,T,LAMBDA,lambda,GAMMA,gamma,sigma1_dB,sigma2_dB,std_shadow,LOSflag,v)

%   returns the real impulse response of the
%   indoor multipath channel as modeled by Intel in
%   IEEE 802.15-02/279r0-SG3a and is based on the
%   Saleh-Valenzuela model with lognormal fading.
%   T is the time resolution (nsec).

sigma1 = 10^(sigma1_dB/20);
sigma2 = 10^(sigma2_dB/20);
Omega0 = 1;

if (strcmp(LOSflag,'on'))
    T1 = [0];
else
    T1 = [randn^2/(2*LAMBDA)+randn^2/(2*LAMBDA)];
end

while T1(end) <= 10*GAMMA
    dT = randn^2/(2*LAMBDA)+randn^2/(2*LAMBDA);
    T1 = [T1 T1(end)+dT];
end

for n = 1:length(T1)
    tau(n) = {[T1(n)]};
    while tau{n}(end) <= tau{n}(1)+10*gamma
        dt = randn^2/(2*lambda)+randn^2/(2*lambda);
        tau(n) = {[tau{n} tau{n}(end)+dt]};
    end
    
    mu(n) = {(10*log(Omega0)-10*T1(n)/GAMMA-10*(tau{n}-T1(n))/gamma)/log(10)-(sigma1^2+sigma2^2)*log(10)/20};
    beta(n) = {10.^((sigma1*randn(size(mu{n}))+sigma2*randn(size(mu{n}))+mu{n})/20)};
    alpha(n) = {beta{n}.*sign(rand(size(beta{n}))-0.5)};  % Intel version, +/-1
    %alpha(n) = {beta{n}.*exp(-j*2*pi*rand(size(beta{n})))};  % cell array of channel coefficients
end

maxtime = max(tau{1});  % find maximum arrival time index in cell array tau
for n = 2:length(T1)
    if max(tau{n}) > maxtime
        maxtime = max(tau{n});
    end
end

N = 32;  % oversampling factor
Nfs = N/T;
h = zeros(1,floor(maxtime*Nfs)+1);  % initialize impulse response vector
for n = 1:length(T1)
    tauN{n} = floor(tau{n}*Nfs);  % quantized time indices
end

for n = 1:length(T1)
    for m = 1:length(tau{n})
        h(tauN{n}(m)+1-tauN{1}(1)) = h(tauN{n}(m)+1-tauN{1}(1))+alpha{n}(m);  % only line different from multipath_sim2.m
    end
end

h = N*resample(h,1,N);
maxtime = ceil((10*GAMMA+10*gamma)/T);  % maximum arrival time for channel
h = h(1:maxtime);  % concatenate h to maximum channel length

E = sum(h.*conj(h));  % compute total channel energy
h = h./sqrt(E);  % normalize total energy to 1

fac = 10^(std_shadow*randn/20);
h = h.*fac;

if(length(h)<v)
    h = [h zeros(1,v-length(h))];
elseif(length(h)>v)
    h = h(1:v);
end
%length(h)
sys = real(h);

% end mdlOutputs

