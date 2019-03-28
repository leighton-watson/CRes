function [S, f, s, t] = sourceFunction(A, L, srcStyle, resParams)
% [S, f, s, t] = sourceFunction(A, L, srcStyle, resParams)
%
% Compute source function in the time and frequency domains.
% 
% INPUT
% A = source amplitude
% L = source width
% srcStyle = style of source mechanism ('Gauss' or 'Brune')
% resParams = parameters to ensure consistency with resonance1d simulation
%
% OUTPUT
% S = source function in frequency domain
% f = frequency vector
% s = source function in time domain
% t = time vector

T = resParams(1);
N = resParams(2);
dt = T/N;
t = [0:N-1]*dt;
tshift = T/4; % shift location of center of Gaussian away from zero

if strcmp(srcStyle,'Gauss')
    s = A*exp(-0.5*(t-tshift).^2/L^2); % source function in time domain
    
elseif strcmp(srcStyle,'Brune')
    s = A.*(t-tshift).*exp(-(t-tshift)./L) .* heaviside(t-tshift);
    
else
    error('Incorrect source style. Please select from ''Gauss'' or ''Brune''');
    
end
      
omega = [0:N/2 (-N/2+1):-1]/N*2*pi/dt; % angular frequency vector
f = omega/(2*pi); % frequency vector
S = fft(s)*dt; % source function in frequency domain

