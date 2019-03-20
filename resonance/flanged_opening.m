function Z = flanged_opening(f,rho,c,a)
  
  % returns acoustic impedance Z of opening though 
  % rigid flange (e.g., solid half-space); note that
  % properties rho,c are those of air into which sound
  % is radiated, not necessarily same as air below opening
  
  % omega = angular frequency []
  % rho = density [kg/m^3]
  % c = sound speed [m/s]
  % a = radius of opening/outlet [m]
  
  % equations referenced are from Rossing and Fletcher (2004)
  
  omega = 2*pi*f; % angular frequency
  k = omega/c; % wavenumber
  
  ka = k*a;
  
  R = ka.^2/2 - ka.^4/(2^2*3) + ka.^6/(2^2*3^2*4); % (8.29)
  X = 1/pi*(2^3*ka/3 - 2^5*ka.^3/(3^2*5) + 2^7*ka.^5/(3^2*5^2*7)); % (8.30)
    
  S = pi*a^2; 
  Z0 = rho*c/S;
 
  Z = Z0*(R+1i*X);
  
 
  
  
  
  
  