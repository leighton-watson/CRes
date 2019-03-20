%% TEST SOURCE FOURIER TRANSFORMS %%
%
% Compute fourier transforms of source numerically and compare to
% analytical expressions to validate analytical expressions

clear all;
clc;
path(pathdef);

%% Source Function Numerical %%

amp = [1]; % amplitude
L = [0.3]; % width
srcStyle = {'Gauss'}; % style
T = 10;
N = 2000;
resParams = [T N]; % resonance parameters

[S, f, s, ts] = sourceFunction(amp, L, srcStyle, resParams);

%% Plot Numerical Source Function %%

figure(1); clf;
subplot(2,1,1);
plot(ts, s); 

subplot(2,1,2);
plot(f(1:N/2+1), abs(S(1:N/2+1)));
xlim([0 5])

%% Analytical Source Function

omega = 2*pi*f(1:N/2+1);
Sa = amp*L*exp(-omega.^2./(2*L.^2));
subplot(2,1,2); hold on;
plot(omega/(2*pi), Sa,'-');
plot(omega, Sa,'-');
%plot(f(1:N/2+1), Sa,'-');
