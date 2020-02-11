%% MAKE FIG %%

clear all
clc
cmap = get(gca,'ColorOrder');

% set default plotting options
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',16);

% add file paths
addpath ../source/resonance
addpath ../source/SBPoperators
addpath Functions/

figHand = figure(1); clf;
%set(figHand,'Position',[100 100 1000 280]);


%% Resonance 1D Parameters %%

T = 100; % total time [s]
N = 2000; % number of grid points in time (formulas below assume even N)
dt = T/N; % time step
t = [0:N-1]*dt; % vector of time grid points [s]

Nyquist = 1/(2*dt); % Nyquist frequency [Hz]
freq = [0, Nyquist]; % range of frequencies [Hz]
Nf = N/2+1; % number of frequency samples

style = 'baffled piston'; % model of sound radiation ('baffled piston' or 'monopole')
order = 8; % order of numerical method

%% Source Function %%

amp = 1; % amplitude
L = 0.3; % width (standard deviation) of Gaussian (s)
tshift = T/4; % location of center of Gaussian
s = amp*exp(-0.5*(t-tshift).^2/L^2); % source time function

% Fourier transform
omega = [0:N/2 (-N/2+1):-1]/N*2*pi/dt; % angular frequency (1/s)
% note storage order: starts with omega=0, increases to
% Nyquist angular frequency omegaNyq=pi/dt (frequency fNyq=1/(2*dt)),
% then starts again at negative frequencies (excluding -omegaNyq and 0)
fs = omega/(2*pi); % frequency (Hz)

% now Fourier transform
S = fft(s)*dt;
% note multiplication by dt for correct amplitude and units!
S_pos = S(1:N/2+1);
f = fs(1:N/2+1);

%% Geometry %%
% crater properties
radius = 100; % crater radius
depth = 300; % crater depth
z = [depth:-1:0]'; % depth vector
shape = [z radius*ones(size(z))]; % crater geometry matrix


dmin = 20; % min depth
dmax = depth; % max depth
dz = 10; % depth increment
depths = dmin:dz:dmax; % range of depths

line_style = {'-','-','-'};
alpha = [1 0.9 0.8];
Tmin = [100 400 100]; % temp at top
Tmax = [100 400 500]; % temp at bottom
for j = 1:length(Tmin)
    
    
    % temperature profile
    T_profile = (Tmax(j)-Tmin(j))/depth.*z + Tmin(j);
    
    A = resonance1d_Temp(shape, depth, freq, Nf, style, order, Tmin(j), Tmax(j));
    
    transFunc = A.P;   
    spectra = transFunc .* S_pos;
    
    B.f = A.f;
    B.P = spectra;
    
    [f0(j) Q(j)] = resPeakProps(B, 'sim');
    
    % plot signal in frequency domain
    subplot(1,3,3);
    if j == 1
        spectra_max = max(abs(spectra));
    end
    h3 = plot(A.f, abs(spectra)./spectra_max,'LineStyle',line_style{j});
    h3.Color(4) = alpha(j);
    set(h3.Parent,'YTick',[]);
    xlim([0 1.5]); hold on;
    ylim([0 1.2])
    xlabel('Frequency (Hz)');
    ylabel('$\Delta p(\omega,r)$','interpreter','latex')
    
    % invert to time domain
    sig_pos = A.P(1:N/2+1) .* S(1:N/2+1);
    sig_full = [sig_pos conj(sig_pos(end-1:-1:2))]; % reflect about f = 0. Take complex conjugate for negative frequencies
    sig_full(N/2+1) = real(sig_full(N/2+1)); % entry at Nyquist must be real
    sig_time = ifft(sig_full,'symmetric')/dt; % inverse Fourier transform to obtain Greens function in time domain
    
    % plot signal in time domain
    subplot(1,3,[1 2]);
    if j == 1
        p_max = max(abs(sig_time));
    end
    h4 = plot(t-23, sig_time./p_max,'LineStyle',line_style{j});
    h4.Color(4) = alpha(j);
    set(h4.Parent,'YTick',[]);
    xlim([0 50]); hold on;
    ylim([-0.9 1.2])
    xlabel('Time (s)');
    ylabel('$\Delta p(t,r)$','interpreter','latex')
    
end


