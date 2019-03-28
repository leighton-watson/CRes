%%

clear all
clc
cmap = get(gca,'ColorOrder');

% set default plotting options
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',16);

% add file paths
addpath resonance
addpath SBPoperators
addpath geometry

figHand = figure(1); clf;
set(figHand,'Position',[100 100 1800 800]);

%% Resonance 1D Parameters %%
%
% Specify parameters for the solver.


T = 50; % total time [s]
N = 500; % number of grid points in time (formulas below assume even N)
dt = T/N; % time step
t = [0:N-1]*dt; % vector of time grid points [s]

Nyquist = 1/(2*dt); % Nyquist frequency [Hz]
freq = [0, Nyquist]; % range of frequencies [Hz]
Nf = N/2+1; % number of frequency samples

M = problemParameters(); % property values are defined in resonance/problemParameters.m
style = 'baffled piston'; % model of sound radiation ('baffled piston' or 'monopole')
order = 8; % order of numerical method

%% Source Function %%
%
% Here, a Gaussian source is assumed. But alternative source models can be
% used by replacing s = amp*exp(-0.5*(t-tshift).^2/L^2) with the desired
% source model in the time domain.


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
S = fft(s)*dt; % note multiplication by dt for correct amplitude and units!

S_pos = S(1:N/2+1); % positive sided Fourier spectra
f = fs(1:N/2+1); % vector of positive frequencies


%% Geometry %%
%
% Load the crater geometry. The crater geometry is required to be in a .mat
% file or a text file and is ordered as depth in first column starting from
% z=0 at the rim of the crater. The second column is the corresponding
% radius. The depth discretization must be regularly spaced. Depth
% increases downwards.


dmin = 20; % min depth
dmax = 140; % max depth
dz = 10; % depth increment
depths = dmin:dz:dmax; % range of depths


geom_str = {'Johnson2018','Richardson2014'};
line_style = {'-',':'};
alpha = [1 0.5];

for j = 1:length(geom_str)
    
    load(geom_str{j}); % load geometry
    shape = flipud(geometry); % reformat geometry matrix to be read by solver
        
    % plot crater geometry
    subplot(3,5,[1 6 11]);
    h = plot(shape(:,2), shape(:,1),'Color','k','LineStyle',line_style{j});
    h.Color(4) = alpha(j);
    set(gca,'YDir','Reverse');
    hold on; box on;
    h = plot(-shape(:,2), shape(:,1),'Color','k','LineStyle',line_style{j});
    h.Color(4) = alpha(j);
    ylim([min(depths) max(depths)]);
    xlabel('radius (m)');
    ylabel('depth (m)');
    grid on
    set(h.Parent,'YMinorGrid','on');
    xlim([-75 75]);
    
    
    %% Iterate Over Depths %%
    %
    % Iteratre through different depth values. Calculate the resonant
    % frequency and quality factor at each depth.
    
    for i = 1:length(depths)
        
        depth = depths(i);
        str = strcat('Depth=',num2str(depth),'m');
        disp(str); % display string to keep track of iterations while they run
        
        A = resonance1d(shape, depth, freq, Nf, style, order, M); % solve transfer function
        
        transFunc(i,:) = A.P;
        
%         % normalize transfer function by amplitude of first peak
%         [Pmax, ipeak] = extrema(A.P); % find amplitude and index of each maxima
%         ipeak = sort(ipeak); % sort indices into ascending order
%         Pmax_first = A.P(ipeak(1)); % find amplitude of the first peak
%         %transFunc(i,:) = transFunc(i,:)./Pmax_first;
        
        spectra_save(i,:) = transFunc(i,:) .* S_pos;
        spectra = spectra_save(i,:);
        
        B.P = spectra;
        B.f = A.f;
        
        [ymax,imax,ymin,imin] = extrema(abs(spectra));
        [imin,min_sort_index] = sort(imin);
        ymin = ymin(min_sort_index);
        [imax,max_sort_index] = sort(imax);
        ymax = ymax(max_sort_index);
        
        flow(i) = interp1(abs(spectra(1:imax(1))),A.f(1:imax(1)),ymax(1)/sqrt(2));
        fhigh(i) = interp1(abs(spectra(imax(1):imin(2))),A.f(imax(1):imin(2)),ymax(1)/sqrt(2));
        if isnan(fhigh(i))
            fhigh(i) = A.f(imax(1)) + (A.f(imax(1))-flow(i));
        end
        f0(i) = (flow(i)+fhigh(i))/2;
        bw(i) = fhigh(i)-flow(i);
        Q(i) = f0(i)/bw(i);
    end
    
    % plot resonant frequency
    subplot(3,5,[2 7 12]);
    h1 = plot(f0, depths,'Color','k','LineStyle',line_style{j});
    set(gca,'YDir','Reverse');
    h1.Color(4) = alpha(j);
    xlabel('frequency (Hz)');
    set(h1.Parent,'YTickLabel',{});
    set(h1.Parent,'YMinorGrid','on');
    set(h1.Parent,'XTick',[0.4 0.5 0.6 0.7 0.8]);
    hold on;
    grid on;
    ylim([min(depths) max(depths)])
    xlim([0.4 0.8])
    
    
    % plot quality factor
    subplot(3,5,[3 8 13]);
    h2 = plot(Q, depths,'Color','k','LineStyle',line_style{j});
    set(gca,'YDir','Reverse');
    h2.Color(4) = alpha(j);
    xlabel('quality factor');
    set(h2.Parent,'YTickLabel',{});
    set(h2.Parent,'YMinorGrid','on');
    set(h2.Parent,'XTick',[0 2 4 6 8 10]);
    hold on;
    grid on
    ylim([min(depths) max(depths)])
    
    
    
    %% Example signal in time and frequency domain for several depths
    
    depth_plot = [50 80 120];
    
    for i = 1:length(depth_plot)
        
        A = resonance1d(shape, depth_plot(i), freq, Nf, style, order, M);
        
        transFunc = A.P;
        
%         % normalize transfer function by amplitude of first peak
%         [Pmax, ipeak] = extrema(A.P); % find amplitude and index of each maxima
%         ipeak = sort(ipeak); % sort indices into ascending order
%         Pmax_first = A.P(ipeak(1)); % find amplitude of the first peak
%         transFunc_norm = transFunc./Pmax_first;
        
        spectra = transFunc .* S_pos;
        
        subplot(3,5,5*i);
        h3 = plot(A.f, abs(spectra)./max(abs(spectra)),'Color',cmap(i,:),'LineStyle',line_style{j});
        xlim([0 2]); hold on;
        set(h3.Parent,'YTick',[]);
        set(h3.Parent,'XTick',[0 0.5 1 1.5 2]);
        grid on
        h3.Color(4) = alpha(j);
        
        % invert to time domain
        sig_pos = A.P(1:N/2+1) .* S(1:N/2+1);
        sig_full = [sig_pos conj(sig_pos(end-1:-1:2))]; % reflect about f = 0. Take complex conjugate for negative frequencies
        sig_full(N/2+1) = real(sig_full(N/2+1)); % entry at Nyquist must be real
        sig_time = ifft(sig_full,'symmetric')/dt; % inverse Fourier transform to obtain Greens function in time domain
        
        subplot(3,5,5*i-1);
        h4 = plot(t-13, sig_time./max(sig_time),'Color',cmap(i,:),'LineStyle',line_style{j});
        xlim([0 10]); hold on;
        set(h4.Parent,'YTick',[]);
        grid on
        ylim([-1.8 1.5])
        h4.Color(4) = alpha(j);
        
        %%% add depth indicators to other plots %%%
        % geometry
        idx = find(A.geometry(:,1)==depth_plot(i));
        subplot(3,5,[1 6 11]);
        %h5 = 
        plot([-A.geometry(idx,2) A.geometry(idx,2)],[depth_plot(i) depth_plot(i)],...
            'Color',cmap(i,:),'LineStyle','-');
        %h5.Color(4) = alpha(j);
        
        % frequency
        idx2 = find(depths==depth_plot(i));
        subplot(3,5,[2 7 12]);
        %h6 = 
        plot([0 f0(idx2)],[depth_plot(i) depth_plot(i)],...
            'Color',cmap(i,:),'LineStyle','-');
        %h6.Color(4) = alpha(j);
        
        % quality factor
        subplot(3,5,[3 8 13]);
        %h7 = 
        plot([0 Q(idx2)],[depth_plot(i) depth_plot(i)],...
            'Color',cmap(i,:),'LineStyle','-');
        %h7.Color(4) = alpha(j);
        
    end
    
    subplot(3,5,14);
    xlabel('time (s)');
    
    subplot(3,5,15);
    xlabel('frequency (Hz)');
    
end


