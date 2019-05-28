%% MAKE FIG %%
% VALIDATION ANALYTIC 
%
% Compare numerical simulation of infrasound signal from a pipe with the
% analytical solutions for resonant frequency 

clear all; clc; %close all;
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',20);
cmap = get(gca,'ColorOrder');

j = 1; % counter

% add paths
addpath ../source/SBPoperators/
addpath ../source/resonance/
addpath Functions/

% create figure handle
%figHand = figure(1); clf;
figure(1);
%set(figHand,'Position',[100 100 1200 600]);

% ambient properties
gamma = 1; % isothermal to match analytical solutions
r = 1000; % distance from source to receiver
M = problemParameters_validation(gamma,r); % load problem parameters

% parameters for numerical model
freq = [0 1]; % frequency range [Hz]
Nf = 500; % number of frequency samples
order = 8; % order of discretization operators
style = 'baffled piston'; % acoustic source description

%%% Crater Properties %%%
a = 100;
depth = [100 150 200 250 300 350 400 450 500];
depth_plot = [200 300 400];

for k = 1:length(a)
    for i = 1:length(depth)
        depth(i)
        radius = a(k);
        % create geometry matrix
        z = [depth(i):-1:0]'; % depth vector
        r = radius*ones(size(z)); % radius vector
        shape = [z r]; % geometry matrix
        
        %%% Numerical Solution %%%
        
        % compute infrasound signal
        A = resonance1d(shape, depth(i), freq, Nf, style, order, M);
        
        % calculate resonant properties
        [Y,I] = max(abs(A.P));
        Is(i) = I;
        Ys(i) = Y;
        flow(i) = interp1(abs(A.P(1:I)),A.f(1:I),Y*sqrt(2)/2);
        fhigh(i) = interp1(abs(A.P(I:2*I)),A.f(I:2*I),Y*sqrt(2)/2);
        flow99(i) = interp1(abs(A.P(1:I)),A.f(1:I),Y*.99);
        fhigh99(i) = interp1(abs(A.P(I:2*I)),A.f(I:2*I),Y*.99);
        fpeak(i) = (flow99(i)+fhigh99(i))/2;
        f0(i) = (flow(i)+fhigh(i))/2;
        bw(i) = fhigh(i)-flow(i);
        Q(i) = f0(i)/(bw(i));
        
        
        %%% Analytical Solution %%%
        Leff = depth(i) + 8/(3*pi)*radius; % effective length of crater
        c = M.cC; % speed of sound in crater
        f0_analyticaleff(i) = c/(4*Leff); % resonant frequency
        f0_analytical(i) = c/(4*depth(i));
        
        if ismember(depth(i),depth_plot)
            subplot(1,4,[2 3]);
            h1 = plot(A.f, abs(A.P)./max(abs(A.P)),'Color',cmap(j,:));
            hold on;
            plot([f0_analyticaleff(i) f0_analyticaleff(i)],[0 1],'Color',cmap(j,:),'LineStyle','--')
            
            subplot(1,4,1);
            h2 = plot([0 1],[depth(i) depth(i)],'Color',cmap(j,:),'LineStyle',':');
            h2.Color(4) = 0.5;
            hold on;

            j = j + 1;
        end
        
    end
    
    %% Plotting %%
    subplot(1,4,[2 3]);
    xlabel('Frequency (Hz)');
    set(h1.Parent,'YTick',[]);
    set(h1.Parent,'XTick',[0.1 0.15 0.2 0.25 0.3 0.35 0.4]);
    set(h1.Parent,'XTickLabel',{'0.1','','0.2','','0.3','','0.4'});
    xlim([0.1 0.4])
    grid on
    ht = text(0.11,0.95,'(b)');
    set(ht,'FontSize',30);
    set(ht,'FontWeight','bold');
    ylabel('$\Delta T(\omega,r)$','interpreter','latex')
    
    subplot(1,4,1);
    plot(f0,depth,'Color',0.2*[1 1 1]);
    set(gca,'YDir','Reverse');
    xlabel('Resonant Frequency (Hz)');
    ylabel('Depth (m)');
    hold on;
    plot(f0_analyticaleff, depth,'Color',0.6*[1 1 1],'LineStyle','--');
    plot(f0_analytical, depth,'Color',0.6*[1 1 1],'LineStyle',':');
    set(h2.Parent,'YTick',[100 150 200 250 300 350 400 450 500]);
    set(h2.Parent,'YTickLabel',{'100','','200','','300','','400','','500'});
    grid on
    xlim([0.1 0.4])
    ht = text(0.12,120,'(a)');
    set(ht,'FontSize',30);
    set(ht,'FontWeight','bold');
        
end

%% ka as a function of depth %%
L = linspace(100,500);
a = 100;
%Leff = L + 8*a/(3*pi);
Leff = L;
ka = pi*a./(2*Leff);
subplot(1,4,4);
h3 = plot(ka,L,'k');
set(gca,'YDir','Reverse');
ylabel('Depth (m)');
xlabel('$ka$','interpreter','latex');
xlim([0.25 1]);
set(h3.Parent,'XTick',[0.25 0.5 0.75 1]);
set(h3.Parent,'YTick',[100 150 200 250 300 350 400 450 500]);
set(h3.Parent,'YTickLabel',{'100','','200','','300','','400','','500'});

