%% FIG S3 %%
%
% Generate figure showing inversion results for different crater
% temperatures for before the onset of the flank eruption. Inversion 
% results have been previously calculated and are saved in invOutput/.
%
% Written by Leighton Watson
% March 11, 2020
% leightonwatson@stanford.edu // leightonmwatson@gmail.com

clear all; clc;
cmap = get(gca,'ColorOrder');
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',18);

path(pathdef)
addpath data/
addpath invOutput/
addpath ../source/resonance/
addpath ../source/SBPoperators/
addpath ../source/inv/

% clear figures
figHand1 = figure(1); clf;
set(figHand1,'Position',[100 100 1200 600]);
figHand2 = figure(2); clf;
set(figHand2,'Position',[100 100 1200 400]);

%%% load data %%%
dataStr = {'Etna2018Phase1.mat'};

%%% crater temperature %%%
outputStr = {'InvOut_Etna2018Phase1_Nit100000_R0100m_D200m_T100C.mat',...
    'InvOut_Etna2018Phase1_Nit100000_R0100m_D200m_T200C.mat',...
    'InvOut_Etna2018Phase1_Nit100000_R0100m_D200m_T300C.mat'};

% iterate over all output files

std_flag = 1; % plot standard deviation of geometry

for i = 1:length(outputStr)
    
    %%% load data %%%
    load(outputStr{i});
    
    %% crater geometry %%
    
    figure(1); subplot(1,2,1);
    burnIn = ceil(count/10); % remove the first 10% of successful samples to remove effect of initial conditions
    xGeom = x(:,1:length(geomParams)); % extract output geometry parameters
    xBurnIn = xGeom(burnIn:end,:); % remove burn-in values
    xMean = mean(xBurnIn); % average parameters
    geomParamsMmean = xMean.*geomParams; % convert to physical values
    shapeF = geomFunction(geomParamsMmean);
    depthF = shapeF(1,1);
    h = plot(shapeF(:,2), shapeF(:,1),'Color',cmap(i,:)); hold on;
    h.Parent.YTick = [0 25 50 75 100 125 150 175 200];
    plot(-shapeF(:,2), shapeF(:,1),'Color',cmap(i,:));
    plot([-shapeF(1,2) shapeF(1,2)],[depthF depthF],'Color',cmap(i,:));
    set(gca,'YDir','Reverse'); 
    xlabel('Radius (m)'); ylabel('Depth (m)');
    xlim([-150 150])
    
    avg_radius = mean(shapeF(:,2));
    crater_vol(i) = depthF * pi * avg_radius^2;
    
    if std_flag == 1 % plot standard deviation of geometry
    
        %standard deviation of geometry (positive)
        xStd = std(xBurnIn); % compute standard deviation once burn in has been removed
        geomStd1 = (xMean + xStd) .* geomParams; % convert to physical values
        shapeStd1 = geomFunction(geomStd1);
        depthStd1 = shapeStd1(1,1);
        
        plot(shapeStd1(:,2), shapeStd1(:,1),...
            'Color', cmap(i,:),'LineStyle',':'); hold on;
        plot(-shapeStd1(:,2), shapeStd1(:,1),...
            'Color', cmap(i,:),'LineStyle',':');
        plot([-shapeStd1(1,2) shapeStd1(1,2)],[depthStd1 depthStd1],...
            'Color',cmap(i,:),'LineStyle',':');
        
        %standard deviation of geometry (negative)
        geomStd2 = (xMean - xStd) .* geomParams; % convert to physical values
        shapeStd2 = geomFunction(geomStd2);
        depthStd2 = shapeStd2(1,1);
        
        plot(shapeStd2(:,2), shapeStd2(:,1),...
            'Color', cmap(i,:),'LineStyle',':'); hold on;
        plot(-shapeStd2(:,2), shapeStd2(:,1),...
            'Color', cmap(i,:),'LineStyle',':');
        plot([-shapeStd2(1,2) shapeStd2(1,2)],[depthStd2 depthStd2],...
            'Color',cmap(i,:),'LineStyle',':');
        
    end
    
    %% spectra %%
    
    figure(1); subplot(1,2,2);
    if i == 1
        for j = 1:length(dataStr)
            load(dataStr{j});
            h = plot(dataF, abs(dataAmp)./max(abs(dataAmp)),'Color','k');
            h.Color(4) = 0.5;
            h.Parent.YTick = [];
            hold on; xlim([0 3]);
            xlabel('Frequency (Hz)'); 
        end
        
    end
    
    % source
    [S, ~, ~, ~] = sourceFunction(1, srcParams, srcStyle, discrParams); % compute source in frequency domain
    % resonance 1d properties
    style = 'baffled piston'; % acoustic radiation model ('monopole' or ' baffled piston')
    order = 4; % order of numerical scheme (4, 6 or 8)
    N = discrParams(2); % number of grid points (in time)
    Nf = discrParams(3); % number of frequency samples
    Nyquist = discrParams(4); % Nyquist frequency
    dt = discrParams(5); % time step
    freq = [0 Nyquist]; % frequency vector
    % resonance 1d
    res = resonance1d(shapeF, depthF, freq, Nf, style, order, M); % compute transfer function
    sim.f = res.f; % frequency vector
    sim.P = (res.P(1:N/2+1).*S(1:N/2+1))./max(res.P(1:N/2+1).*S(1:N/2+1)); % convolve transfer function and source function and normalize amplitudes
    % convert to time domain
    sim.pFreq = [sim.P conj(sim.P(end-1:-1:2))];
    sim.pFreq(N/2+1) = real(sim.pFreq(N/2+1));
    sim.pTime = ifft(sim.pFreq,'symmetric')/dt;
    % apply filtering
    filterband = [filterProps(1) filterProps(2)]; % frequency band to filter
    filterorder = filterProps(3); % order of butterworth filter
    Fs = filterProps(4);
    sim.pTimeFilt = bandpass_butterworth(sim.pTime,filterband,Fs,filterorder);
    % convert back to frequency domain
    L = length(sim.pTimeFilt);
    NFFT = 2^nextpow2(L);
    sim.pFreqFilt = fft(sim.pTimeFilt,NFFT)/L;
    sim.pFreqAbsNorm = abs(sim.pFreqFilt(1:NFFT/2+1))./max(abs(sim.pFreqFilt(1:NFFT/2+1)));
    sim.fFilt = Fs/2*linspace(0,1,NFFT/2+1);
    % plot
    plot(sim.fFilt(1:N/2+1), abs(sim.pFreqAbsNorm(1:N/2+1)),'Color',cmap(i,:),'LineStyle',':');
    
    %% parameter histograms %%
    
    depths = [0 50 100 150];
    param_hist = zeros(length(xBurnIn),length(depths));
    
    tic
    figure(2);
    for j = 1:length(xBurnIn)
        disp(strcat('Percentage complete:',num2str(j/length(xBurnIn)*100),'%'));
        geom_trunc = xBurnIn(j,:) .* geomParams;
        shape_intp = geomFunction(geom_trunc);
        for m = 1:length(depths)
            if shape_intp(1,1) < depths(m)
                param_hist(j,m) = NaN;
            else
                param_hist(j,m) = shape_intp(end-depths(m),2);
            end
        end
    end
    
        
    for m = 1:length(depths)
        
        if m == 1 % plot depth
            subplot(1, length(depths)+1,1);
            h = histogram(xBurnIn(:,1).*geomParams(1),'Normalization','probability');
            hold on; box on;
            h.Parent.YTick = {};
            h.BinWidth = 1;
            h.FaceColor = cmap(i,:);
            h.EdgeColor = cmap(i,:);
            h.EdgeAlpha = 0;
            xlim([50 250])
        end
        
        subplot(1,length(depths)+1,m+1);
        h = histogram(param_hist(:,m),'Normalization','probability');
        hold on; box on;
        h.Parent.YTick = {};
        h.BinWidth = 1;
        h.FaceColor = cmap(i,:);
        h.EdgeColor = cmap(i,:);
        h.EdgeAlpha = 0;
        
    end
    toc

end

