%% MAKE FIG FORWARD %%
%
% Make figure showing the crater geometry and the transfer function as a
% function of depth and frequency (psuedo-spectrogram)

clear all; clc;
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',18);
cmap = get(gca,'ColorOrder');

%%%%%%%%%%%%%%%%
%%% Geometry %%%
%%%%%%%%%%%%%%%%

% CRes geometry
load EtnaGeomCRes.mat
figHand = figure(1); clf;
set(figHand,'Position',[100 100 300 600]);

plot(shape(:,2), shape(:,1),'Color',cmap(1,:));
hold on;
set(gca,'YDir','Reverse');
plot(-shape(:,2), shape(:,1),'Color',cmap(1,:));
ylim([-20 300]); grid on; hold on;
xlabel('Radius (m)');
ylabel('Depth (m)');
xlim([-100 100])
%hline(58)

% DEM geometry
W = readmatrix('w_section.txt');
% shift profile so that outlet height is zero and crater is centered at zero
z_shift = max(W(:,2))-10;
r_shift = 160;
W0(:,1) = W(:,1)-r_shift;
W0(:,2) = -(W(:,2)-z_shift);
plot(W0(:,1), W0(:,2),'k');

%%%%%%%%%%%%%%%%%
%%% SPECTOGRAM %%%
%%%%%%%%%%%%%%%%%

load CRES_Output_200C.mat
figHand = figure(2); clf;
set(figHand,'Position',[100 100 600 500]);
[FF,DD] = meshgrid(f,DEPTH_VECTOR);
h = surf(FF,DD,abs(TRANS_FUNC)');
hold on;
view(2);
shading interp
ylabel('Depth (m)');
xlabel('Frequency (Hz)');
xlim([0 5]);
ylim([min(DEPTH_VECTOR) max(DEPTH_VECTOR)])
set(gca,'YDir','Reverse');
hold on;
colorbar; caxis([0 8])


load CRES_OUTPUT_ALL.mat
plot3(F0_200C,DEPTH,10*ones(size(DEPTH)),'r:');
