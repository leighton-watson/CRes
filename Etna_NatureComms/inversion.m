%% Inversion %%
% Run inversion for the different temperature profiles. The code for the
% inversion procedure is in "inversion_func.mat".

clear all; clc;
load CRES_OUTPUT_ALL.mat;
figure(1); clf;

% 200 C
out1 = inversion_func(F0_200C,DEPTH);

for j = 1:5
    figure(1); subplot(3,1,1);
    plot(out1.t, out1.D(j,:),'.','MarkerSize',10);
    hold on;   
end
xlabel('Time (hrs) since 00:00 on Feb 19, 2021');
ylabel('Depth (m)');
set(gca,'YDir','Reverse'); grid on;
title('200 C'); ylim([40 220]);
save('inversionDepth_200C.mat','out1');
% 400 C
out2 = inversion_func(F0_400C,DEPTH);

for j = 1:5
    figure(1); subplot(3,1,2);
    plot(out2.t, out2.D(j,:),'.','MarkerSize',10);
    hold on;
end
xlabel('Time (hrs) since 00:00 on Feb 19, 2021');
ylabel('Depth (m)');
set(gca,'YDir','Reverse'); grid on;
title('400 C');
save('inversionDepth_400C.mat','out2');
% 400 C to 200 C
out3 = inversion_func(F0_400Cto200C,DEPTH);

for j = 1:5
    figure(1); subplot(3,1,3);
    plot(out3.t, out3.D(j,:),'.','MarkerSize',10);
    hold on;
end
xlabel('Time (hrs) since 00:00 on Feb 19, 2021');
ylabel('Depth (m)');
set(gca,'YDir','Reverse'); grid on;
title('400 C to 200 C');
save('inversionDepth_400Cto200C.mat','out3');

%% compare between temperature profiles for the different stations
figure(2); clf;
for i = 1:5;
    subplot(3,2,i);
    plot(out1.t,out1.D(i,:),'.','MarkerSize',10);
    hold on;
    plot(out2.t,out2.D(i,:),'.','MarkerSize',10);
    plot(out3.t,out3.D(i,:),'.','MarkerSize',10);
    xlabel('Time (hrs) since 00:00 on Feb 19, 2021');
    ylabel('Depth (m)');
    set(gca,'YDir','Reverse'); grid on;
    ylim([50 250])
end

figure(3); clf;
for i = 1:5
    subplot(3,2,i);
    plot(out1.t,out2.D(i,:)-out1.D(i,:));
    xlabel('Time (hrs) since 00:00 on Feb 19, 2021');
    ylabel('Depth (m)');
    set(gca,'YDir','Reverse'); grid on;
    M = nanmean(out2.D(i,:)-out1.D(i,:));
    hline(M);
    ylim([0 50])
end



