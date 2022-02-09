%% Forward Disp Output %%
%
% Quality control the peak frequency outputs from the CRes simulations for
% three different temperature profiles: (1) 200 C, (2) 400 C, (3) 400 C
% varying linearly to 200 C.

% indexes of jagged points to remove from peak frequency
R200C = [6 7 8 9 10 11 12 20 21 22 23 24 33 34 35 44 45 54 63 226 227 242 248 250 251]';
R400C = [11 12 22 218 219 220 235 236 237 238 239 240 242 243 244 245 246 247 248 249 250 251]';
R400Cto200C = [1 2 3 13 14 15 16 26 27 28 38 48 224 225 226 242 248 250 251]';

%%%%%%%%%
% 200 C %
%%%%%%%%%
load('CRES_Output_200C.mat');
F0_init = f0_sim(:,1);
D_init = DEPTH_VECTOR;

figure(1); clf;
subplot(1,3,1);
plot(F0_init,D_init);
set(gca,'YDir','Reverse');
xlabel('Freq (Hz)');
ylabel('Depth (m)');
title('200 C');
hold on; grid on;

F0 = F0_init;
F0(R200C) = [];
D = D_init;
D(R200C) = [];
plot(F0,D);

Dintp = D_init;
F0intp = pchip(D,F0,Dintp);
plot(F0intp,Dintp);

% save peak frequency
DEPTH = Dintp;
F0_200C = F0intp;

%%%%%%%%%
% 400 C %
%%%%%%%%%
load('CRES_Output_400C.mat');
F0_init = f0_sim(:,1);
D_init = DEPTH_VECTOR;

figure(1); 
subplot(1,3,2);
plot(F0_init,D_init);
set(gca,'YDir','Reverse');
xlabel('Freq (Hz)');
ylabel('Depth (m)');
title('400 C');
hold on; grid on;

F0 = F0_init;
F0(R400C) = [];
D = D_init;
D(R400C) = [];
plot(F0,D);

Dintp = D_init;
F0intp = pchip(D,F0,Dintp);
plot(F0intp,Dintp);

% save peak frequency
DEPTH = Dintp;
F0_400C = F0intp;


%%%%%%%%%%%%%%%%%%
% 400 C to 200 C %
%%%%%%%%%%%%%%%%%%
load('CRES_Output_400Cto200C.mat');
F0_init = f0_sim(:,1);
D_init = DEPTH_VECTOR;

figure(1); 
subplot(1,3,3);
plot(F0_init,D_init);
set(gca,'YDir','Reverse');
xlabel('Freq (Hz)');
ylabel('Depth (m)');
title('400 C to 200 C');
hold on; grid on;

F0 = F0_init;
F0(R400Cto200C) = [];
D = D_init;
D(R400Cto200C) = [];
plot(F0,D);

Dintp = D_init;
F0intp = pchip(D,F0,Dintp);
plot(F0intp,Dintp);

% save peak frequency
DEPTH = Dintp;
F0_400Cto200C = F0intp;

%% Compare peak frequencies for different temperature profiles %%

figure(2); clf;
plot(F0_200C,DEPTH);
hold on;
plot(F0_400C,DEPTH);
plot(F0_400Cto200C,DEPTH);
set(gca,'YDir','Reverse');
xlabel('Freq (Hz)');
ylabel('Depth (m)');
hold on;
legend('200 C','400 C','400 C to 200 C','Location','SouthEast');
grid on;

% save formatted forward outputs
save('CRES_OUTPUT_ALL.mat','DEPTH','F0_200C','F0_400C','F0_400Cto200C');