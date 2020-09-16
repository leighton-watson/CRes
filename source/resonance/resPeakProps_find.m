function [f0, Q] = resPeakProps_find(f,G)
% computes resonant frequency and quality factor of fundamental resonant mode.
%
% INPUT
% f = frequency vector
% G = amplitude spectra
%
% OUTPUT
% f0 = resonant frequency
% Q = quality factor
[Y,I] = max(abs(G)); 
Is = I;
Ys = Y;
flow_idx = find(G>Y*sqrt(2)/2,1,'first');
flow = f(flow_idx);
fhigh_idx = find(G>Y*sqrt(2)/2,1,'last');
fhigh = f(fhigh_idx);

flow99_idx = find(G>Y*.99,1,'first');
flow99 = f(flow99_idx);

fhigh99_idx = find(G>Y*.99,1,'last');
fhigh99 = f(fhigh99_idx);
%flow = interp1(abs(G(1:I)),f(1:I),Y*sqrt(2)/2);
%fhigh = interp1(abs(G(I:2*I)),f(I:2*I),Y*sqrt(2)/2);
%flow99 = interp1(abs(G(1:I)),f(1:I),Y*.99);
%fhigh99 = interp1(abs(G(I:2*I)),f(I:2*I),Y*.99);

fpeak = (flow99+fhigh99)/2;
f0 = (flow+fhigh)/2; % resonant frequency
bw = fhigh-flow; % bandwidth
Q = f0/(bw); % quality factor