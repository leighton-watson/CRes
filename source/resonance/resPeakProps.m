function [f0, Q] = resPeakProps(f,G)

[Y,I] = max(abs(G)); 
Is = I;
Ys = Y;
flow = interp1(abs(G(1:I)),f(1:I),Y*sqrt(2)/2);
fhigh = interp1(abs(G(I:2*I)),f(I:2*I),Y*sqrt(2)/2);
flow99 = interp1(abs(G(1:I)),f(1:I),Y*.99);
fhigh99 = interp1(abs(G(I:2*I)),f(I:2*I),Y*.99);

fpeak = (flow99+fhigh99)/2;
f0 = (flow+fhigh)/2; % resonant frequency
bw = fhigh-flow; % bandwidth
Q = f0/(bw); % quality factor