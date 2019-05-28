function [f0, Q] = resPeakPropsInversion(A, type)
% [f0, Q] = resPeakProps(A)
%
% Computes the resonant frequency, f0, and quality factor, Q, of the
% dominant peak specified in A where A is the output from resonance1d

if strcmp(type,'sim')

    P = abs(A.P); % excess pressure transfer function
    f = A.f; % frequency vector
    
    [ymax,imax,ymin,imin] = extrema(abs(P));
    [imin,min_sort_index] = sort(imin);
    ymin = ymin(min_sort_index);
    [imax,max_sort_index] = sort(imax);
    ymax = ymax(max_sort_index);
    
    flow = interp1(abs(P(1:imax(1))),A.f(1:imax(1)),ymax(1)/sqrt(2));
    fhigh = interp1(abs(P(imax(1):imin(2))),A.f(imax(1):imin(2)),ymax(1)/sqrt(2));
    if isnan(fhigh)
        fhigh = A.f(imax(1)) + (A.f(imax(1))-flow);
    end
    
    f0 = (flow+fhigh)/2;
    bw = fhigh-flow;
    Q = f0/bw;

elseif strcmp(type,'data')
    
    
    [Y,I] = max(abs(A.P));
    flow = interp1(abs(A.P(1:I)),A.f(1:I),Y*sqrt(2)/2);
    fhigh = interp1(abs(A.P(I:end)),A.f(I:end),Y*sqrt(2)/2);
    %fhigh(i) = interp1(abs(A.P(I:2*I)),A.f(I:2*I),Y*sqrt(2)/2);
    flow99 = interp1(abs(A.P(1:I)),A.f(1:I),Y*.99);
    %fhigh99(i) = interp1(abs(A.P(I:2*I)),A.f(I:2*I),Y*.99);
    fhigh99 = interp1(abs(A.P(I:end)),A.f(I:end),Y*.99);
    f0 = (flow99+fhigh99)/2;
    %f0 = (flow+fhigh)/2;
    bw = fhigh-flow;
    Q = f0/(bw);
    
end

