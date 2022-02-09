function out = inversion_func(peak_freq,depth)

% load peak frequency of data
% t = time vector
% f0smooth = peak frequency for the five different stations
load etna2021_data_peak_freq.mat;

N = length(t);
DEPTH_F0 = zeros(5,N);

%%%%%%%%%%%%%%%%%%%%%%%%
%%% invert for depth %%%
%%%%%%%%%%%%%%%%%%%%%%%%

% iterate over stations
figure(2); clf;
for j = 1:5
    j
    
    % iterate over time
    for k = 1:length(t)
        
        fval = f0smooth(j,k); % peak frequency value from data at station j at time k
        if isnan(fval)
            DEPTH_F0(j,k) = NaN;
        else
            tmp_f0 = abs(peak_freq - fval); % subtract freq value of interest from simulations
            [~,~,tmpmin_f0,imin_f0] = extrema(tmp_f0); % find local minimun
            
            if length(imin_f0) == 1
                DEPTH_F0(j,k) = depth(imin_f0);
            else
                DEPTH_F0(j,k) = depth(imin_f0(1));
            end
        end
        
        if DEPTH_F0(j,k) == max(depth)
            DEPTH_F0(j,k) = NaN;
        end
        
    end
    
end

out.t = t;
out.D = DEPTH_F0;