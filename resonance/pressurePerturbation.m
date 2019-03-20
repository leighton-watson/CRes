function P = pressurePerturbation(input, style, M)
% P = pressurePerturbation(input)
%
% input = output structure from resonance1d
% style = 'baffled piston' or 'monopole' description of sound radiation

    r = M.r; % distance from vent to receiver
    rho = M.rhoA; % density of air
    shape = input.geometry; % geometry [depth, radius]

    f = input.f; % frequency
    omega = 2*pi*f; % angular frequency

    if strcmp(style,'baffled piston')

        %%% BAFFLED PISTON %%%
        % Compute the pressure perturbation at a distance r from the vent assuming
        % a circular piston embedded in an infinite plane baffle (equation 7.30 from
        % Rossing and Fletcher, 2004)

        c = M.cA; % speed of sound

        a = shape(end,2); % radius of vent

        u = input.vOutlet; % outlet velocity transfer function
        k = 2*pi*f/c; % wave number
        theta = pi/2; % angle of vector from source to receiver. Measured from vertical

        p = 1/(2)*1i.*omega*rho.*u.*a.^2 .* ((exp(-1i*k*r))/r) .* ...
            ((2*besselj(1,k*a*sin(theta)))./(k*a*sin(theta))); % excess pressure
        P = p;

    elseif strcmp(style,'monopole')

        %%% MONOPOLE %%%
        % Monopole source description (after Johnson and Miller, 2014)

        ventRadius = shape(end,2);
        ventArea = pi*ventRadius^2;

        vOutlet = input.vOutlet; % outlet velocity transfer function
        vDot = vOutlet.*1i.*omega; % rate of change of outlet velocity transfer function
        deltaP = rho*vDot*ventArea / (2*pi*r); % pressure perturbation
        P = deltaP;
        
    else
        error('Incorrect style of sound radiation. Please choose from --baffled piston-- or --monopole--');        

    end
end






