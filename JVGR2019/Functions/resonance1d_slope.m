function output = resonance1d_slope(geometry, depth, freq, Nf, style, order, M, slope_angle)
% output = resonance1d(geometry, depth, freq, Nf, style, order, M)
%
% evaluates the resonant frequencies of a crater with arbitrary
% axisymmetic geometry
%
% INPUTS
% geometry = structure containing the depth discretization (geometry.x) and 
% radius of the crater at those depths (geometry.radius)
% depth = depth of lava lake surface in crater
% freq = maximum and minimum frequency of interest
% Nf = number of frequency samples
% style = sound radiation description ('baffled piston' or 'monopole')
% order = internal order of accuracy
% M = model parameters
%
% OUTPUTS
% save outputs into the structure output
% output.geometry = crater geometry in the same format as geometry
% output.depth = depth of lava lake surface
% output.f = frequency vector
% output.T = transfer function
% output.pOutlet = outlet pressure transfer function
% output.vOutlet = outlet velocity transfer function
% output.P = far-field pressure perturbation transfer function


    % truncate shape file to the desired depth
    z = geometry(:,1);
    iz = find(z==depth);
    z = z(iz:end);

    % set default options
    if (nargin <= 5 || isempty(order)); order = 8; end;
    if (nargin <= 4 || isempty(style)); style = 'baffled piston'; end
    if (nargin <= 3 || isempty(Nf)); Nf = 100; end; 
    
    % geometry
    radius = geometry(iz:end,2); 
    S = pi*radius.^2;
   
    %dx = 1;
    dx = abs(z(2)-z(1));
    Nx = length(z)-1;
    [D,Hinv,~,~,~,~,~] = SBPoperatorsRussian(Nx,dx,order);
        
    % problem parameters
    % M = problemParameters();
    % modelParameters
    
    % SAT penalties
    SAT = [M.cC(1) M.cC(end)]*Hinv; % penalty relaxation rate
    
    
    % Compute the transfer function between the inlet velocity and the
    % properties (velocity and pressure) at the outlet

    % frequencies
    f = linspace(freq(1),freq(end),Nf);

    % initialize matrices
    m = Nx+1;
    I = speye(size(D)); % identity matrix
    b = zeros(2*m,1); % vector of forcing terms
    q = zeros(2*m,length(f)); % solution matrix
    D1 = D; D2 = D;
    
    %for i = 1:length(z)
    for i = 1:Nx
        D1(i,:) = D(i,:) * S(i);
        D2(i,:) = D(i,:) / S(i);
    end

    %Z0 = M.rhoC.*M.cC./S(1); %characteristc impedance of pipe
    %ZN = M.rhoC.*M.cC./S(end);
    
    Z0 = M.rhoC(1).*M.cC(1)./S(1); %characteristc impedance of pipe at base
    ZN = M.rhoC(end).*M.cC(end)./S(end); %characteristc impedance of pipe at outlet

    for i = 1:length(f) %(length(f)=Nf)

        omega = 2*pi*f(i); % angular frequency

        % impedance
        Zoutlet = flanged_opening(f(i),M.rhoA,M.cA,radius(end)); %outlet impedance

        % momentum equation
        A1 = 1i.*omega.*I;
        A2 = 1./M.rhoC.*D1;
        A1(1,1) = A1(1,1) + SAT(1);
        A1(end,end) = A1(end,end) + SAT(end) - SAT(end) .* (ZN./(Zoutlet+ZN));
        A2(end,end) = A2(end,end) - SAT(end) ./ (Zoutlet+ZN);

        % mass equation
        B1 = M.pC.*D2;
        B2 = 1i.*omega.*I;
        B1(1,1) = B1(1,1) + SAT(1) .* Z0;
        B1(end,end) = B1(end,end) - SAT(end) .* (Zoutlet*ZN)./(Zoutlet+ZN);
        B2(end,end) = B2(end,end) + SAT(end) - SAT(end) .* (Zoutlet./(Zoutlet+ZN));

        % form system of governing equations
        A = [A1 A2;B1 B2];
 

        % boundary conditions
        b(1) = SAT(1).*S(1); % inlet velocity
        b(Nx+2) = SAT(1).*Z0.*S(1); % inlet pressure

        % solve system of equations
        q(:,i) = A\b;

    end

    presOutTransfer = q(end,:); %outlet pressure transfer function
    velOutTransfer = q(Nx+1,:)./S(end); %outlet velocity transfer function

    %%% SAVE OUTPUTS %%%

    output.geometry = [z, radius]; % crater geometry
    output.depth = depth; % depth of lava lake surface
    
    output.f = f; % frequency
    output.T = q; % transfer function

    output.pOutlet = presOutTransfer; % outlet pressure transfer function
    output.vOutlet = velOutTransfer; % outlet velocity transfer function
   
    % what are these used for?
    %output.p = (presOutTransfer); % outlet pressure transfer function
    %output.v = (velOutTransfer); % outlet velocity transfer function
    
    output.P = pressurePerturbation_slope(output, style, M, slope_angle); % far-field pressure perturbation transfer function
    if output.f(1) == 0
        output.P(1) = 0; % set DC component equal to zero
    end
end
