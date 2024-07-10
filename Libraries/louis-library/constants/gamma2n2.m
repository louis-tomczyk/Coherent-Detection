function n2 = gamma2n2(gamma_f,lambda,Aeff)

    % Inputs
    %   GAMMA_F [1/W/km]
    %   LAMBDA  [nm]
    %   AEFF    [um^2]
    % Output
    %   n2      [m^2/W]

    gamma_f = gamma_f*1e-3; % [1/W/m]
    lambda  = lambda*1e-9;  % [m]
    Aeff    = Aeff*1e-12;   % [m^2]

    n2      = gamma_f*lambda*Aeff/(2*pi);

end