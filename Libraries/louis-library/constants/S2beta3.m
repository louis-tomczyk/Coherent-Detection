function beta3 = S2beta3(S,D,lambda)

% Inputs
%   S       [ps/nm²/km]
%   D       [ps/nm/km]
%   LAMBDA  [nm]
% Output
%   BETA3   [s³/rad²/m]

    S       = S*1e3;                % [s/m²/m]
    D       = D*1e-6;               % [s/m/m]
    lambda  = lambda *1e-9;         % [m]
    c       = 299792458;            % [m/s]
    beta3   = (lambda/2/pi/c).^2.*(2*lambda*D+lambda.^2*S); % [s³/rad²/m]
  
end