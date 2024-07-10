function beta2 = D2beta2(D,lambda)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : D2BETA2
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-05-19
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Convert the Dispersion paramter of a fibre into
%   the propagation order parameter Beta2
%
% ----- INPUTS -----
%   D       [ps/nm/km]  fibre dispersion parameter
%   LAMBDA  [nm]        wavelength
%
% ----- OUTPUTS -----
%   BETA2   [s²/rad²/m] second order dispersion
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    lambda  = lambda *1e-9;         % [m]
    D       = D*1e-6;               % [s/m/m]
    c       = 299792458;            % [m/s]

    beta2   = -lambda.^2/2/pi/c*D;  % [s²/rad²/m]
  
end