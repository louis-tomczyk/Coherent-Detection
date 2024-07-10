function gamma_f = n22gamma(n2,lambda,aeff)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : n22gamma
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-05-23
%   Version         : 1.0
%
% ----- Main idea -----
%   Convert the non linear refractive index N2 into the non linear
%   parameter GAMMA
%
% ----- INPUTS -----
%   N2      [m^2/W] Non linear refractive index
%   LAMBDA  [nm]    Wavelength
%   AEFF    [um^2]  Effective Fiber Area
%
% ----- OUTPUTS -----
%   GAMMA_F [1/W/m] Non linear parameter
%       the '_f' is to avoid any confusion with the native function "gamma"
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------
    
    lambda  = lambda*1e-9;              % [m]
    aeff    = aeff*1e-12;               % [m²]
    gamma_f = 2*pi*n2/(lambda*aeff);    % [1/W/m]
end