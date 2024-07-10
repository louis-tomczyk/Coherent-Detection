function N_symb_mixed = number_mixed_symbols(Dacc,lambda,baudrate)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : NUMBER_MIXED_SYMBOLS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2023-01-10
%   Version         : 1.0.1
%
% ----- Main idea -----
%   Estimate the minimum number of mixed symbols due to the symbol rate
%   in Optical Telecommunication Systems by taking into account also
%   the fibre dispersion.
%
% ----- INPUTS -----
%   DACC    (scalar)[ps/nm] The ACCumulated Dispersion in a fibre only
%                           due to chromatic dispersion (Polarisation
%                           Mode Dispersion not considered)
%
%   LAMBDA  (scalar)[nm]    The wavelength considered
%   BAUDRATE(scalar)[Gbd]   The symbol rate
%
% ----- OUTPUTS -----
%   N_SYMB_MIXED    (scalar)The number of symbols mixed
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    Dacc    = Dacc*1e-3;                    % [s/m]
    c       = 299792458;                    % [m/s]
    Rs      = baudrate*1e9;                 % [bd==1/s]

    delta_t         = Dacc*lambda.^2*Rs/c;  % [s]
    N_symb_mixed    = delta_t*Rs;
end