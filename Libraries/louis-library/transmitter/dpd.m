function [Epd,length_pd] = dpd(Ein,ft,D)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : DPD - Digital Pre Dispersion
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-01-12
%   Version         : 1.1.2
%
% ----- Main idea -----
%   Pre distorsion of the waveform by applying a total chromatic dispersion 
%   D. Done by propagating only in linear regime in a fake fibre.
%
% ----- INPUTS -----
%   EIN:    (structure) Fields to be compensated
%               - LAMBDA    [nm]            Wavelength
%               - FIELD     [sqrt(mW)]      Normalised electric fields
%   FT:     (structure) fiber parameters
%               - ALPHADB   [dB/km]         Power attenuation
%               - DISP      [ps/nm/km]      Dispersion
%               - SLOPE     [ps/nm²/km]     Slope of the dispersion
%               - LENGTH    [m]             Length
%               - n2        [W²/m]          Non linear index
%               - pmdpar    [ps/sqrt(km)]   Polarization Mode Dispersion
%   D:      (scalar)[ps/nm]                 Accumulated chromatic dispersion
%
% ----- OUTPUTS -----
%  EPD:     [structure] Chromatic Dispersion Compensated field
%           structure containing
%               - LAMBDA    [nm]            Wavelength
%               - FIELD     [sqrt(mW)]      Normalised electric fields
%
% ----- BIBLIOGRAPHY -----
%   Functions           : FIBER
%   Author              : Paolo SERENA
%   Author contact      : serena@tlc.unipr.it
%   Date                : 2021
%   Title of program    : Optilux
%   Code version        : 2021
%   Type                : Optical simulator toolbox - source code
%   Web Address         : https://optilux.sourceforge.io/
% ---------------------------------------------

    ft_pd           = ft;
    ft_pd.length    = 1e3;              % [m]
    ft_pd.disp      = D;                % [ps/nm/km]
    ft_pd.n2        = 0;                % [m²/W]
    ft_pd.alphadB   = 0;                % [dB/km]
    ft_pd.alphaLin  = 0;                % [dB/km]
    ft_pd.slope     = 0;
    ft_pd.beta2     = 0;
    ft_pd.beta3     = 0;
    ft_pd.gf        = 0;
    ft_pd.coupling  = "none";
    ft_pd.ismanakov = false;
    
    Epd             = fiber(Ein,ft_pd);
    length_pd       = D/ft.disp*1e3;    % [m]
end