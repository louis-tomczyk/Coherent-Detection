function ft = set_ft(las,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : SET_FT
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-03-15
%   Version         : 1.5.2
%
% ----- Main idea -----
%   Set the FT structure from the kind of fibre, the laser and axis 
%   parameters
%   
%   !!! CAUTION !!! Optilux versions 2009/2022, in 2009 
%                   ISMANAKOV & COUPLING non existing fields
% 
% ----- INPUTS -----
%   FIBRE_TYPE  (string) The kind of the fibre
%   LAS:        (structure) containing the Laser parameters with at least
%                   - LAM       [nm]    Wavelengths
%                   - N_POLAR   []      Number of polarisations
%                   - PLIN      [mW]    Output power of the laser
%   AXIS:       (structure) containing the Axis parameters
%                   - NSYMB     []      Number of symbols to send
%                   - NT        []      Number of bits per symbols
%                   - SYMBRATE  [GBAUD] Symbol rate
%
% ----- OUTPUTS -----
%   D_PROP  (array) containing:
%               - BETA2 [s²/rad²/m] Second order of chromatic dispersion
%               - BETA3 [s³/rad³/m] Third order of chromatic dispersion
%               - LD    [m]         Typical Length Dispersion
%                           ( Amplitude(z=LD) = Amplitude(z=0)/sqrt(2) )
%   FT:   (structure) containing the Fibre parameters
%               - ALPHAdB   [dB/km]         Attenuation - log scale
%               - ALPHALIN  [1/m]           Attenuation - lin scale
%               - LENGTH    [m]             Fibre span length
%               - DISP      [ps/nm/km]      Dispersion
%               - SLOPE     [ps/nm²/km]     Slope
%               - AEFF      [µm²]           Effective area
%               - GF        [1/W/km]        Non linear param
%               - N2        [m²/W]          Non linear index
%               - PMDPAR    [ps/sqrt(km)]   polarisation mode dispersion
%               - NPLATES   []              Number of plates to simulate PMD
%               - ISMANAKOV []              Boolean to decide to Solve or not
%                                           Manakov equation
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    if nargin == 1
        ft = struct();
    else
        ft = varargin{1};
    end
    
    if isfield(ft,'type') == 0
        ft.type = 'SMF';
    end

    ft.lambda   = las.lam;
    fibre_type  = convertStringsToChars(ft.type);

    available_types = string(["B2B","noDISP","noNL","noATT","SMF","DSF","NZDSF","LEAF"]);
    if size(find(strcmp(available_types,ft.type)),2) == 0
        fibre_type = input("Fibre type not available.\nPlease choose among 'SMF / DSF / NZDSF / LEAF'\n>>",'s');
        ft         = set_ft(fibre_type,las);
    end

    switch fibre_type
        case "B2B"
            ft.alphadB  = 0;        % [dB/km]       attenuation
            ft.disp     = 0;        % [ps/nm/km]    dispersion
            ft.slope    = 0;        % [ps/nm²/km]   slope

            ft.aeff     = 80;                               % [µm²]    effective area
            ft.gf       = 0;                                % [1/W/km] non linear param
            ft.n2       = gamma2n2(ft.gf,las.lam,ft.aeff);  % [m²/W]   nonlinear index

        case "noDISP"
            ft.alphadB  = 0.22;     % [dB/km]       attenuation
            ft.disp     = 0;        % [ps/nm/km]    dispersion
            ft.slope    = 0;        % [ps/nm²/km]   slope

            ft.aeff     = 80;                               % [µm²]    effective area
            ft.gf       = 1.31;                             % [1/W/km] non linear param
            ft.n2       = gamma2n2(ft.gf,las.lam,ft.aeff);  % [m²/W]   nonlinear index

        case "noNL"
            ft.alphadB  = 0.22;     % [dB/km]       attenuation
            ft.disp     = 16.7;     % [ps/nm/km]    dispersion
            ft.slope    = 0.058;    % [ps/nm²/km]   slope

            ft.aeff     = 80;                               % [µm²]    effective area
            ft.gf       = 0;                                % [1/W/km] non linear param
            ft.n2       = gamma2n2(ft.gf,las.lam,ft.aeff);  % [m²/W]   nonlinear index

        case "noATT"
            ft.alphadB  = 0;        % [dB/km]       attenuation
            ft.disp     = 16.7;     % [ps/nm/km]    dispersion
            ft.slope    = 0.058;    % [ps/nm²/km]   slope

            ft.aeff     = 80;                                % [µm²]    effective area
            ft.gf       = 1.31;                              % [1/W/km] non linear param
            ft.n2       = gamma2n2(ft.gf,las.lam,ft.aeff);   % [m²/W]   nonlinear index

        case "SMF"
            ft.alphadB  = 0.206;     % [dB/km]       attenuation
            ft.disp     = 17.0;     % [ps/nm/km]    dispersion
            ft.slope    = 0.058;    % [ps/nm²/km]   slope

            ft.aeff     = 80;                               % [µm²]    effective area
            ft.gf       = 1.31;                             % [1/W/km] non linear param
            ft.n2       = gamma2n2(ft.gf,las.lam,ft.aeff);  % [m²/W]   nonlinear index

        case "DSF"
            ft.alphadB  = 0.2;      % [dB/km]       attenuation
            ft.disp     = 0;        % [ps/nm/km]    dispersion
            ft.slope    = 0.08;     % [ps/nm²/km]   slope

            ft.aeff     = 50;                               % [µm²]    effective area
            ft.gf       = 2.1;                              % [1/W/km] non linear param
            ft.n2       = gamma2n2(ft.gf,las.lam,ft.aeff);  % [m²/W]   nonlinear index

        case "NZDSF"
            ft.alphadB  = 0.2;      % [dB/km]       attenuation
            ft.disp     = 4.5;      % [ps/nm/km]    dispersion
            ft.slope    = 0.045;    % [ps/nm²/km]   slope

            ft.aeff     = 50;                               % [µm²]    effective area
            ft.gf       = 2.1;                              % [1/W/km] non linear param
            ft.n2       = gamma2n2(ft.gf,las.lam,ft.aeff);  % [m²/W]   nonlinear index

        case "LEAF"
            ft.alphadB  = 0.2;      % [dB/km]       attenuation
            ft.disp     = 4;        % [ps/nm/km]    dispersion
            ft.slope    = 0.09;     % [ps/nm²/km]   slope

            ft.aeff     = 72;                               % [µm²]    effective area
            ft.gf       = 1.5;                              % [1/W/km] non linear param
            ft.n2       = gamma2n2(ft.gf,las.lam,ft.aeff);  % [m²/W]   nonlinear index
    end

    ft.beta2 = D2beta2(ft.disp,ft.lambda);
    ft.beta3 = S2beta3(ft.slope,ft.disp,ft.lambda);
    
    %--- linear effects
    % power
    ft.alphaLin = log(10)/10*ft.alphadB*1e-3;   % attenuation [1/m]
    if isfield(ft,'length') == 0
        ft.length   = 100e3;                        % [m]
    end

    % polarisation
    if las.n_polar == 1
        ft.coupling = 'none';
        ft.ismanakov= false;    % Solve Manakov equation
        ft.nplates  = 1;
        ft.pmdpar   = 0;        % [ps/sqrt(km)] polarisation mode dispersion  --- 0.04 (SMF)
    else
        if isfield(ft,'pmdpar') == 0
            ft.ismanakov    = true;
            ft.nplates      = 1;
            ft.pmdpar       = 0;
        elseif isfield (ft,"pmdpar") == 1 && ft.pmdpar == 0
            ft.ismanakov    = false;
            ft.nplates      = 100;
            ft.pmdpar       = 0.04;
        end

        if isfield(ft,'coupling') == 0
            ft.coupling     = 'pol';
            ft.ismanakov    = false;
            ft.nplates      = 100;
            ft.pmdpar       = 0.04;
        elseif isfield (ft,"coupling") == 1 && strcmp(ft.coupling,'none') == 1
            ft.ismanakov    = true;
            ft.nplates      = 1;
            ft.pmdpar       = 0;
        end
    
    if strcmp(fibre_type,'B2B') == 1
        ft.coupling     = 'none';
        ft.ismanakov    = false;
        ft.nplates      = 1;
        ft.pmdpar       = 0;
    end

    % chromatic dispersion
    beta2       = D2beta2(ft.disp,ft.lambda);          % [s²/rad²/m]
    beta3       = S2beta3(ft.slope,ft.disp,ft.lambda); % [s³/rad³/m]

    % values are extremely low (~ 1e-26/1e-40, might return '0' in shell,
    % but value stored are NOT null)
    ft.beta2    = beta2;
    ft.beta3    = beta3;
    
    if nargin == 3
        Axis    = varargin{2};
        ft.LD   = Axis.dt^2/abs(beta2)*1e-18;   % [m]     characteristic dispersion length
        ft.LNL  = 1/(ft.gf*las.Plin*1e-3)*1e3;  % [m]     characteristic non linearity length
    end
    
    ft          = sort_struct_alphabet(ft);

end