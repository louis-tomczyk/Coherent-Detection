function las = set_las(tx,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : SET_LAS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-03-15
%   Version         : 1.4.4
%
% ----- Main idea -----
%   Set the LAS structure from the laser transmitter parameters
% 
% ----- INPUTS -----
%   LAS         (structure) containing the LASer parameters
%                   - N_POLAR   (scalar)[]  Number of polarisations
%                   - LAM       (scalar)[nm]Laser wavelength
%                   - TYPE      (string)    Type of laser
%                   Can be:
%                       - HPFL: (High Power Fiber Laser) for high noise
%                               level laser
%                       - LRL:  (Low RIN Laser) for low noise level laser
%                       - IDEAL:no noise
%   TX          (structure) containing the Transmitter parameters
%                   - PLIN     ()[mW] Power in signal:
%                                   (scalar) if 1 polarisation
%                                   (array)  if 2 polarisations
%                   - NCH      (scalar)[] Number of channels
%
% ----- OUTPUTS -----
%   FT    (structure) containing the Fibre parameters
%               - LAM       (scalar)[nm]        Laser wavelength
%               - LINEWIDTH (scalar)[GHz]       Laser linewidth
%               - N_POLAR   (scalar)[]          Number of polarisations
%               - N0        (scalar)[dB/GHz]    Noise power density
%               - PLIN      (scalar)[mW]        Power per channel
%               - PDBM      (scalar)[dBm]       Power per channel
%               - TYPE      (string)[]          Type of the laser, can be:
%                                               - HPFL (High Power Fibre
%                                                       Laser)
%                                               - LRL   (Low Rin Laser)
%                                               - IDEAL (noisefree)
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    if nargin == 1
        las.lam     = 1550;
        las.n_polar = 1;
        las.type    = 'LRL';
    else
        las = varargin{1};
    end

    if isfield(las,'lam') == 0
        las.lam     = 1550; % [nm]      carrier wavelength
    end

    if isfield(las,'n_polar') == 0
        las.n_polar = 1;
    end

    if isfield(las,'type') == 0
        las.type        = "ideal";
        las.linewidth   = 0;
        las.n0          = -inf;
    else
        laser_type  = convertStringsToChars(las.type);
        switch laser_type
            case "HPFL"
                las.linewidth   = 1e-3; % 1  [MHz]
                las.n0          = -50;
            case 'LRL'
                las.linewidth   = 1e-5; % 10  [kHz]
                las.n0          = -80;
            otherwise
                las.linewidth   = 0;
                las.n0          = -inf;
        end
    end

    las.Plin    = tx.Plin/tx.Nch;
    las.PdBm    = 10*log10(las.Plin);

    las         = sort_struct_alphabet(las);

end