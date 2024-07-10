function tx = set_tx(tx)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : TX - Transmitter
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-03-02
%   Version         : 1.2.2
%
% ----- Main idea -----
%   Set the transmitter paramters
% 
% ----- INPUTS -----
%   TX          (structure) containing the fields:
%               - MODOFR        (string)[]      Modulation format
%               - NCH           (scalar)[]      Number of channels (WDM)
%               - PDBM          (array)[dBm]    Powers in the channels
%               - PD            (array)[pd/nm]{optional}
%                                               Predispersion
%               - SEED          (bool)[]        Constant bit sequence
%               - EMPH          (string)[]{optional}
%                                               Emphasizer filter
%               - LAMBDA        (scalar)[nm]{optional}    
%                                               Central wavelength
%               - ROLLOFF       (scalar)[]{optional}      
%                                               Pulse roll-off
%               - PULSE_SHAPE   (string)[]{optional}      
%                                               Shape of the pulse
%               - DELTAF        (scalar)[GHz]{optional}
%                                               Frequency spacing (WDM)
%
% ----- OUTPUTS -----
%   TX          (structure) filled with (if not as input)
%               - PLIN          (array)[mW]     idem
%               - PD            (array)[pd/nm]  Predispersion
%               - SEED          (bool)[]        Constant bit sequence
%               - EMPH          (string)[]      Emphasizer filter
%               - LAMBDA        (scalar)[nm]    Central wavelength
%               - ROLLOFF       (scalar)[]      Pulse roll-off
%               - PULSE_SHAPE   (string)[]      Shape of the pulse
%               - DELTAF        (scalar)[GHz]   Frequency spacing (WDM)
%               - PAT           (array)[]       Symbol sequence to send
%               - PATBIN        (array)[]       Bits for each symbol
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    tx.emph = 'asin';   % digital-premphasis type

    if isfield(tx,'lambda') == 0
        tx.lambda   = 1550;      % [nm]
    end
   
    if isfield(tx,'rolloff') == 0
        tx.rolloff  = 0.1;      % pulse roll-off
    end

    if isfield(tx,'pulse_shape') == 0
        tx.pulse_shape  = "rootrc";
    end

    if isfield(tx,'deltaf') == 0
        tx.deltaf = 37.5; % [GHz]
    end

    tx.spac = tx.deltaf*tx.lambda^2/299792458; % [nm]
    tx.Plin = 10^(tx.PdBm/10);

    if isfield(tx,'seed') == 0
        tx.seed = 0;
    end

    if isfield(tx,'pd') == 0
        tx.pd = 0; % [ps/nm]
    end

    if isfield(tx,'Nch') == 0
        tx.Nch = 1;
    end

    % Number of points in constellation

    if isfield(tx,'modfor') == 0
        tx.modfor = 'qpsk';
    end

    switch tx.modfor
        case "bpsk"
            tx.length_alphabet  = 2;
        case "qpsk"
            tx.length_alphabet  = 4;
        case "ook"
            tx.length_alphabet  = 2;
        otherwise % [num]qam
            tx.length_alphabet  = str2double(tx.modfor(1:end-3));
    end

    % Number of Bits Per Symbol
    tx.Nbps = log2(tx.length_alphabet);

    tx = sort_struct_alphabet(tx);
end