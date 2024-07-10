function amp = set_topology(tx,ft,amp,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : SET_TOPOLOGY
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.com
%   Date            : 2023-01-30
%   Version         : 1.0
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------

    amp = set_amp_int(amp,ft,tx);

    if isfield(amp,'losses') == 1
        amp = rmfield(amp,'losses');
    end

    for k=1:amp.Nspan
        amp.losses.(strcat('span',num2str(k))) = [0,0];
    end
    amp.losses.ids = 1:amp.Nspan;

    % if losses in the link
    if nargin > 3
        losses = varargin{1};
        assert(size(losses,2) == 3, ...
            "'LOSSES matrix should be [span n0,distance,loss]")

        losses(:,2) = losses(:,2)*1e3;

        for k = 1:size(losses,1)
            amp.losses.(strcat('span',num2str(losses(k,1)))) = losses(k,2:3);
        end
    end

    amp         = sort_struct_alphabet(amp);
    amp.losses  = sort_struct_alphabet(amp.losses);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : SET_AMP_INT - AMPLIFIERS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-03-15
%   Version         : 1.4
%
% ----- Main idea -----
%   Set the AMP structure from the amplifier type and fibre parameters
% 
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

function amp = set_amp_int(amp,ft,tx)
    
    %%% NF (Noise Figure)
    if strcmp(amp.type,"noAWGN") == 1
        amp.f = 0;
    elseif strcmp(amp.type,'ideal') == 1
        amp.f = 3;
    elseif strcmp(amp.type,'classic') == 1
        amp.f = 5;
    end

    %%% G/P.C (Gain/Constant) modes - PC by default
    if isfield(amp,'mode') == 0
        amp.mode    = "PC"; % Constant Power
    else
        if strcmp(amp.mode,'GC') == 1

            if isfield(amp,'outpower') == 1
                amp = rmfield(amp,'outpower');
            end
            amp.gain= ft.length*ft.alphadB*1e-3;

        end
    end

    if strcmp(amp.mode,'PC') == 1
        amp.outpower    = tx.PdBm;
        amp.gain_dB_th  = ft.length*ft.alphadB*1e-3;
    end

    amp = sort_struct_alphabet(amp);