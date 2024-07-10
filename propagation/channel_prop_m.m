function Eoutlink = channel_prop_m(Ein,ft,amp)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : CHANNEL_PROP
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@enst.fr
%   Date            : 2023-01-12
%   Version         : 1.1.1
%
% ----- Main idea -----
%   Simulate a link of a telecommunication network with a given number of
%   spans in a constant mode operation with identical spans.
%   The names of the input and output fields is because we first propagate
%   the field into a fibre before this function.
%
% ----- INPUTS -----
%   EIN:    (structure) containing the Fields to propagate
%               - LAMBDA    (scalar)[nm]: wavelength
%               - FIELD     (array)[sqrt(mW)]: normalised electric fields
%   FT:     (structure) fiber parameters
%               - ALPHADB   (scalar)[dB/km]         Power attenuation
%               - DISP      (scalar)[ps/nm/km]      Dispersion
%               - SLOPE     (scalar)[ps/nm²/km]     Slope of the dispersion
%               - LENGTH    (scalar)[m]             Length
%               - n2        (scalar)[W²/m]          Non linear index
%               - pmdpar    (scalar)[ps/sqrt(km)]   Polarization Mode Dispersion
%   AMP:     (structure) containing the amplifiers parameters
%               - NSPAN     (scalar)[]              Number of spans in the
%                                                   link
%               - GAIN      (scalar)[dB]            Gain of each amplifier
%                                                   if constant mode
%               - F         (scalar)[dB]            Noise Figure
%   METHOD   (structure or sring)                   Select the propagation method:
%                                                   'SSFM' or 'RP1'
%
% ----- OUTPUTS -----
%  EOUTLINK:  (structure)                           Propagated field
%               - LAMBDA    (scalar) [nm]: wavelength
%               - FIELD     (array) [sqrt(mW)]: normalised electric fields
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    %%% maintenance: checking arguments
    argnames = string(nargin);
    for k=1:nargin
        argnames(k)  = inputname(k);
    end
    
    is_arg_missing('Ein',argnames);
    is_arg_missing('ft',argnames);
    is_arg_missing('amp',argnames);


    % if there are not losses
    if isfield(amp,'losses') == 0
        Eoutfib  = fiber(Ein,ft);    
        Eoutlink = amp_EDFA(Eoutfib,ft,amp);

    % else, if there are losses
    else
        for k = 1:length(amp.losses.ids)

            % if some of the losses are in the first span
            if sum(amp.losses.ids==1)==1

                dloss       = amp.losses.span1(:,1);
                loss_dB     = amp.losses.span1(:,2);
                if strcmp(amp.mode,'GC')==1
                    Eoutfib = span_with_loss_m(Ein,ft,dloss,loss_dB,amp);
                elseif strcmp(amp.mode,'PC')==1
                    PC_dBm  = get_power(Ein,struct('unit',"dBm",'polar','tot'));
                    Eoutfib = span_with_loss_m(Ein,ft,dloss,loss_dB,amp,PC_dBm);
                end
                amp.Nspan   = amp.Nspan-1;
                Eoutlink    = amp_EDFA(Eoutfib,ft,amp);
            % else if there are not in the first span 
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Eout,amp] = amp_EDFA(Eout,ft,amp)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : AMP_EDFA - Erbium Doped Fibre Amplifier
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@enst.fr
%   Date            : 2022-08-06
%   Version         : 1.1.1
%
% ----- Main idea -----
%   Simulate a link of a telecommunication network with a given number of
%   spans in a constant mode operation with identical spans.
%   The names of the input and output fields is because we first propagate
%   the field into a fibre before this function.
%
% ----- INPUTS -----
%   EOUT:    (structure) containing the Fields to be compensated
%               - LAMBDA    (scalar)[nm]: wavelength
%               - FIELD     (array)[sqrt(mW)]: normalised electric fields
%   FT:     (structure) fiber parameters
%               - ALPHADB   (scalar)[dB/km]         Power attenuation
%               - DISP      (scalar)[ps/nm/km]      Dispersion
%               - SLOPE     (scalar)[ps/nm²/km]     Slope of the dispersion
%               - LENGTH    (scalar)[m]             Length
%               - n2        (scalar)[W²/m]          Non linear index
%               - pmdpar    (scalar)[ps/sqrt(km)]   Polarization Mode Dispersion
%   AMP:     (structure) containing the amplifiers parameters
%               - NSPAN     (scalar)[]              Number of spans in the
%                                                   link
%               - GAIN      (scalar)[dB]            Gain of each amplifier
%                                                   if constant mode
%               - F         (scalar)[dB]            Noise Figure
%
% ----- OUTPUTS -----
%  AMP:     (structure) to which is added:
%           - LENGTH_LINK   (scalar)[m]         total length of the link
%  EOUT:    [structure] Chromatic Dispersion Compensated field
%               - LAMBDA    (scalar)[nm]            Wavelength
%               - FIELD     (array)[sqrt(mW)]       Normalised electric fields
%
% ----- BIBLIOGRAPHY -----
%   Functions           : FIBER - AMPLIFLAT
%   Author              : Paolo SERENA
%   Author contact      : serena@tlc.unipr.it
%   Date                : 2021
%   Title of program    : Optilux
%   Code version        : 2021
%   Type                : Optical simulator toolbox - source code
%   Web Address         : https://optilux.sourceforge.io/
% ---------------------------------------------

    ns      = 0;
    while ns<2*amp.Nspan
        % factor 2 as 1 span = propagation + amplification
        
        % we amplify only after propagation, so at each even steps
        if mod(ns,2) == 0
            Eout = ampliflat(Eout,amp);
%             fprintf('span number --- %i \n',floor(ns/2)+1)
    
        % we re-propagate after amplification, so at each odd steps
        elseif mod(ns,2) == 1 && ns ~= 2*amp.Nspan-1
            Eout = fiber(Eout,ft);
        end
    
        ns = ns+1;
    end

    amp.length_link = amp.Nspan*ft.length;
