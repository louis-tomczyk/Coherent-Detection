function [Eoutlink,varargout] = channel_prop(Ein,ft,amp,method,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : CHANNEL_PROP
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   ArXivs          : 2023-01-12
%   Date            : 2024-01-29: including 'only-MIMO' option
%   Version         : 1.2
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
    try
        is_arg_missing('Ein',argnames);
    catch
        is_arg_missing('Epd',argnames);
    end
    is_arg_missing('ft',argnames);
    is_arg_missing('amp',argnames);

    %%% method selection
    % if METHOD is a structure
    if isstruct(method)==1
        if strcmp(method.calc,'SSFM') == 1
            Eoutfib     = fiber(Ein,ft);
        elseif strcmp(method.calc,'RP1') == 1
            RPparams    = method.params;
            Axis        = method.axis;
            Eoutfib     = RP1flex(Ein,ft,Axis,RPparams);
        end
    % if METHOD is a string
    else
        if strcmp(method,'SSFM') == 1
            Eoutfib     = fiber(Ein,ft);
            varargout{1}= [];
        
        elseif strcmp(method,'RP1') == 1
            RPparams    = method.params;
            Axis        = method.axis;
            Eoutfib     = RP1flex(Ein,ft,Axis,RPparams);
            varargout{1}= [];

        elseif strcmp(method,'only-MIMO') == 1
            Axis        = varargin{1};

            [V,h11,h12,h21,h22,Y] = pmd_emulator(Ein.field,Axis.freq,ft.pmdpar,0,1);
            Eoutfib         = Ein;
            Eoutfib.field   = Y;

            ft.Hplates  = V;
            ft.Hmimo.HH = h11;
            ft.Hmimo.HV = h12;
            ft.Hmimo.VH = h21;
            ft.Hmimo.VV = h22;

            varargout{1}=ft;

        elseif strcmp(method,'CD-MIMO') == 1
            Axis            = varargin{1};
            [Eoutfib,ft]    = CD(Ein,Axis,ft);

            [V,h11,h12,h21,h22,Y] = pmd_emulator(Eoutfib.field,Axis.freq,ft.pmdpar,0,1);
            Eoutfib         = Ein;
            Eoutfib.field   = Y;

            ft.Hplates  = V;
            ft.Hmimo.HH = h11;
            ft.Hmimo.HV = h12;
            ft.Hmimo.VH = h21;
            ft.Hmimo.VV = h22;
            
            varargout{1}=ft;
        end
    end
    
    if amp.Nspan ~= 0
        Eoutlink = amp_EDFA(Eoutfib,ft,amp);
    else
        Eoutlink = Eoutfib;
    end
    Eoutlink.Nsps = Ein.Nsps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Eout,amp] = amp_EDFA(Eout,ft,amp)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : AMP_EDFA - Erbium Doped Fibre AmplifierS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
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
