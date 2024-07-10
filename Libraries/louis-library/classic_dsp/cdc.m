function Ecdc = cdc(Ein,Axis,ft,rx)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : CDC - CHROMATIC DISPERSION COMPENSATION
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   ArXivs          : 2023-01-06: creation
%                   : 2024-01-17: including Yves' code
%   Date            : 2024-01-23: adding the "fibre-adc" option
%   Version         : 2.1
%
% ----- Main idea -----
%   Chromatic compensation using ideal Dispersion Compensating Fiber (DCF).
%   The chromatically dispersed field is propagated into a DCF having the
%   exact opposite dispersion characteristics (D [ps/nm/km], S [ps/nm²/km])
%   but without any (non) linear effects affecting the propagation.
%
%   If the method is "fibre-adc", then we allow cdc after the downsampling
%   to Nsps = 2
%
% ----- INPUTS -----
%   EIN:    [structure] containing the Fields to be compensated
%           structure should be organised as:
%               - LAMBDA [nm]: wavelength
%               - FIELD [sqrt(mW)]: normalised electric fields
%   FT:     [structure] fiber parameters
%           structure should at least contain:
%               - ALPHADB   [dB/km]         Power attenuation
%               - DISP      [ps/nm/km]      Dispersion
%               - SLOPE     [ps/nm²/km]     Slope of the dispersion
%               - LENGTH    [m]             Length
%               - n2        [W²/m]          Non linear index
%               - pmdpar    [ps/sqrt(km)]   Polarization Mode Dispersion
%   DL:     [scalar] [m] length of the fiber you want CDC
%
% ----- OUTPUTS -----
%  ECDC:    [structure] Chromatic Dispersion Compensated field
%           structure containing
%               - LAMBDA    (scalar)    [nm]: wavelength
%               - FIELD     (array)     [sqrt(mW)]: normalised electric fields
%
% ----- BIBLIOGRAPHY -----
%   Functions   : FIBER - PBC - PBS
%   Author              : Paolo SERENA
%   Author contact      : serena@tlc.unipr.it
%   Date                : 2021
%   Title of program    : Optilux
%   Code version        : 2021
%   Type                : Optical simulator toolbox - source code
%   Web Address         : https://optilux.sourceforge.io/
% ---------------------------------------------
   
flag    = check_cdc(ft,rx);

if flag == 1
    fc      = set_fc(ft,rx);
    Ecdc    = CDC(Ein,Axis,fc,rx);
else
    Ecdc    = Ein;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------
% ----- CONTENTS -----
%   Ecdc = check_cdc(Ein,ft,rx)
%   varargout = set_fc(ft,rx)
%   Ecdc = CDC(Ein,Axis,fc,rx)
%   - H_CD = CD_filter(Axis,fc,rx)
% ---------------------------------------------

function flag = check_cdc(ft,rx)

    flag = 1;
    if isfield(rx,"CDC") == 0
        flag = 0;

    end
    
    % if no dispersion or no length, then return the input
    if ft.disp == ft.slope && ft.disp == 0 || rx.CDC.len_tot == 0    
        flag = 0;
    end
% ---------------------------------------------

function fc = set_fc(ft,rx)
    
    if strcmp(rx.CDC.method,"FIR") ~= 1
    
        % DCF properties
        fc.type     = "ideal DCF";
        fc.length   = rx.CDC.len_tot;
        fc.lambda   = ft.lambda;
    
        fc.alphadB  = 0;
        fc.alphaLin = 0;
    
        fc.disp     = -ft.disp;
        fc.slope    = -ft.slope;
        fc.beta2    = -ft.beta2;
        fc.beta3    = -ft.beta3;
        
        fc.coupling = 'none';
        fc.pmdpar   = 0;
        fc.nplates  = 1;
        fc.ismanakov= false;
        
        fc.n2       = 0;
        fc.gf       = 0;
        fc.aeff     = ft.aeff;
    else
        fc = ft;
    end

    fc  = sort_struct_alphabet(fc);
% ---------------------------------------------

function Ecdc = CDC(Ein,Axis,fc,rx)

if strcmp(rx.CDC.method,"fibre-optilux") == 1
    
    if size(Ein.field,2) == 2*size(Ein.lambda,1)% n_polar = 2;

        % split the polarisations
        % compensate each polarisation
        % recombine both polarisations
        
        [Einx,Einy] = sep_XYfields(Ein);
        Ecdcx       = fiber(Einx,fc);
        Ecdcy       = fiber(Einy,fc);
        Ecdc        = merge_XYfields(Ecdcx,Ecdcy);
        
    
    else% n_polar = 1;
        Ecdc   = fiber(Ein,fc);
    end

elseif strcmp(rx.CDC.method,"fibre-adc") == 1

    % filter definition in frequency domain
    % filtering in frequency domain
    % back to the time domain
    % checking power levels

    H_CDC       = CDC_filter_in(Axis,fc,rx);
    Ecdcft      = H_CDC.*FFT(Ein,"array");
    Ecdc        = iFFT(Ecdcft,'struct');

    Diff_pows   = sum(abs(get_power(Ein)-get_power(Ecdc)));
    ErrRel_pows = Diff_pows/sum(get_power(Ein));
    assert(ErrRel_pows<1e-10,"power mismatch")
else
    Ecdc = GVD_mitigation_OS(Eadc,ft,Axis,rx);
end
% ---------------------------------------------

function H_CDC = CDC_filter_in(Axis,fc,rx)
    
    n_polar = size(rx.CMA.txpolars,2);

    f       = Axis.ADC.freq*1e9;
    b2      = fc.beta2;
    L       = rx.CDC.len_tot;

    H_CDC   = exp(1i*b2/2*(2*pi*f).^2*L);
    H_CDC   = repmat(H_CDC,[n_polar,1]).';
% ---------------------------------------------

