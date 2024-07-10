function Efilt = filt_elec(Ein,Axis,varargin)
%%
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : FILT_ELEC
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.com
%   ArXivs          : 
%   Date            : 2024-01-26: creation
%   Version         : 1.0
%
% ----- Main idea -----
%   Change the units of the Axis parameters
%
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------
%%

    if strcmpi(Axis.type,'optilux') == 1
        Axis = convert_scales_axis(Axis,"SI");
        Axis.ADC = convert_scales_axis(Axis.ADC,"SI");
    end

    % filter definition in frequency domain
    % filtering in frequency domain
    % back to the time domain
    % checking power levels

    if nargin == 2
        rx = struct();
        rx.CMA.n_polar = 2;
        rx.OEF.elec_rolloff = 0.07;
        rx.OEF.elec_bw  = Axis.symbrate;
        rx.OEF.elec_filt = "rect";
    end

    if isfield(rx.OEF,"elec_filt") == 0 || strcmp(rx.OEF.elec_filt,"rect") == 1
        H_filt  = rect_filter(Axis,rx);
    elseif  strcmp(rx.OEF.elec_filt,"RRC") == 1
        H_filt  = RRC_filter(Axis,rx);
    end

    Eft         = H_filt.*FFT(Ein,"array");        
    Efilt       = iFFT(Eft,'struct');
    Efilt.Nsps  = Ein.Nsps;

    % if wanna check the spectra
%     plot_spectrum(Ein,Axis)
%     plot_spectrum(Efilt,Axis)

    Diff_pows   = sum(abs(get_power(Ein)-get_power(Efilt)));
    ErrRel_pows = Diff_pows/sum(get_power(Ein));

    assert(ErrRel_pows<1e-1,"power mismatch")
%-----------------------------------------------------



function H_Rect = rect_filter(Axis,rx)

    f                           = Axis.ADC.freq;

    f_thresh_low                = (1+rx.OEF.elec_rolloff)*Axis.symbrate/2;
    f_thresh_high               = (1+rx.OEF.elec_rolloff)*Axis.symbrate/2;

    [~,f_low_neg]               = find(f >= - f_thresh_low); 
    [~,f_high_neg]              = find(f >= - f_thresh_high);
    [~,f_low_pos]               = find(f >= f_thresh_low);
    [~,f_high_pos]              = find(f >= f_thresh_high);

    f_low_neg                   =  f_low_neg(1);
    f_high_neg                  =  f_high_neg(1);
    f_low_pos                   =  f_low_pos(1);
    f_high_pos                  =  f_high_pos(1);

    H_Rect                      = zeros(Axis.ADC.Nsamp,1);
    H_Rect(1:f_high_neg)        = 0.0;
    H_Rect(f_high_pos:end)      = 0.0;
    H_Rect(f_low_neg:f_low_pos) = 1.0;

    H_Rect                      = repmat(H_Rect,[1,rx.CMA.n_polar]);
%-----------------------------------------------------

function H_RRC = RRC_filter(Axis,rx)
    f       = Axis.freq;

    % raw Raised Cosine
    NORM    = pi/(rx.OEF.elec_rolloff*rx.OEF.elec_bw);
    ARG     = abs(f)-(1-rx.OEF.elec_rolloff)*rx.OEF.elec_bw/2;
    H_RC    = (1+cos(NORM*ARG))/2;

    % Cutting frequencies
    f_thresh_low    = (1-rx.OEF.elec_rolloff)*Axis.symbrate/2;
    f_thresh_high   = (1+rx.OEF.elec_rolloff)*Axis.symbrate/2;

    [~,f_low_neg]   = find(f >= - f_thresh_low); 
    [~,f_high_neg]  = find(f >= - f_thresh_high);
    [~,f_low_pos]   = find(f >= f_thresh_low);
    [~,f_high_pos]  = find(f >= f_thresh_high);

    f_low_neg       =  f_low_neg(1);
    f_high_neg      =  f_high_neg(1);
    f_low_pos       =  f_low_pos(1);
    f_high_pos      =  f_high_pos(1);

    % modified RC
    H_RC(1,1:f_high_neg)        = 0.0;
    H_RC(1,f_high_pos:end)      = 0.0;
    H_RC(1,f_low_neg:f_low_pos) = 1.0;

    % Root Raised Cosine
    H_RRC   = transpose(sqrt(H_RC));
    H_RRC   = repmat(H_RRC,[1,rx.CMA.n_polar]);
%-----------------------------------------------------
