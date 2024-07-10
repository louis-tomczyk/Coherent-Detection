function varargout = rx_filt(Elink,Axis,rx,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : RX_FILT - Receiver side filtering
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   ArXivs          : 2022-09-09 : creation
%   Date            : 2024-01-23 : simplifying and increase robustness
%   Version         : 2.0
%
% ----- Main idea -----
%   Filter the input field to select the channel using the OPT_FILT
%   Then isolate the wanted channel from others using the ELEC_FILT.
%
% ----- INPUTS -----
%   ELINK   (structure) containing the Fields to be filtered
%               - LAMBDA    [nm]        wavelength
%               - FIELD     [sqrt(mW)]  Normalised electric fields
%   AXIS    (structure)     See SET_AXIS function for details
%   RX      (structure)     See SET_RX function for details
%
% ----- BIBLIOGRAPHY -----
%   Functions           : Nyquist_WDM_DP_QPSKand16QAM_Vfinal
%   Author              : Yves JAOUEN
%   Author contact      : yves.jaouen@telecom-paris.fr
%   Date                : Unknown
%   Title of program    : Code DP_QPSKand16QAM
%   Code version        : 2022
%   Type                : Optical simulator toolbox - source code
%   Web Address         : NA
% -----------------------
%   Articles
%   Author              :   Junyi WANG, Chongjin XIE, Zhongqi PAN
%   Title               :   Generation of Spectrally Efficient Nyquist-WDM
%                           QPSK Signals Using Digital FIR or FDE Filters
%                           at Transmitters
%   Jounal              :   Journal of Lightwave Technology
%   Volume - N°         :   30 - 23 
%   Date                :   2012-01-23
%   DOI                 :   10.1109/JLT.2012.2226207
% ---------------------------------------------

%     Eoptfilt    = opt_filt(Elink,Axis,rx);
    Eelecfilt   = elec_filt(Elink,Axis,rx);
    
    if nargin == 3
        [Eadc,Axis]     = adc(Eoptfilt,Axis,rx);
        Erx             = fn(Eadc);
        varargout{1}    = Erx;
        varargout{2}    = Axis;
    elseif nargin == 4
        if strcmp(varargin{1},"no adc")
            varargout{1}    = Eelecfilt;
        else
            [Eadc,Axis]     = adc(Eelecfilt,Axis,rx);
            Erx             = fn(Eadc);
            varargout{1}    = Erx;
            varargout{2}    = Axis;
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------
% ----- CONTENTS -----
% --- from YVES JAOUEN
%   Eoptfilt = opt_filt(Ein,Axis,rx)
%   Eelecfilt = elec_filt(Eoptfilt,Axis,rx)
%   - H_RRC = RRC_filter(Axis,rx)
% ---------------------------------------------

function Eoptfilt = opt_filt(Ein,Axis,rx)

    % filter definition in frequency domain
    % filtering in frequency domain
    % back to the time domain
    % checking power levels

    Hopt        = opt_filter(Axis,rx);
    Efiltft     = Hopt.*FFT(Ein,"array");
    Eoptfilt    = iFFT(Efiltft,'struct');

    Diff_pows   = sum(abs(get_power(Ein)-get_power(Eoptfilt)));
    ErrRel_pows = Diff_pows/sum(get_power(Ein));

    assert(ErrRel_pows<5e-2,"power mismatch")
%-----------------------------------------------------

function Efilt = elec_filt(Ein,Axis,rx)

    % filter definition in frequency domain
    % filtering in frequency domain
    % back to the time domain
    % checking power levels


    if isfield(rx.OEF,"elec_filt") == 0 || strcmp(rx.OEF.elec_filt,"rect") == 1
        H_filt  = rect_filter(Axis,rx);
    elseif  strcmp(rx.OEF.elec_filt,"RRC") == 1
        H_filt  = RRC_filter(Axis,rx);
    end

    Eft         = H_filt.*FFT(Ein,"array");        
    Efilt       = iFFT(Efiltft,'struct');

    Diff_pows   = sum(abs(get_power(Ein)-get_power(Efilt)));
    ErrRel_pows = Diff_pows/sum(get_power(Einfilt));

    assert(ErrRel_pows<1e-1,"power mismatch")
%-----------------------------------------------------

% =========================== %
% ---- SUB SUB FUNCTIONS ---- %
% =========================== %


function H_opt = opt_filter(Axis,rx)

    n_polar     = rx.CMA.txpolars;

    if abs(Axis.freq(1))>1e5
        scale = 1e-9;
    else
        scale = 1;
    end

    Hopt    = exp(-(Axis.freq*scale/rx.OEF.opt_bw).^(2*rx.OEF.opt_order));
    H_opt   = repmat(Hopt,[n_polar,1]).';
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

    n_polar = rx.CMA.txpolars;
    H_RRC   = repmat(H_RRC,[1,n_polar]);
%-----------------------------------------------------




function H_Rect = rect_filter(Axis,rx)

    f       = Axis.freq;

    % Cutting frequencies
    f_thresh_low    = Axis.symbrate/2;
    f_thresh_high   = Axis.symbrate/2;

    [~,f_low_neg]   = find(f >= - f_thresh_low); 
    [~,f_high_neg]  = find(f >= - f_thresh_high);
    [~,f_low_pos]   = find(f >= f_thresh_low);
    [~,f_high_pos]  = find(f >= f_thresh_high);

    f_low_neg       =  f_low_neg(1);
    f_high_neg      =  f_high_neg(1);
    f_low_pos       =  f_low_pos(1);
    f_high_pos      =  f_high_pos(1);

    H_Rect(1,1:f_high_neg)        = 0.0;
    H_Rect(1,f_high_pos:end)      = 0.0;
    H_Rect(1,f_low_neg:f_low_pos) = 1.0;

    n_polar = rx.CMA.txpolars;
    H_Rect   = repmat(H_Rect,[1,n_polar]);
%-----------------------------------------------------
