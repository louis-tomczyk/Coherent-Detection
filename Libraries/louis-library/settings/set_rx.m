function Rx = set_rx(varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : SET_DSP - Digital Signal Processing
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   ArXivs          : <2022-10-03: creation
%   Date            : 2024-01-18: adding Phenix platform management
%   Version         : 2.0.1
%
% ----- Main idea -----
%   Set the DSP structure from the laser transmitter parameters and
%   optional dsp parameters
% 
% ----- INPUTS -----
%   DSP:        (strcture) containing DSP parameters - OPTIONAL
%               - WR        (bool)[]    = 1 if decision step
%               - CPC_AVG   (scalar)[]  window averaging length
%
% ----- OUTPUTS -----
%   DSP         (structure) containing the fields:
%               - CDC       (bool)[]    Chromatic Dispersion Compensation
%               - PMDC      (bool)[]    Polarisation Mode DC
%               - CPC       (bool)[]    Carrier Phase C
%               - CPC_AVG   (scalar)[]  Window averaging length
%               - WR        (bool)[]    Waveform Reconstruction
%
% ----- BIBLIOGRAPHY -----
%   Category            : Book
%   Author              : GANG-DING Peng
%   Title               : Handbook of Optical Fibers
%   Author contact      : NA
%   Date                : 2019
%   Editor              : Springer
%   DOI                 : 10.1007/978-981-10-7087-7
%   ISBN                : 978-981-10-7085-3
%   Pages/Equations     : 176 / 50-51
%   Title of program    : NA
%   Code version        : NA
%   Type                : NA
%   Web Address         : NA
% ---------------------------------------------

    %% MAINTENANCE
    argnames = string(nargin);
    for k=1:nargin
        argnames(k)  = inputname(k);
    end

    is_arg_missing('Axis',argnames);
    is_arg_missing('dsp',argnames);
    is_arg_missing('tx',argnames);
    is_arg_missing('las',argnames);
    is_arg_missing('ft',argnames);
    is_arg_missing('amp',argnames);

    Axis= varargin{argnames == 'Axis'};
    dsp = varargin{argnames == 'dsp'};
    tx  = varargin{argnames == 'tx'};
    las = varargin{argnames == 'las'};
    ft  = varargin{argnames == 'ft'};
    amp = varargin{argnames == 'amp'};

    if nargin == 6
        Rx      = struct();
    else
        Rx = varargin{argnames == 'rx'};
    end

    if isfield(Rx,'data') == 0
        Rx.data         = "simu";
    end

    %% Optical Electrical Filtering

    if isfield (dsp,'rxfilt') == 0
        dsp.rxfilt          = 1;

        Rx.OEF.opt_bw       = 1*Axis.symbrate;
        Rx.OEF.opt_order    = 16;

        Rx.OEF.elec_bw      = 1*Axis.symbrate;
        Rx.OEF.elec_rolloff = tx.rolloff;
    end

    Rx      = sort_struct_alphabet(Rx);
    
    %% ANALOG TO DIGITAL CONVERTER

    if isfield (dsp,'adc') == 0
        dsp.adc         = 1;
        Rx.ADC.Nsps     = 2;
        Rx.ADC.method   = "zeropadding";
    end

    Rx.ADC  = sort_struct_alphabet(Rx.ADC);
    Rx      = sort_struct_alphabet(Rx);

    %% CHROMATIC DISPERSION COMPENSATION
    if strcmp(ft.type,'B2B') == 1 
        dsp.cdc     = 0;
        amp.Nspan   = 0;
    end

    if dsp.cdc == 1
        if isfield(Rx,"CDC") == 0    
            Rx.CDC = struct();
        end
        if isfield(Rx.CDC,"method") == 0
            Rx.CDC.method = "fibre-adc";
        end
        length_pd       = tx.pd/ft.disp*1e3;
        length_link     = amp.Nspan*ft.length;
        Rx.CDC.len_tot  = length_pd+length_link;
    end

    if dsp.cdc == 1
        if sum(strcmp(argnames,'rx')) == 1
            rx = varargin{argnames == 'rx'};
            if isfield(Rx.CDC,'method') == 0
                rx.CDC.method = "fibre-optilux";
            else
                rx.CDC.method = "FIR";
            end
        end

        Rx.CDC = sort_struct_alphabet(Rx.CDC);
    end
    
    Rx = sort_struct_alphabet(Rx);


    %% POLARISATION MODE DISPERSION COMPENSATION
    if isfield(dsp,'mimo') == 1
        %%%
        if isfield(Rx,"CMA") == 0
            CMAparams = struct();
        else
            CMAparams = Rx.CMA;
        end

        %%%
        if isfield(CMAparams,'method') == 0
            CMAparams.method    = "simu";
        else
            CMAparams.method = Rx.CMA.method;
        end

        %%%
        if isfield(CMAparams,"Ntrain") == 0
            CMAparams.Ntrain    = Axis.Nsymb;
        else
            CMAparams.Ntrain    = Rx.CMA.Ntrain;
        end

        %%%
        if isfield(CMAparams,"method_init") == 0
            CMAparams.method_init    = "dirac";
        else
            CMAparams.method_init    = Rx.CMA.method_init;
        end

        CMAparams.Rv        = [];
        CMAparams.step      = 2; % 2 samples per symbol

        CMAparams.R         = [1,1];
        CMAparams.lr_train  = 1e-3;
        CMAparams.lr_eq     = 1e-3;

        lll     = ft.length*amp.Nspan;
        Ncd     = 2*pi*abs(ft.beta2)*lll*(tx.Nbps*Axis.symbrate*1e9)^2;
        Npmd    = 2*pi*abs(ft.beta2)*sqrt(lll)*(tx.Nbps*Axis.symbrate*1e9)^2;
        Ntaps   = ceil((Ncd+Npmd)/10);

        % round up to next dozen
        Ntaps   = ceil(Ntaps/10)*10;

        if Ntaps < 20
            Ntaps = 21;
        end
        % ensure odd number or taps
        if mod(Ntaps,2) ==0
            Ntaps = Ntaps+1;
        end

        if isfield(CMAparams,'flag') == 0
            CMAparams.flag.thresh = 1/50;
        end
        
        CMAparams.taps      = Ntaps;
        CMAparams.txpolars  = las.n_polar;
        CMAparams.phizero   = 0;        % rotation angle between tx and rx field
        CMAparams.eps       = 1e-3;
        CMAparams.plot      = 0;

        if isfield(CMAparams,"singularity") == 0
            CMAparams.singularity.method    = "rand";
            CMAparams.singularity.thresh    = 0.3;
            CMAparams.singularity.loops     = 1;
        else
            CMAparams.singularity           = Rx.CMA.singularity;
        end

        if isfield(CMAparams.singularity,"thresh") == 0
            CMAparams.singularity.thresh    = 0.3;
        end

        CMAparams.singularity.flag  = false;
        CMAparams.Nsps_in           = 2;
        CMAparams.Nsps_out          = 1;
        CMAparams.plot_final        = 1;

        if strcmp(CMAparams.method,"phenix") == 1
            CMAparams.Rv   = []; % for MMA
            CMAparams.step = rat(CMAparams.Nsps_in);
        end


        Rx.CMA              = CMAparams;
        Rx.CMA              = sort_struct_alphabet(Rx.CMA);
        Rx.CMA.singularity  = sort_struct_alphabet(Rx.CMA.singularity);
    end

    Rx = sort_struct_alphabet(Rx);

    %% CARRIER PHASE COMPENSATION
    if sum(strcmp(argnames,'rx')) == 1
        rx = varargin{argnames == 'rx'};
        if isfield(rx,'CPC') == 0
            Rx.CPC.navg = 51;
        else
            Rx.CPC.navg = rx.CPC.navg;
        end
    else
        Rx.CPC.navg = 51;
        Rx.CPC = sort_struct_alphabet(Rx.CPC);
    end

    Rx  = sort_struct_alphabet(Rx);
    
    %% NON LINEAR NOISE ESTIMATION
    if isfield(dsp,"nlne") == 1
        Rx.NLN.plot     = 1;
        Rx.NLN.stats    = "var";
    end

    Rx = sort_struct_alphabet(Rx);

end