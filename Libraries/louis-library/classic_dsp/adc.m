function [Eadc,Axis]   = adc(Ein,Axis,rx)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : ADC - Analog to Digital Converter
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   ArXivs          : 2024-01-11: cleaning Yves' code
%                   : 2024-01-17: including Élie's code
%   Date            : 2024-01-22: adding other singularity management algorithm
%   Version         : 2.2
%
% ----- Main idea -----
%   Polarisation Mode Dispersion compensation using CMA algorithm.
%
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
%   Functions           : ZEROPADDING
%   Author              : Elie AWWAD
%   Author contact      : elie.awwad@telecom-paris.fr
%   Date                : NA
%   Title of program    : PHENIX
%   Code version        : 
%   Type                : Optical simulator toolbox - source code
%   Web Address         : NA
% ---------------------------------------------

    if strcmp(Axis.type,"optilux") == 1
        Axis = convert_scales_axis(Axis);
    end

    if strcmpi(rx.data,"exp") == 1 && isfield(Axis,"ADC") == 1
        ADC = Axis.ADC;
    end

    ADC.symbrate = Axis.symbrate;
    
    if size(Ein.field,2) == 1
        Xout    = Ein.field(:,1);
        Rx_ADC  = downsample(Xout, floor(Axis.Nsps/rx.ADC.Nsps));
    else
             
        if strcmp(rx.ADC.method,'downsampling') == 1
            Eadc = ds(Ein,Axis.Nsps,2,"replace");

        elseif strcmp(rx.ADC.method,'zeropadding') == 1
            Xout        = Ein.field(:,1);
            Yout        = Ein.field(:,2);

            if strcmp(rx.data,"exp") == 1
                if isfield(Axis,"ADC") == 0 || isfield(Axis.ADC,"fs_osc") == 0
                    ratio   = Axis.symbrate*2/Axis.fs;
                else
                    ratio   = Axis.symbrate*2/Axis.ADC.fs_osc;
                end
            else
                ratio   = Axis.symbrate*2/Axis.fs;
            end

            Rx_ADC_X    = zeropadding(Xout, ratio);
            Rx_ADC_Y    = zeropadding(Yout, ratio);
            Rx_ADC      = [Rx_ADC_X,Rx_ADC_Y];

            % managing the number of symbols after resampling
            if isfield(rx.ADC,'Nsps') == 0
                rx.ADC.Nsps = 2;
            end
            
            Nsymb_tmp   = length(zeropadding(Xout, ratio))/rx.ADC.Nsps;
            Nsymb_ADC   = Axis.Nsymb;

            if Nsymb_tmp >= Axis.Nsymb
                if mod(Nsymb_ADC,2) ~= 0
                    Nsymb_ADC = Nsymb_ADC-1;
                end
                ADC.Nsymb   = Nsymb_ADC;
                ADC.Nsamp   = rx.ADC.Nsps*ADC.Nsymb;
                Rx_ADC      = Rx_ADC(1:ADC.Nsamp,:);
            end
        end
    end

    Eadc.lambda = Ein.lambda;
    Eadc.field  = Rx_ADC;

    % new axis
    ADC.fs      = ratio*Axis.fs;
    ADC.dt      = 1/ADC.fs;
    if abs(ADC.fs - round(ADC.fs,0))<1e-5
        ADC.fs  = round(ADC.fs,0);
    end
    ADC.df      = ADC.fs/ADC.Nsamp;

    ADC.time    = 0:ADC.dt:(ADC.Nsamp-1)*ADC.dt;
    ADC.freq    = (-ADC.fs/2:ADC.df:ADC.fs/2-ADC.df)+ADC.df/2;

    ADC.tmax    = ADC.time(end);

    % +df/2 in order to center the freq-axis

    ADC.Nsps        = 2;
    Axis.ADC        = ADC;
    Eadc.Nsps       = ADC.Nsps;
    Axis.ADC.type   = "SI";


    Axis        = convert_scales_axis(Axis,"optilux");    
    Axis.ADC    = convert_scales_axis(Axis.ADC);

    assert(Axis.type==Axis.ADC.type, "units mismatch between TX and ADC axis")
end


%% TRASH
% ROUNDING ERRORS CAUSED BY LINSPACE
%     ADC.time    = linspace(0,ADC.Nsamp*ADC.dt,ADC.Nsamp);
%     ADC.freq    = linspace(-1,1,ADC.Nsamp)*ADC.fs/2;
