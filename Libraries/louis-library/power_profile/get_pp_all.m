function varargout = get_pp_all(PPEparams,Axis,las,tx,ft,amp)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : PP - POWER PROFILE
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.@telecom-paris.com
%   Date            : 2023-01-24
%   Version         : 1.2
%
% ----- Main idea -----
%   Get the optical power profile in a point-to-point telecommunication
%   link thanks to the following steps:
%   1 - source + modulation
%   2 - pre distorsion
%   3 - Propagation after span 0 (amp+[fibre+amp]*Nspan)
%   4 - Chromatic Dispersion Compensation
%   5 - Chromatic dispersion for PPE
%   6 - Power Profile Estimation
% 
%   Accumulted dispersiona at each step
%   1 - 0
%   2 - pd+disp*length*Nspan
%   3 - pd
%   4 - 0
%   5 - disp*length*Nspan
%   6 - 0
% 
% ----- INPUTS -----
%   PPEparams:  (structure) containing the Power Profile Estimator 
%               parameters.
%               * FINDLOC --- See MOVING_AVERAGE.m for details.
%                   - AV_METHOD  (string)
%                   - AV_PERIOD  (scalar)[]
%
%               * LINK --- See SET_PPEPARAMS.m function for details
%                   - DL             (scalar)[m] 
%                   - NSTEPS_FIBRE   (scalar)[]
%                   - NSTEPS_PD      (scalar)[]
%                   - NSTEPS_TOT     (scalar)[]
%
%               * METHOD --- See PPE.m function for details.
%                   ** PP --- 
%                       - Q  (string)
%                       - M  (string)
%                   -  WF (string)
%                   
%               * PHYS --- See SET_PPEPARAMS.m function for details
%                   - PD         (scalar)[ps/nm]
%                   - ATT_FACTOR (scalar)[]
%                   - NL_FACTOR  (scalar)[1/W/km] or [rad/W]
%
%               * PLOT --- See SET_PPEPARAMS.m function for details
%                   - NORM  (boolean)
%                   - PLOT  (boolean)
%                   - DIST  (array)[m]
%
%               * REPET --- See SET_PPEPARAMS.m function for details
%                   - NTRIES (scalar)[]
%                   - SEED   (boolean)[]
%                   - PARFOR (boolean)[]
%
%   AXIS   (structure) containing the Axis parameters
%          See SET_AXIS.m function for details
%               - NSYMB     (scalar)[]
%               - NT        (scalar)[]
%               - SYMBRATE  (scalar)[GBAUD]
%               - NSAMP     (scalar)[]
%               - FS        (scalar)[GHZ]
%               - DF        (scalar)[GHZ]
%               - FREQ      (array)[GHZ]
%               - DT        (scalar)[ns]
%               - TIME      (array)[ns]
%                 
%   LAS   (structure) containing the Laser parameters
%               - LAM       (scalar)[nm]      carrier wavelength
%               - LINEWIDTH (scalar)[GHz]     laser linewidth
%               - N0        (scalar)[dB/GHz]  AWGN
%               - N_POLAR   (scalar)[]        number of polarisations
%               - PDBM      (scalar)[dBm]     power
%               - PLIN      (scalar)[mW]      power
%
%   FT   (structure) containing the Fibre parameters
%        See SET_FT.m function for details
%               - ALPHAdB   (scalar)[dB/km]
%               - ALPHALIN  (scalar)[1/m]
%               - LENGTH    (scalar)[m]
%               - DISP      (scalar)[ps/nm/km]
%               - SLOPE     (scalar)[ps/nm²/km]
%               - AEFF      (scalar)[µm²]
%               - GF        (scalar)[1/W/km]
%               - N2        (scalar)[m²/W]
%               - PMDPAR    (scalar)[ps/sqrt(km)]
%               - NPLATES   (scalar)[]
%               - ISMANAKOV (boolean)[]
%
%   AMP  (structure) containing the amplifiers parameters
%               - GAIN  (scalar)[dB] Power gain at each span
%               - NSPAN (scalar)[]   Number of spans in the link
%               - F     (scalar)[dB] Noise added by each amplifier
%
% ----- OUTPUTS -----
% VARARGOUT (cell array)
% If NTRIES > 1
%   -----------------------------------------------------------------------
%   | ---------------  ---------------   ---------------  --------------- |
%   | |  Ntry = 1   |  |             |   |             |  |Mean_PP_trunc| |
%   | |     ***     |  |  Mean_PP    |   |Mean_PP_trunc|  |     _av     | | 
%   | |  Ntry = M   |  |lgt=PD+FIBRE |   | lgt = FIBRE |  | lgt = FIBRE | | 
%   | |lgt=PD+FIBRE |  |             |   |             |  |             | |
%   | ---------------  ---------------   ---------------  --------------- |
%   ----------------------------------------------------------------------
% 
% Else
% 
%   -----------------------------------------------------
%   | ---------------  ---------------  --------------- |
%   | |             |  |             |  |             | |
%   | |     PP      |  |  PP_trunc   |  | PP_trunc_av | | 
%   | |lgt=PD+FIBRE |  | lgt = FIBRE |  | lgt = FIBRE | | 
%   | |             |  |             |  |             | |
%   | ---------------  ---------------  --------------- |
%   -----------------------------------------------------
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

%% output initialisation
if PPEparams.repet.Ntries > 1
    % --- back propagation only in the link
    if strcmp(PPEparams.link.method_BP,'no_pd') == 1
        varargout = cell(1,3);

        % PP_all
        varargout{1,1}{1} = zeros(1,PPEparams.link.nsteps_fibre);
        
        % Mean_PP
        varargout{1,2}{1} = zeros(1,PPEparams.link.nsteps_fibre);

        % Mean_PP_av
        if strcmp(PPEparams.findloc.av_method,'mirror') == 1
            varargout{1,3}{1} = zeros(1,PPEparams.link.nsteps_fibre);
        else
            varargout{1,3}{1} = zeros(1,PPEparams.link.nsteps_fibre-PPEparams.findloc.av_period);
        end
    else % --- back propagation in link + predisp
        varargout   = cell(1,4);

        % PP_all
        varargout{1,1}{1} = zeros(PPEparams.repet.Ntries,PPEparams.link.nsteps_BP);

        % Mean_PP
        varargout{1,2}{1} = zeros(1,PPEparams.link.nsteps_BP);

        % Mean_PP_trunc
        varargout{1,3}{1} = zeros(1,PPEparams.link.nsteps_fibre);

        % Mean_PP_trunc_av
        if strcmp(PPEparams.findloc.av_method,'mirror') == 1
            varargout{1,4}{1} = zeros(1,PPEparams.link.nsteps_fibre);
        else
            varargout{1,4}{1} = zeros(1,PPEparams.link.nsteps_fibre-PPEparams.findloc.av_period);
        end
    end
else % Ntries == 1
    % --- back propagation only in the link
    if strcmp(PPEparams.link.method_BP,'no_pd') == 1
        varargout   = cell(1,2);

        % PP_all
        varargout{1,1}{1} = zeros(1,PPEparams.link.nsteps_fibre);

        % PP_av
        if strcmp(PPEparams.findloc.av_method,'mirror') == 1
            varargout{1,2}{1} = zeros(1,PPEparams.link.nsteps_BP);
        else
            varargout{1,2}{1} = zeros(1,PPEparams.link.nsteps_BP-PPEparams.findloc.av_period);
        end
    else % --- back propagation in link + predisp
        varargout   = cell(1,3);

        % PP_all
        varargout{1,1}{1} = zeros(1,PPEparams.link.nsteps_BP);

        % PP_trunc
        varargout{1,2}{1} = zeros(1,PPEparams.link.nsteps_fibre);

        % PP_trunc_av
        if strcmp(PPEparams.findloc.av_method,'mirror') == 1
            varargout{1,3}{1} = zeros(1,PPEparams.link.nsteps_fibre);
        else
            varargout{1,3}{1} = zeros(1,PPEparams.link.nsteps_tot-PPEparams.findloc.av_period);
        end
    end
end

%% parameters init
length_pd       = PPEparams.phys.pd/ft.disp*1e3;    % [m]
length_link     = amp.Nspan*ft.length;              % [m]
length_total    = length_pd+length_link;            % [m]

% those three lines SHOULD NOT be removed, it is for the PARFOR loop
ft_tmp          = ft;
ft_tmp.length   = length_link;
D_in_ppe        = PPEparams.link.D_in_ppe;

%% raw power profiles
if PPEparams.repet.parfor == 0
    inigstate(Axis.Nsamp,Axis.fs);

    PP_all  = zeros(PPEparams.repet.Ntries,PPEparams.link.nsteps_BP);
    
    for k = 1:PPEparams.repet.Ntries

        %%% TRANSMITTER
        Ein         = TX(Axis,las,tx);
        Epd         = dpd(Ein,ft,PPEparams.phys.pd);

        %%% CHANNEL
        Elink       = channel(Epd,ft,amp);

        %%% RECEIVER
        Ecdc        = cdc(Elink,ft,length_total);

%         Einds   = ds(Ein,Axis.Nsps,1,1);
%         Ecdcds  = ds(Ecdc,Axis.Nsps,1,1);
%         [N0,c]  = getN0_MMSE(Ein,Ecdc)
%         plotfield(Ecdc,'--p-')
%         plot_const(Einds,Ecdcds)

        Einppe      = dpd(Ecdc,ft,D_in_ppe);
        [PP,Ecdcres]= get_ppe(Einppe,Ein,ft_tmp,PPEparams);

        %%% PROCESSING
        PP_all(k,:) = fliplr(PP);
        
        if PPEparams.plot.standardise == 1
            PP_all(k,:) = standardise(PP_all(k,:));
        end

%         Pin     = tx.Plin;
%         Ppd     = get_power(Epd);
%         Pout    = get_power(Elink);
%         Pres    = get_power(Ecdcres);
%         Mp      = [Pin,Ppd,Pout,Pres];
% 
%         Std_in  = get_stats(abs(Ein.field),'std','b');
%         Std_pd  = get_stats(abs(Epd.field),'std','b');
%         Std_out = get_stats(abs(Elink.field),'std','b');
%         Std_res = get_stats(abs(Ecdcres.field),'std','b');
%         Mstd    = [Std_in,Std_pd,Std_out,Std_res];
% 
%         writematrix(Mp,strcat('power --- ',num2str(PPEparams.phys.pd),' --- ',string(datetime)));
%         writematrix(Mstd,strcat('std --- ',num2str(PPEparams.phys.pd),' --- ',string(datetime)));
    end

    varargout{1,1}{1} = PP_all';

    if PPEparams.plot.plot == 1
        fields = pack_structs(Ein,Epd,Elink,Ecdc,Einppe,Ecdcres);

        file_name_w = sprintf("constellations --- fibre %s --- Pin = %i --- " + ...
            "dnu = %.0e --- n0 = %.0e --- Rs = %i --- %s --- RO = %.f --- %s.png",...
            ft.type,las.PdBm,las.linewidth*1e9,las.n0, ...
            Axis.symbrate,tx.pulse_shape,tx.rolloff,tx.modfor);
        custom_constellation_plot(fields,file_name_w)
        
        pause(1)
    end
else
    PP_all  = zeros(PPEparams.repet.Ntries,PPEparams.link.nsteps_BP);
    Nsamp   = Axis.Nsamp;
    fs      = Axis.fs;
    predisp = PPEparams.phys.pd;

    parfor k = 1:PPEparams.repet.Ntries

        inigstate(Nsamp,fs);

        Ein         = TX(Axis,las,tx);
        Epd         = dpd(Ein,ft,predisp);

        Elink       = channel(Epd,ft,amp);

        Ecdc        = cdc(Elink,ft,length_total);                    
        Einppe      = dpd(Ecdc,ft,D_in_ppe);

        PP          = get_ppe(Einppe,Ein,ft_tmp,PPEparams);
        PP_all(k,:) = fliplr(PP);

        if PPEparams.plot.standardise == 1
            PP_all(k,:) = standardise(PP_all(k,:));
        end
    end
    varargout{1,1}{1} = PP_all';
end

%% power profiles pre-processing

if PPEparams.repet.Ntries > 1
    % --- back propagation only in the link
    if strcmp(PPEparams.link.method_BP,'no_pd') == 1
        Mean_PP             = mean(PP_all);
        Mean_PP_av          = moving_average(Mean_PP,PPEparams.findloc);
        varargout{1,2}{1}   = Mean_PP;
        varargout{1,3}{1}   = Mean_PP_av;
    else % --- back propagation in link + predisp
        Mean_PP             = mean(PP_all);
        Mean_PP_trunc       = Mean_PP(PPEparams.link.nsteps_pd+1:end);
        Mean_PP_trunc_av    = moving_average(Mean_PP_trunc,PPEparams.findloc);
        varargout{1,2}{1}   = Mean_PP;
        varargout{1,3}{1}   = Mean_PP_trunc;
        varargout{1,4}{1}   = Mean_PP_trunc_av;
    end
else % Ntries = 1
% --- back propagation only in the link
    if strcmp(PPEparams.link.method_BP,'no_pd') == 1
        PP_av               = moving_average(PP_all,PPEparams.findloc);
        varargout{1,2}{1}   = PP_av;
    else % --- back propagation in link + predisp
        PP_trunc            = PP_all(PPEparams.link.nsteps_pd+1:end);
        PP_trunc_av         = moving_average(PP_trunc,PPEparams.findloc);
        varargout{1,2}{1}   = PP_trunc;
        varargout{1,3}{1}   = PP_trunc_av;
    end
end

