function PPEparams = set_PPEparams_demo(tx,ft,amp,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name      SET_PPEPARAMS
%   Author             louis tomczyk
%   Institution        Telecom Paris
%   Email              louis.tomczyk@telecom-paris.fr
%   Date               2023-03-21
%   Version            2.2.1
%
% ----- Main idea -----
%   Set the POWER PROFILE ESTIMATOR strcucture of parameters
% 
% ----- INPUTS -----
%   PPEparams   (structure) containing the PPE parameters.
%               * FINDLOC --- Find Location of extrema in the  pp by averaging the  pp.
%                   - AV_METHOD     (string)    Method to average the profile using a moving average
%                   - AV_PERIOD     (scalar)[]  Number of elements for the moving average
%                                               See MOVING_AVERAGE.m
%
%               * METHOD --- Mathematical part of the  ppe.
%                   ** PP --- Set parameters for the  ppe
%                       - Q         (string)    QUANTITY   for the Mathematical operation
%                       - M         (string)    MATHEMATICAL operation
%                   -  WF           (string)    WAVEFORM reconstruction method,See PPE.m
%
%               * PHYS --- Physical parameters for the  ppe
%                   -   PD          (scalar)    PREDISPERSION to be applied to the modulated laser
%                                   [ps/nm]     before channel.
%
%               * PLOT --- Ploting parameters
%                   - NORM          (boolean)   center and reduce variable
%                   - PLOT          (boolean)   plot intermediate constellations
%
%               * REPET --- Repetition parameters
%                   - NTRIES        (scalar)[]  number of realisations to emulate real data flow
%
%   FT          (structure) containing the Fibre parameters, see SET_FT.m
%   AMP         (structure) containing the amplifiers parameters, see SET_AMP.m/SET_TOPOLOGY.m
%
% ----- OUTPUTS -----
%   PPEparams   (structure) containing the PPE parameters structure to which will/might be added 
%               * PHYS ---
%                   - ATT_FACTOR    (scalar)[]  tuning attenuation in fibre
%                   - NL_FACTOR     (scalar)    the non linear parameter
%                                   [1/W/km] or [rad/W]
%               * PLOT ---
%                   - DIST          (array)[m]  propagation axis
%
%               * REPET ---
%                   - SEED          (boolean)   enable PRBS
%                   - PARFOR        (boolean)   enable parallel computing
%
%               * LINK --- number of iterations related parameters
%                   - DL            (scalar)[m] elementary step in the link
%                   - NSTEPS_FIBRE  (scalar)[]  number of steps related tothe 'true' propagation 
%                                               into the fibres
%                   - NSTEPS_PD     (scalar)[]  number of steps related to the 'fake' propagation
%                                               into a PreDispersion fibre
%                   - NSTEPS_TOT    (scalar)[]  sum of the 2 previous
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    if nargin >= 4
        PPEparams = varargin{1};
    else
        PPEparams = struct();
    end

    %% POWER PROFILE ESTIMATOR
    %%% Method
    if isfield(PPEparams,"method") == 0             % q-mod / m-cor is best
        PPEparams.method.pp.q       = 'mod';        % choose among pow-mod-field
        PPEparams.method.pp.m       = 'cor';        % choose among pow-cov-cor
        PPEparams.method.wf         = 'tanimura';
    else
        if strcmp(PPEparams.method.pp.q,'mod')  == 0 && ...
           strcmp(PPEparams.method.pp.m,'cor')  == 0 && ...
           strcmp(PPEparams.method.wf,'mod')    == 0
        disp("CAUTION: it is not the Tanimura estimator that is going to be used.")
        end
    end

    %%% Physics
    if isempty(strfind(PPEparams.method.wf,'tanimura'))
        PPEparams.phys.nl_factor    = ft.gf;
    else
        PPEparams.phys.nl_factor    = 0.01;
    end

    PPEparams.phys.pd = tx.pd;
    if strcmp(PPEparams.method.wf,'tanimura') ~= 1 && isfield(PPEparams.phys,'att_factor') == 0
        PPEparams.phys.att_factor   = 1.9958;
    end

    %%% plot and mathematical operations
    if isfield(PPEparams,'plot') == 0
        PPEparams.plot = struct();
    end

    if isfield(PPEparams.plot,"standardise") == 0
        PPEparams.plot.standardise = 0;
    end

    if strcmp(PPEparams.method.pp.m,"pow") == 1
        PPEparams.method.pp.q = "pow";
    end

    %%% link
    if isfield(PPEparams,"link") == 0
        PPEparams.link.dl = 1e3;    % [m]
    end
    
    ratio_length    = 1e3/PPEparams.link.dl;
    Nspans_tot      = amp.Nspan;

    if amp.Nspan == 0
        PPEparams.link.nsteps_fibre = floor(ft.length/PPEparams.link.dl);
        PPEparams.plot.dist = linspace(0,ft.length*1e-3,PPEparams.link.nsteps_fibre);
    else
        PPEparams.link.nsteps_fibre = floor(Nspans_tot*ft.length/PPEparams.link.dl);
        PPEparams.plot.dist = linspace(0,Nspans_tot*ft.length*1e-3,PPEparams.link.nsteps_fibre);
    end

    PPEparams.link.nsteps_span = floor(ft.length/PPEparams.link.dl);

    if isfield(PPEparams.plot,"plot") == 0
        PPEparams.plot.plot = 0;
    end

    %%% back propagation
    PPEparams.link.nsteps_pd        = floor(PPEparams.phys.pd/ft.disp*ratio_length);
    PPEparams.link.nsteps_tot       = floor(PPEparams.link.nsteps_fibre+PPEparams.link.nsteps_pd);

    if isfield(PPEparams.link,"method_BP") == 1
        if strcmp(PPEparams.link.method_BP,'no_pd')
            PPEparams.link.D_in_ppe = Nspans_tot*ft.length*1e-3*ft.disp;
        else
            PPEparams.link.D_in_ppe = Nspans_tot*ft.length*1e-3*ft.disp+PPEparams.phys.pd;
        end
    else
        PPEparams.link.method_BP    = "pd";
        PPEparams.link.D_in_ppe     = Nspans_tot*ft.length*1e-3*ft.disp+PPEparams.phys.pd;
    end

    PPEparams.link.nsteps_BP        = floor(PPEparams.link.D_in_ppe/ft.disp*ratio_length);

    %%% statistics
    % seed values 
    %   - 1  if we always want the same symbol flow
    %   - 0  if we want different random symbol flows
    
    if isfield(PPEparams,'repet') == 0
        PPEparams.repet.Ntries = 10;
    end

    if PPEparams.repet.Ntries == 1
        PPEparams.repet.seed        = 1;
        PPEparams.repet.plot        = 1;
    else
        PPEparams.repet.seed        = 0;
        PPEparams.repet.plot        = 0;
    end

    if PPEparams.repet.Ntries > 10
        PPEparams.repet.parfor      = 1;
        PPEparams.repet.plot        = 0;
    else
        PPEparams.repet.parfor      = 0;
    end

    %% POWER PROFILE STUDY
    if isfield(PPEparams,'findloc') == 0
        PPEparams.findloc.av_method     = 'mirror';
        PPEparams.findloc.av_period     = 10;
        PPEparams.findloc.uncertainty   = 5e3;        % [m]
    end

    %% REFERENCE POWER PROFILE
    if isfield(PPEparams.plot,"ref") == 0
        PPEparams.plot.ref              = struct();
    end
    if isfield(PPEparams.plot.ref,'what') == 0
        PPEparams.plot.ref.what         = 'heavy';
    end
   
    if isfield(PPEparams.plot.ref,'std') == 0
        PPEparams.plot.ref.std          = 0;
    end

    %% ANOMALIES
    if isfield(PPEparams,'anomalies') == 0
        PPEparams.anomalies.ai = struct();
    end

    if isfield(PPEparams.anomalies,'get_ai') == 0
        PPEparams.anomalies.get_ai = 0;
        PPEparams.anomalies = rmfield(PPEparams.anomalies,'ai');
    else
%         if PPEparams.anomalies.get_ai == 1
            if isfield(PPEparams.anomalies,'method') == 0
                PPEparams.anomalies.method      = 'pd';
            end
        
            if tx.pd == 0
                PPEparams.anomalies.method      = 'correction';
                if isfield(PPEparams.anomalies.ai,'correction') == 0 || ...
                    PPEparams.anomalies.ai.correction.getAFC == 0
        
                    PPEparams.anomalies.ai.correction.getAFC= 0;
                    Rs  = varargin{2}.symbrate;
                    PPEparams.anomalies.ai.correction.AFC   = 0.173*(Rs^(-1.79));% OECC
                else
                    PPEparams.anomalies.ai.correction.getAFC= 1;            
                end
            else
                PPEparams.anomalies.method          = 'pd';
            end
        
            if isfield(PPEparams.anomalies,'ai') == 0
                PPEparams.anomalies.ai                  = struct();
            end
        
            % resolution
            if isfield(PPEparams.anomalies.ai,'resolution') == 0
                PPEparams.anomalies.ai.resolution.round = get_decimals(PPEparams.link.dl/1000);
                PPEparams.anomalies.ai.resolution.zloss = 20;   % [km]
            end
        
            % calibration
            if isfield(PPEparams.anomalies.ai,'calibration') == 0
                PPEparams.anomalies.ai.calibration.getCal       = 1;
                PPEparams.anomalies.ai.calibration.lossdB       = 3;
                PPEparams.anomalies.ai.calibration.PrefAmpLin   = tx.Plin;
                PPEparams.anomalies.ai.calibration.PrefAmpdBm   = tx.PdBm;
            end
        
            if isfield(PPEparams.anomalies.ai.calibration,'calC') == 1
                PPEparams.anomalies.ai.calibration.getCal       = 0;
            end
        
            if PPEparams.anomalies.ai.calibration.getCal == 0
                assert(isfield(PPEparams.anomalies.ai.calibration,'calC') == 1,...
                    ' == SET_PPEparams ==\n Calibration factor is missing')
            end
            if isfield(PPEparams.anomalies.ai,'PrefAmpLin') == 0
                PPEparams.anomalies.ai.calibration.PrefAmpLin   = tx.Plin;
                PPEparams.anomalies.ai.calibration.PrefAmpdBm   = tx.PdBm;
            end
        
            % apodisation
            % In case the AI is too noisy, it helps preventing
            % localisation if we use the maximum of the derivative
        
            if isfield(PPEparams.anomalies.ai,'apodisation') == 0
                PPEparams.anomalies.ai.apodisation      = struct();
            end
        
            if isfield(PPEparams.anomalies.ai.apodisation,'helpApo') == 0
                PPEparams.anomalies.ai.apodisation.helpApo = 0;
            end
        
            if PPEparams.anomalies.ai.apodisation.helpApo == 1
                PPEparams.anomalies.ai.apodisation.makeApo   = 1;
            end
                
            if isfield(PPEparams.anomalies.ai.apodisation,'makeApo') == 0
                PPEparams.anomalies.ai.apodisation.makeApo   = 0;
            end
        
            if PPEparams.anomalies.ai.apodisation.helpApo == 1        
                lossy_span = 0;
                for k = 1:length(fieldnames(amp.losses))-1
                    tmp_help = amp.losses.(sprintf('span%i',k));
                    if tmp_help(2)~= 0
                        lossy_span = k;
                        break
                    end
                end
                if lossy_span == 0
                    lossy_span = [];
                end
                assert(isempty(lossy_span) ~= 1, " == SET_PPEparams ==\n There is no loss inserted ")
                PPEparams.anomalies.ai.apodisation.zloss = (lossy_span-1)*ft.length*1e-3+tmp_help(1)*1e-3;
            end
        
            % estimation
            if isfield(PPEparams.anomalies.ai,'estimation') == 0
                PPEparams.anomalies.ai.estimation = struct();
                PPEparams.anomalies.ai.estimation.zloss.method = "maxAI";
            end
        
            if isfield(PPEparams.anomalies.ai.estimation,'method') == 0
                PPEparams.anomalies.ai.estimation.zloss.method = "maxAI";
            end
        
            if isfield(PPEparams.anomalies.ai.estimation,'surimposition') == 0
                PPEparams.anomalies.ai.estimation.surimposition    = "mean";
            end
%         end
    end
    
    %% SORTING THE STRUCTURES
    PPEparams   = sort_struct_alphabet(PPEparams);  %

%     if PPEparams.anomalies.get_ai == 1
        PPEparams.anomalies = sort_struct_alphabet(PPEparams.anomalies);% ---
        PPEparams.anomalies.ai  = sort_struct_alphabet(PPEparams.anomalies.ai); % --- ---
        PPEparams.anomalies.ai.apodisation  = sort_struct_alphabet(PPEparams.anomalies.ai.apodisation); % --- --- ---
        PPEparams.anomalies.ai.calibration  = sort_struct_alphabet(PPEparams.anomalies.ai.calibration); % --- --- ---
        PPEparams.anomalies.ai.resolution   = sort_struct_alphabet(PPEparams.anomalies.ai.resolution);  % --- --- ---
%     else
%         PPEparams = rmfield(PPEparams,"anomalies");
%     end

    PPEparams.findloc   = sort_struct_alphabet(PPEparams.findloc);  % ---
    PPEparams.link      = sort_struct_alphabet(PPEparams.link);     % ---
    PPEparams.method    = sort_struct_alphabet(PPEparams.method);   % ---
    PPEparams.method.pp     = sort_struct_alphabet(PPEparams.method.pp);    % ---
    PPEparams.phys      = sort_struct_alphabet(PPEparams.phys);     % ---
    PPEparams.plot      = sort_struct_alphabet(PPEparams.plot);     % ---
    PPEparams.plot.ref      = sort_struct_alphabet(PPEparams.plot.ref);     % ---
    PPEparams.repet     = sort_struct_alphabet(PPEparams.repet);    % ---




