function [PPEparams,ppref,pps,AI]  = get_anomalies_alix(PPEparams,ppref,pploss,amp,ft,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_ANOMALIES
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-02-08
%   Version         : 1.1
%
% ----- MAIN IDEA -----
%   Estimate the losses from power profile estimations
%
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------

    %%% correction if needed
    if nargin > 5
        Axis                = varargin{1};
        las                 = varargin{2};
        tx                  = varargin{3};

        [~,pprefc,pplossc,PPEparams]= correct_pp(PPEparams,pploss,ppref,amp,ft,Axis,las,tx);
    else
        [~,pprefc,pplossc]    = correct_pp(PPEparams,pploss,ppref,amp);
    end

    %%% Profiles scaling for surimposition
    pps         = get_pps(PPEparams,pplossc,pprefc);

    %%% Anomaly indicator
    AI          = get_AI(PPEparams,pprefc,pps);

    %%% Get loss localisation
    PPEparams   = get_loss_lochat(PPEparams,AI);

    %%% Calibration if needed
    PPEparams   = get_calibration(PPEparams,ft);

    %%% Estimate the loss
    PPEparams   = get_loss_valhat(PPEparams,ft);

    %%% Sort and export PPEparams
    PPEparams   = sort_struct_alphabet(PPEparams);

    zloss       = PPEparams.anomalies.ai.estimation.zloss.value;
    Apeak       = PPEparams.anomalies.ai.estimation.Apeak;
    calF        = PPEparams.anomalies.ai.calibration.calC*PPEparams.anomalies.ai.calibration.PrefAmpLin;

    export_name_base = sprintf('PPEparams - zloss %i km - Apeak %.4e - calF %.4e', ...
                               zloss,Apeak,calF);

    if strcmp(PPEparams.anomalies.method,'correction') == 1
        AFC     = PPEparams.anomalies.ai.correction.AFC;
        export_name_method = sprintf(' - AFC %.4e',AFC);
    elseif strcmp(PPEparams.anomalies.method,'pd') == 1
        pd      = PPEparams.phys.pd;
        export_name_method = sprintf(' - pd %i ps-nm',pd);
    end

    if PPEparams.anomalies.ai.calibration.getCal == 0
        lossdB  = PPEparams.anomalies.ai.estimation.lossdB;
        if isreal(lossdB)
            export_name_cal = sprintf(' - %.2f dB',lossdB);
        else
            export_name_cal = sprintf(' - %.2f dB(C)',lossdB);
        end
    else
        export_name_cal = '';
    end

    export_name = strcat(export_name_base,export_name_method,export_name_cal);
    export_structure(PPEparams,export_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------
% ----- CONTENTS -----
%   correct_pp
%   get_pps
%   get_AI
%   get_loss_lochat
%   get_calibration
%   get_loss_valhat
% ---------------------------------------------

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : CORRECT_PP
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-02-08
%   Version         : 1.1
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------
function [z1,ppref,pploss,PPEparams] = correct_pp(PPEparams,pploss,ppref,amp,ft,varargin)

    if strcmp(PPEparams.anomalies.method,'correction') == 1
        
        assert(nargin == 8," ==--- CORRECT_PP ---==\n AXIS and/or LAS and/or structure(s) are/is missing")
        assert(isstruct(varargin{1}) == 1 && isfield(varargin{1},'symbrate') == 1,...
            ' ==--- CORRECT_PP ---==\n varargin{1} should be the AXIS structure')
        assert(isstruct(varargin{2}) == 1 && isfield(varargin{2},'lam') == 1,...
            ' ==--- CORRECT_PP ---==\n varargin{2} should be the LAS structure')
        assert(isstruct(varargin{3}) == 1 && isfield(varargin{3},'modfor') == 1,...
            ' ==--- CORRECT_PP ---==\n varargin{3} should be the TX structure')

        Axis    = varargin{1};
        las     = varargin{2};
        tx      = varargin{3};

        % z-axis
        z1      = PPEparams.plot.dist(2:end);

        % dealing with corrected profile if there is as the correction gives a profile with 
        % (nsteps_fibre-1) points

        dyloss  = pploss(2)-pploss(1);
        dyref   = ppref(2)-ppref(1);
 
        pploss  = [pploss(1)+dyloss,pploss];
        ppref   = [ppref(1)+dyref,ppref];

        % correction function
        [AFC,FC]= get_FC(PPEparams,ppref,Axis,las,tx,ft,amp);
        PPEparams.anomalies.ai.correction.getAFC    = 1;
        PPEparams.anomalies.ai.correction.AFC       = AFC;

        % profile corrections
        ppref   = get_ppc(ppref,FC);
        pploss  = get_ppc(pploss,FC);
    
    else
        PPEparams.anomalies.ai.correction.getAFC    = 0;
        z1 = PPEparams.plot.dist;
    end

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_PPS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-01-31
%   Version         : 1.0
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------
function pps = get_pps(PPEparams,pploss,ppref)

    % method to get the offset
    % - for low number of realisations (<100), best is the "mean" method
    % - see ``MSE and DPAPR1.ods'' file

    surimp = PPEparams.anomalies.ai.estimation.surimposition;
    if strcmp(surimp,"start")==1
        tmp_loss    = pploss(1);
        tmp_ref     = ppref(1);
    elseif strcmp(surimp,"end")==1
        tmp_loss    = pploss(end);
        tmp_ref     = ppref(end);
    elseif strcmp(surimp,"mean")==1
        tmp_loss    = mean(pploss);
        tmp_ref     = mean(ppref);
    end

    % computing the needed offset
    diff_tmp        = tmp_ref-tmp_loss;
    offset          = abs(diff_tmp);

    % applying the offset
    if diff_tmp>0
        pps         = pploss+offset;
    else
        pps         = pploss-offset;
    end

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_AI
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-01-31
%   Version         : 1.0
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------
function AI = get_AI(PPEparams,ppref,pps,varargin)

    assert(nargin >= 3,' ==--- GET_AI ---==\n resolution cannot be set','err');

    % get the AI
    AI      = abs(ppref-pps);

    NAI     = length(AI);   
    papr    = get_pdapr('papr',AI);

    if papr<5
        PPEparams.anomalies.ai.apodisation.makeApo = 1;
        sprintf(" ==--- GET_AI ---==\n Apodisation performed")
    end

    if PPEparams.anomalies.ai.apodisation.makeApo == 1
        % setting the width of the window
        if nargin == 3
            winWidth = floor(3*PPEparams.link.nsteps_span/2);
        else
            winWidth = varargin{1};
        end
        
        % computing the window
        if PPEparams.anomalies.ai.apodisation.makeApo == 1
            if PPEparams.anomalies.ai.apodisation.helpApo == 0
                [~,Xm]  = max(AI);
                my_apo  = my_gausswin(NAI,Xm,winWidth);
            else
                sprintf(" ==--- GET_AI ---==\n Apodisation performed with - helpApo = %i", ...
                    PPEparams.anomalies.ai.apodisation.helpApo)
                 my_apo  = my_gausswin(NAI,PPEparams.anomalies.ai.apodisation.zloss,winWidth);
            end
        end
        % reshaping to match the sizes
        [nr,nc] = size(AI);
        my_apo  = reshape(my_apo,[nr,nc]);
    
        % applying the window
        AI      = AI.*my_apo;
    end

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_LOSS_LOCHAT
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-02-06
%   Version         : 1.0
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------
function PPEparams = get_loss_lochat(PPEparams,AI)

    if strcmp(PPEparams.anomalies.ai.estimation.zloss.method,"maxAI") == 1
        [Apeak,zhat]    = max(AI);
    elseif strcmp(PPEparams.anomalies.ai.estimation.zloss.method,"maxDiffAI") == 1
        [Apeak,zhat]    = max(diff(AI));
    end

    if isempty(PPEparams.anomalies.ai.estimation.zloss.helpLocVal) == 1
        zhat    = round(PPEparams.plot.dist(zhat),PPEparams.anomalies.ai.resolution.round.n); % [km]
    else
        zhat    = PPEparams.anomalies.ai.estimation.zloss.helpLocVal;
    end
    
    PPEparams.anomalies.ai.estimation.Apeak         = Apeak;
    PPEparams.anomalies.ai.estimation.zloss.value   = zhat;

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_CALIBRATION
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-02-17
%   Version         : 1.2
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% PTL, May
% ----------------------------------------------
function PPEparams = get_calibration(PPEparams,ft)

    if PPEparams.anomalies.ai.calibration.getCal == 1

        zloss       = PPEparams.anomalies.ai.calibration.zloss;             % [km]
        Apeak       = PPEparams.anomalies.ai.esimation.Apeak;

        span_nb     = get_span_n0(zloss,PPEparams);
        Prefamp     = PPEparams.anomalies.ai.calibration.PrefAmpLin;        % [mW]
        lossdB      = PPEparams.anomalies.ai.calibration.lossdB;            % [dB]

        T0          = 10^(-lossdB/10);                                      % []
        lossLin     = 1-T0;                                                 % []
        FibAttLin   = 10^(-ft.alphadB*(zloss-span_nb*ft.length*1e-3)/10);   % []
        C           = Apeak/(lossLin*Prefamp*FibAttLin);                    % [1/W]
   
        PPEparams.anomalies.ai.calibration.AIpeak   = Apeak;
        PPEparams.anomalies.ai.calibration.calC     = C;
        PPEparams.anomalies.ai.calibration.FibAttdB = -10*log10(FibAttLin);
        PPEparams.anomalies.ai.calibration.FibAttLin= FibAttLin;
        PPEparams.anomalies.ai.calibration.T0       = T0;
        PPEparams.anomalies.ai.calibration.lossLin  = lossLin;

        PPEparams.anomalies.ai.calibration = sort_struct_alphabet(PPEparams.anomalies.ai.calibration);
    else
        return
    end

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name      GET_LOSS_VALHAT
%   Author             louis tomczyk
%   Institution        Telecom Paris
%   Email              louis.tomczyk@telecom-paris.fr
%   Date               2023-02-13
%   Version            2.1
%
% ----- Main idea -----
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------
function PPEparams = get_loss_valhat(PPEparams,ft)

    if PPEparams.anomalies.ai.calibration.getCal == 0
        zhat        = PPEparams.anomalies.ai.estimation.zloss.value;
        Apeak       = PPEparams.anomalies.ai.estimation.Apeak;
        
        span_nb     = get_span_n0(zhat,PPEparams);
        CalF        = PPEparams.anomalies.ai.calibration.calF;
        FibAttdB    = ft.alphadB*(zhat-span_nb*ft.length*1e-3);
        FibAttLin   = 10^(-FibAttdB/10);
        
        lossHatLin  = Apeak/(CalF*FibAttLin);
        T0          = 1-lossHatLin;
        if T0<0
            sprintf([' ==--- GET_LOSS_VALHAT ---==\n ' ...
                'LossHat.lin should be < 1, Apeak = %.2e - CalF*FibAttLin = .2e'], ...
                Apeak,CalF*FibAttLin)
        end

        PPEparams.anomalies.ai.estimation.losslin = lossHatLin;
        PPEparams.anomalies.ai.estimation.lossdB  = -10*log10(T0);
    else
        fprintf(" ==--- GET_LOSS_VALHAT ---==\n No loss estimation as Calibration step : ON\n\n")
        return
    end

