function plot_params = set_plot_params(PPEparams,ft,amp,Axis,tx,varargin)

%offset,what_charac, epsilon,ft,amp,Axis,tx,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : SET_PLOT_PARAMS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-10-13
%   Version         : 1.2.1
%
% ----- MAIN IDEA -----
%  Set the title for the files, specifically for Power Profiles Estimations
%
% ----- INPUTS -----
%   PPEPARAMS   (structure)
%   ITER        (scalar)
%   OFFSET      (scalar)
%   WHAT_CHARAC (string)
%   PDS         (scalar)[ps/nm]
%   EPSILON     (scalar)[]
%   FT          (structure)     See SET_FT function for details
%   AMP         (structure)     See SET_AMP function for details
%   AXIS        (structure)     See SET_AXIS function for details
%   TX          (structure)     See ET_AXIS function for details
%
% ----- OUTPUTS -----
%   PLOT_PARAMS (structure)     Contains: 
%                               - NORM
%                               - TITLE
%                               - ALPHA
%                               - WHAT_CHARAC
%                               - FIBRE
%                               - TX_AMP
%                               - PPE
%                               - SUBTITLE
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------

    input_names     = strings(1,nargin);
    for k = 1:nargin
        input_names(k) = string(inputname(k));
    end    

    if sum('iter' == input_names) == 1,iter = varargin{1}; end
    if sum('offset' == input_names) == 1,offset = varargin{2}; end
    if sum('what_charac' == input_names) == 1,what_charac = varargin{3}; end
    if sum('epsilon' == input_names) == 1,epsilon = varargin{4}; end

    %% mandatory arguments
    plot_params.norm        = PPEparams.plot.norm;
    plot_params.title       = 'power profile estimation';
    plot_params.alpha       = strcat('alpha_',num2str(ft.alphadB));
    
    plot_params.fibre       = sprintf("disp %.1f - gf %.2f - att %.2f --- ", ...
                                      ft.disp,ft.gf,ft.alphadB);
    plot_params.tx_amp      = sprintf("Rs %i - Nsymb %i - %i dBm - %i spans --- ", ...
                                      Axis.symbrate,log2(Axis.Nsymb),tx.PdBm,amp.Nspan);

    %% optional arguments
    if sum('what_charac' == input_names) == 0
        what_charac = "pd";
        epsilon     = PPEparams.phys.nl_factor;
        pds         = PPEparams.phys.pd;
    end

    plot_params.what_charac = what_charac;


    if sum('iter' == input_names) == 0 && sum('what_charac' == input_names) == 0
            plot_params.ppe = sprintf("pd %i - eps %.1e - nrea %i --- ", ...
                                       pds,epsilon,PPEparams.repet.Ntries);
    else
        if strcmp(plot_params.what_charac,'epsilon') == 1
            plot_params.ppe = sprintf("%i --- pd %i - eps %.1e - nrea %i --- ", ...
                                       iter+offset,PPEparams.phys.pd,epsilon(iter), ...
                                       PPEparams.repet.Ntries);
        else
            plot_params.ppe = sprintf("%i --- pd %i - eps %.1e - nrea %i --- ",...
                                       iter+offset,PPEparams.phys.pd,epsilon, ...
                                       PPEparams.repet.Ntries);
        end
    end


    %% merging
    plot_params.subtitle    = strcat(plot_params.ppe,plot_params.fibre, ...
                                plot_params.tx_amp,string(datetime));









