function [AFC,FC] = get_FC(PPEparams,pp_ref,Axis,las,tx,ft,amp)

    % ---------------------------------------------
    % ----- INFORMATIONS -----
    %   Function name   : GET_FC
    %   Author          : louis tomczyk
    %   Institution     : Telecom Paris
    %   Email           : louis.tomczyk@telecom-paris.fr
    %   Date            : 2023-03-02
    %   Version         : 1.1
    %
    % ----- Main idea -----
    % ----- INPUTS -----
    % ----- OUTPUTS -----
    % ----- BIBLIOGRAPHY -----
    % ---------------------------------------------

    assert(PPEparams.plot.standardise == 0,...
     " == GET_FC ==\n" + ...
     "Correction not working if standardized original power profile")
    assert(PPEparams.plot.ref.std == 0,...
     " == GET_FC ==\n" + ...
     "Correction not working if standardized reference power profile")

    if PPEparams.anomalies.ai.correction.getAFC == 1
        AFC = get_AFC(PPEparams,pp_ref,Axis,las,tx,ft,amp);
    else
        AFC = PPEparams.anomalies.ai.correction.AFC;
    end

    z       = PPEparams.plot.dist;
    FC      = z.^(-AFC);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------
% ----- CONTENTS -----
%   get_AFC
%   get_power_fit
% ---------------------------------------------

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_CORR_PARAM
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-08-29
%   Version         : 1.0.2
%
% ----- Main idea -----
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

function AFC = get_AFC(PPEparams,pp_ref,Axis,las,tx,ft,amp)

    tx_flat             = tx;
    tx_flat.pd          = set_pd(90,Axis.symbrate);
    PPEparams_flat      = set_PPEparams(tx_flat,ft,amp,PPEparams);
    pp_flat             = get_pp(PPEparams_flat,Axis,las,tx_flat,ft,amp);

    z                   = PPEparams.plot.dist(2:end);
    ratio               = pp_flat./pp_ref(2:end);
    ratio               = ratio(2:end);

    [fitresult, ~]      = get_power_fit(z, ratio);
    AFC                 = -fitresult.b;
    assert(AFC>0, ...
        sprintf(" ==--- GET_AFC ---==\n" + ...
        "Fit failed, AFC= %.1e",AFC))

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_POWER_FIT
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-02-01
%   Version         : 1.0
%
% ----- Main idea -----
%   Create, export and plot the power fit from data
%
% ----- INPUTS -----
%   X           (array)     The X-axis of the data
%   Y           (array)     The Y-axis of the data
%   varagin{1}  (string){optional}
%                           if "DATA", then the Y-values of the fit are
%                           exported
%                           if "PLOT ON", then the fit is plotted
%  varargin{2}  (string){optional}
%                           Same as varargin{1}
% ----- OUTPUTS -----
%   FITRESULT   (CFIT)      contains: 
%                           - FITRESULT (formal function) whose values can
%                           be accessed with the FEVAL function or just by
%                           putting the x-values as argument
%                           - (structure) containing the coefficients of
%                           the fit: fitresult(x)= a*x^b
%                           that can be accessed via structure call:
%                           fitresult.a - fitresult.b
%   GOF         (structure) contains:
%                           SSE     (scalar) Summed Square of Residuals
%                           RSQUARE (scalar) Determination coefficient
%                           RMSE    (scalar) Root Mean Squared Error
% ----- BIBLIOGRAPHY -----
%   Functions           : Curve Fitter App
%   Author              : Matlab
%   Author contact      : 
%   Date                : 2023-02-01
%   Title of program    : createFit
%   Code version        : 1.0
%   Type                : 
%   Web Address         : 
% ---------------------------------------------

function [fitresult, gof,varargout] = get_power_fit(z, AFC,varargin)

% Fit: 'untitled fit 1'.
[xData, yData]  = prepareCurveData( z*1e3, AFC );

% Set up fittype and options.
ft              = fittype( 'power1' );
opts            = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display    = 'Off';
opts.StartPoint = [0.953068610647315 -0.000495760929578534];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% export and plot
if nargin >= 3
    for k = 3:nargin
        if strcmp(varargin{k-2},'data') == 1
            varargout{1} = feval(fitresult,x);
        elseif strcmp(varargin{k-2},'plot on') == 1
            figure( 'Name', 'power fit' );
                h = plot( fitresult, xData, yData );
                legend( h, 'ratio vs. z', 'power fit','Location', 'NorthEast');
                xlabel( 'x', 'Interpreter', 'none' );
                ylabel( 'y', 'Interpreter', 'none' );
                grid on
        end
    end
end



    
