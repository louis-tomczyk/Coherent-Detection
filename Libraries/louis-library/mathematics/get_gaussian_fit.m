function [fitresult, gof,varargout] = get_gaussian_fit(x,y,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_GAUSSIAN_FIT
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2022-09-10
%   Version         : 1.1
%
% ----- Main idea -----
%   Create, export and plot the gaussian fit from data
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
%                           the fit: fitresult(x)= a1*exp(-((x-b1)/c1)^2)
%                           that can be accessed via structure call:
%                           fitresult.a1 - fitresult.b1 - fitresult.c1
%   GOF         (structure) contains:
%                           SSE     (scalar) Summed Square of Residuals
%                           RSQUARE (scalar) Determination coefficient
%                           RMSE    (scalar) Root Mean Squared Error
% ----- BIBLIOGRAPHY -----
%   Functions           : Curve Fitter App
%   Author              : Matlab
%   Author contact      : 
%   Date                : 2022-09-30
%   Title of program    : createFit
%   Code version        : 1.0
%   Type                : 
%   Web Address         : 
% ---------------------------------------------

% reshaping to columns
[xData, yData]  = prepareCurveData( x, y );

% Set up fittype and options.
ft              = fittype( 'gauss1' );
opts            = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display    = 'Off';
opts.Lower      = [-Inf -Inf 0];

% Fit model to data.
[fitresult,gof] = fit( xData, yData, ft, opts );

% export and plot
if nargin >= 3
    for k = 3:nargin
        if strcmp(varargin{k-2},'data') == 1
            varargout{1} = feval(fitresult,x);
        elseif strcmp(varargin{k-2},'plot on') == 1
            figure( 'Name', 'untitled fit 1' );
                h = plot( fitresult, xData, yData );
                legend( h, 'y vs. x', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
                xlabel( 'x', 'Interpreter', 'none' );
                ylabel( 'y', 'Interpreter', 'none' );
                grid on
        end
    end
end
