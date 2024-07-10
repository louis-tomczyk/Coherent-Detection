function [xdiff,ydiff] = get_derivative(x,y,ordre)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_DERIVATIVE
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-01-28
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Compute recurssively the derivatives of a signal
%
% ----- INPUTS -----
%   X       (array)     x-axis
%   Y       (array)     y-axis
%   ORDRE   (integer)   derivative order
%
% ----- OUTPUTS -----
%   XDIFF   (array)     x-axis
%   YDIFF   (array)     y-axis
%
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------

    dx      = x(2)-x(1);
    xdiff   = x(1:end-1);
    xdiff   = xdiff + dx/2;
    ydiff   = diff(y)/dx;

    if ordre > 1
        for k = 2:ordre
            [xdiff,ydiff] = get_derivative(xdiff,ydiff,1);
        end
    end