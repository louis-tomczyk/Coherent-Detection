function apo = my_gausswin(N,center,res,varargin)
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : MY_GAUSSWIN
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

    if nargin == 3
        Order = 3;
    else
        Order = varargin{1};
    end

    x     = 1:N;
    sigma = res/((8^Order*log(2)).^(1/2/Order));
    apo   = exp(-((x-center)./(sqrt(2)*sigma)).^(2*Order));
    apo   = apo.';