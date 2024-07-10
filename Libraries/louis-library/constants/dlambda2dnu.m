function dnu = dlambda2dnu(dlambda,lambda)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : DLAMBDA2DNU
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-05-10
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Convert the linewidth DLAMBDA expressed in wavelenth
%   into the linewidth DNU expressed in frequency.
%
% ----- INPUTS -----
%   DLAMBDA [nm]
%   LAMBDA  [nm]
%
% ----- OUTPUTS -----
%   DNU     [Hz]
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

  lambda  = lambda*1e-6;            % [m]
  dlambda = dlambda*1e-6;           % [m]

  c       = 299792458;              % [m/s]
  dnu     = c.*dlambda/(lambda.^2); % [Hz]
end
