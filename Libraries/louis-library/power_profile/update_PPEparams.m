function PPEparams = update_PPEparams(PPEparams,ft,amp,what_charac,iter,epsilon,pds)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : UPDATE_PPEPARAMS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.-tomczyk.work@gmail.com
%   Date            : 2022-12-31
%   Version         : 1.0
%
% ----- Main idea -----
%   Set and update the Power Profile Estimator parameters
%
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

if strcmp(what_charac,'eps') == 1
    fprintf("\titer % = %i\n",what_charac,epsilon(iter))
    PPEparams.phys.nl_factor    = epsilon(iter);
    PPEparams.phys.pd           = pds;
else
    fprintf("\titer % = %i\n",what_charac,pds(iter))
    PPEparams.phys.nl_factor    = epsilon;
    PPEparams.phys.pd           = pds(iter);
end

PPEparams   = set_PPEparams(PPEparams,ft,amp);

    