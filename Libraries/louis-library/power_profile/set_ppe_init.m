function PP = set_ppe_init(PPEparams)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : SET PPE INIT
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-02-24
%   Version         : 1.1
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    if PPEparams.repet.Ntries > 1
        % --- back propagation only in the link
        if strcmp(PPEparams.link.method_BP,'no_pd') == 1
            PP  = cell(1,3);
        else % --- back propagation in link + predisp
            PP  = cell(1,4);
        end
    else % Ntries = 1
        % --- back propagation only in the link
        if strcmp(PPEparams.link.method_BP,'no_pd') == 1
            PP  = cell(1,2);
        else % --- back propagation in link + predisp
            PP  = cell(1,3);
        end
    end