function PP_out = get_pp(PPEparams,Axis,las,tx,ft,amp)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_PP - POWER PROFILE
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.com
%   Date            : 2023-01-27
%   Version         : 1.0
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------

     assert(strcmp(PPEparams.link.method_BP,'pd') == 1)

    PP                  = set_ppe_init(PPEparams);
    [PP{:}]             = get_pp_all(PPEparams,Axis,las,tx,ft,amp);

    if PPEparams.repet.Ntries > 1
        [~,~,PP_out]    = set_ppe_output(2,'pd',PP);
    else
        PP_out          = set_ppe_output(1,'pd',PP);
    end
end