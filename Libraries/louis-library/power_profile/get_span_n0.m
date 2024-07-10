function nhat = get_span_n0(z,PPEparams)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name      GET_SPAN_N0
%   Author             louis tomczyk
%   Institution        Telecom Paris
%   Email              louis.tomczyk@telecom-paris.fr
%   Date               2023-02-16
%   Version            2.1
%
% ----- Main idea -----
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    Lspan   = PPEparams.link.dl*PPEparams.link.nsteps_span*1e-3;
    
    nhat    = floor(z/Lspan+1/2);
    zmin    = nhat*Lspan-PPEparams.anomalies.ai.resolution.zloss;
    zmax    = (nhat+1)*Lspan-PPEparams.anomalies.ai.resolution.zloss;

    if z >= zmin && z<zmax
%         fprintf(" == GET_SPAN_N0 ==\n loss location close to the next span\n")
        return
    else
        nhat = floor(z/Lspan);
    end