function ppc = get_ppc(pp_original,FC)
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_PPC - Power Profile Corrected
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

    if FC(1) == Inf
        FC  = FC(2:end);
    end

    Npp     = length(pp_original);
    Nfc     = length(FC);
    NNN     = min(Npp,Nfc);

    if Npp ~= Nfc
        if NNN == Npp
            ppc     = FC(2:end).*pp_original;
            offset  = mean(pp_original)-mean(ppc);
        elseif NNN == Nfc
            ppc     = FC.*pp_original(Npp-NNN+1:end);
            offset  = mean(pp_original(Npp-NNN+1:end))-mean(ppc);
        end
    end

    ppc     = ppc+offset;