function pd_needed = set_pd(confidence_level,symbol_rate)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : SET_PD - SET PRE DISPERSION
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-10-11
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Compute the needed pre-disperson to get flat profile
%
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------

    S = sum(confidence_level==[90,95,99]);
    assert(S~=0,"Confidence level can only be 90,95,or 99")
    
    if confidence_level == 90
        coeff = 1.25;
    elseif confidence_level == 95
        coeff = 0.625;
    elseif confidence_level == 99
        coeff = 0.725;
    else
        return
    end
    
    a = 1.32e7*(confidence_level/100).^34.2;
    b = -3.36*(confidence_level/100)+1.17;
    
    pd_needed = 20*ceil(coeff*a*symbol_rate.^b);

end