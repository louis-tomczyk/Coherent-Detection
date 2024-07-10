function Ticklabels = ticks2ticklabels(ticks,decimals)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : TICKS2TICKLABELS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.-tomczyk.work@gmail.com
%   Date            : 2022-08-01
%   Version         : 1.0
%
% ----- Main idea -----
%
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    Ticklabels = string(ticks);
    for k = 1:length(Ticklabels)
        tmp = char(Ticklabels(k));
        index_comma = strfind(tmp,'.');
        if isempty(index_comma) == 1
            Ticklabels(k) = tmp;
        else
            if decimals == 0
                Ticklabels(k) = tmp(1:index_comma-1);
            else
                Ticklabels(k) = tmp(1:index_comma +decimals);
            end
        end
    end

    Ticklabels = cellstr(Ticklabels);

end