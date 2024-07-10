function [Ex,Ey] = sep_XYfields(E)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : SEP_XYFIELDS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-08-05
%   Version         : 1.0
%
% ----- Main idea -----
%   Replace the PBS function of OPTILUX
%
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    if size(E.field,2) ~= 2*length(E.lambda)
        error('A beam splitter can be used only in dual-polarization mode.')
    end
    
    Ex.lambda   = E.lambda;
    Ey.lambda   = E.lambda;
    
    Ex.field    = E.field(:,1);
    Ey.field    = E.field(:,2);

end

