function [fields,powers] = field2symbs(E,Axis)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : FIELD2SYMBS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.-tomczyk.work@gmail.com
%   Date            : 2022-08-03
%   Version         : 1.0
%
% ----- Main idea -----
%   Get the the symbols from the input field
%
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    if isstruct(E) == 1
        E   = E.field;
    end

    fields  = zeros(Axis.Nsymb,Axis.Nsps);

    for k = 1:Axis.Nsymb
        fields(k,:) = E((k-1)*Axis.Nsps+1:k*Axis.Nsps);
    end

    powers = empower(fields).';

end