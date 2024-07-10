function E = symbs2field(symbols)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : SYMBS2FIELD
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.-tomczyk.work@gmail.com
%   Date            : 2023-01-12
%   Version         : 1.0.1
%
% ----- Main idea -----
%   Get the the field from the symbols
%
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    [Nsymb,Nsps]= size(symbols);
    E           = zeros(1,Nsymb.*Nsps);

    for k = 1:Nsymb
        E((k-1)*Nsps+1:k*Nsps) = symbols(k,:);
    end

end