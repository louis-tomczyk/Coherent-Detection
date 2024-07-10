function bool = is_integer(in)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : IS_INTEGER
%   Author          : Bruno LUONG
%   Institution     : MathWorks help
%   Email           : 
%   Date            : 2022-09-16
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Check if the input is an integer
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    bool = mod(double(in),1) == 0;

end
